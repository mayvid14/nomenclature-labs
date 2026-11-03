//! 2D coordinates for skeletal formulas: zig-zag main chain or regular-polygon
//! ring, then branches placed into the widest free angle.

use crate::molecule::Molecule;
use crate::namer::molecule_rings;
use std::collections::VecDeque;
use std::f32::consts::PI;

fn deg(d: f32) -> f32 {
    d * PI / 180.0
}

fn is_linear(m: &Molecule, a: usize) -> bool {
    let n = m.neighbors(a);
    n.iter().any(|x| x.1 == 3) || n.iter().filter(|x| x.1 == 2).count() >= 2
}

fn norm(a: f32) -> f32 {
    let mut a = a % (2.0 * PI);
    if a < 0.0 {
        a += 2.0 * PI;
    }
    a
}

/// Lay out the molecule; `backbone` is the main chain or ring in locant order.
pub fn layout(m: &mut Molecule, backbone: &[usize], ring: bool) {
    let n = m.atoms.len();
    if n == 0 {
        return;
    }
    let rings = molecule_rings(m);
    let ring_of = |a: usize| rings.iter().position(|r| r.contains(&a));
    let mut placed = vec![false; n];
    let mut queue: VecDeque<usize> = VecDeque::new();
    let mut placed_rings = vec![false; rings.len()];

    let set = |m: &mut Molecule, a: usize, x: f32, y: f32, placed: &mut Vec<bool>| {
        m.atoms[a].x = x;
        m.atoms[a].y = y;
        placed[a] = true;
    };

    let bb: Vec<usize> = if backbone.is_empty() { vec![0] } else { backbone.to_vec() };
    if ring && bb.len() >= 3 {
        let k = bb.len();
        let r = 0.5 / (PI / k as f32).sin();
        for (i, &a) in bb.iter().enumerate() {
            let ang = deg(-90.0) + 2.0 * PI * i as f32 / k as f32;
            set(m, a, r * ang.cos(), r * ang.sin(), &mut placed);
            queue.push_back(a);
        }
        if let Some(ri) = ring_of(bb[0]) {
            placed_rings[ri] = true;
        }
    } else {
        let mut x = 0.0;
        let mut y = 0.0;
        let mut ang = deg(-30.0);
        for (i, &a) in bb.iter().enumerate() {
            if i > 0 {
                if i >= 2 && !is_linear(m, bb[i - 1]) {
                    ang = if ang < 0.0 { deg(30.0) } else { deg(-30.0) };
                }
                x += ang.cos();
                y += ang.sin();
            }
            set(m, a, x, y, &mut placed);
            queue.push_back(a);
        }
    }

    while let Some(u) = queue.pop_front() {
        let nbrs = m.neighbors(u);
        let new: Vec<(usize, u8)> = nbrs.iter().copied().filter(|x| !placed[x.0]).collect();
        if new.is_empty() {
            continue;
        }
        let (ux, uy) = (m.atoms[u].x, m.atoms[u].y);
        let existing: Vec<f32> = nbrs
            .iter()
            .filter(|x| placed[x.0])
            .map(|x| norm((m.atoms[x.0].y - uy).atan2(m.atoms[x.0].x - ux)))
            .collect();
        let k = new.len();
        let mut angles: Vec<f32> = vec![];
        if existing.is_empty() {
            for i in 0..k {
                angles.push(deg(-30.0) + 2.0 * PI * i as f32 / k as f32);
            }
        } else if existing.len() == 1 {
            let back = existing[0];
            if is_linear(m, u) && k == 1 {
                angles.push(back + PI);
            } else if k == 1 {
                // zig-zag: keep going the same vertical way as the bond we came from
                let c1 = back + deg(120.0);
                let c2 = back - deg(120.0);
                let sb = back.sin();
                if sb.abs() < 0.05 {
                    angles.push(if c1.sin() < 0.0 { c1 } else { c2 });
                } else if c1.sin() * sb > 0.0 {
                    angles.push(c1);
                } else {
                    angles.push(c2);
                }
            } else if k == 2 {
                angles.push(back + deg(120.0));
                angles.push(back - deg(120.0));
            } else {
                angles.push(back + deg(90.0));
                angles.push(back + deg(180.0));
                angles.push(back + deg(270.0));
            }
        } else {
            let mut e = existing.clone();
            e.sort_by(|a, b| a.partial_cmp(b).unwrap());
            // gaps
            let mut gaps: Vec<(f32, f32)> = vec![];
            for i in 0..e.len() {
                let a = e[i];
                let b = if i + 1 < e.len() { e[i + 1] } else { e[0] + 2.0 * PI };
                gaps.push((a, b - a));
            }
            // place one by one into the largest gap, splitting it
            let mut remaining = k;
            // distribute: count per gap proportional to size
            gaps.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
            let mut per = vec![0usize; gaps.len()];
            while remaining > 0 {
                let mut best = 0;
                let mut best_size = -1.0;
                for (i, g) in gaps.iter().enumerate() {
                    let size = g.1 / (per[i] + 1) as f32;
                    if size > best_size {
                        best_size = size;
                        best = i;
                    }
                }
                per[best] += 1;
                remaining -= 1;
            }
            for (i, g) in gaps.iter().enumerate() {
                for j in 0..per[i] {
                    angles.push(g.0 + g.1 * (j + 1) as f32 / (per[i] + 1) as f32);
                }
            }
        }
        for (idx, &(v, _)) in new.iter().enumerate() {
            if placed[v] {
                continue;
            }
            let a = angles[idx.min(angles.len() - 1)];
            let (vx, vy) = (ux + a.cos(), uy + a.sin());
            if let Some(ri) = ring_of(v) {
                if !placed_rings[ri] {
                    placed_rings[ri] = true;
                    let r = &rings[ri];
                    let k = r.len();
                    let rad = 0.5 / (PI / k as f32).sin();
                    let (cx, cy) = (vx + a.cos() * rad, vy + a.sin() * rad);
                    let start = r.iter().position(|&x| x == v).unwrap();
                    let base = a + PI; // angle from centre to v
                    for i in 0..k {
                        let at = r[(start + i) % k];
                        let ang = base + 2.0 * PI * i as f32 / k as f32;
                        set(m, at, cx + rad * ang.cos(), cy + rad * ang.sin(), &mut placed);
                        queue.push_back(at);
                    }
                    continue;
                }
            }
            set(m, v, vx, vy, &mut placed);
            queue.push_back(v);
        }
    }
    // atoms in other components (shouldn't happen) – place in a row
    let mut x = 0.0;
    for a in 0..n {
        if !placed[a] {
            m.atoms[a].x = x;
            m.atoms[a].y = 3.0;
            x += 1.0;
        }
    }
    resolve_clashes(m);
    center(m);
}

/// Nudge atoms that land on top of each other.
fn resolve_clashes(m: &mut Molecule) {
    let n = m.atoms.len();
    for _ in 0..30 {
        let mut moved = false;
        for i in 0..n {
            for j in i + 1..n {
                let dx = m.atoms[j].x - m.atoms[i].x;
                let dy = m.atoms[j].y - m.atoms[i].y;
                let d = (dx * dx + dy * dy).sqrt();
                if d < 0.45 {
                    let (ux, uy) = if d < 1e-4 { (0.3, 0.2) } else { (dx / d, dy / d) };
                    let push = (0.45 - d) * 0.5 + 0.02;
                    m.atoms[j].x += ux * push;
                    m.atoms[j].y += uy * push;
                    m.atoms[i].x -= ux * push;
                    m.atoms[i].y -= uy * push;
                    moved = true;
                }
            }
        }
        if !moved {
            break;
        }
    }
}

pub fn center(m: &mut Molecule) {
    if let Some((x0, y0, x1, y1)) = m.bounds() {
        let (cx, cy) = ((x0 + x1) / 2.0, (y0 + y1) / 2.0);
        for a in &mut m.atoms {
            a.x -= cx;
            a.y -= cy;
        }
    }
}
