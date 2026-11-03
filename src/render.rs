//! Paints skeletal formulas with egui.

use crate::molecule::{Elem, Molecule};
use crate::namer::molecule_rings;
use egui::text::LayoutJob;
use egui::{Align, Color32, FontId, Painter, Pos2, Stroke, TextFormat, Vec2, pos2, vec2};

pub const INK: Color32 = Color32::from_rgb(31, 41, 51);
pub const MUTED: Color32 = Color32::from_rgb(100, 112, 128);
pub const ACCENT: Color32 = Color32::from_rgb(13, 148, 136);
pub const ACCENT_DARK: Color32 = Color32::from_rgb(15, 118, 110);
pub const ERR: Color32 = Color32::from_rgb(220, 38, 38);
pub const BG: Color32 = Color32::from_rgb(243, 245, 248);
pub const CARD: Color32 = Color32::WHITE;
pub const BORDER: Color32 = Color32::from_rgb(222, 226, 232);

pub fn accent_soft(a: u8) -> Color32 {
    Color32::from_rgba_unmultiplied(13, 148, 136, a)
}

pub fn elem_color(e: Elem) -> Color32 {
    match e {
        Elem::C => INK,
        Elem::H => Color32::from_rgb(90, 98, 110),
        Elem::O => Color32::from_rgb(214, 40, 40),
        Elem::N => Color32::from_rgb(37, 99, 235),
        Elem::F => Color32::from_rgb(101, 163, 13),
        Elem::Cl => Color32::from_rgb(22, 150, 60),
        Elem::Br => Color32::from_rgb(154, 52, 18),
        Elem::I => Color32::from_rgb(126, 34, 206),
    }
}

pub struct DrawOpts<'a> {
    pub show_carbons: bool,
    pub parent: &'a [usize],
    pub parent_ring: bool,
    pub show_numbers: bool,
    pub error_atoms: &'a [usize],
    pub hover_atom: Option<usize>,
    pub hover_bond: Option<usize>,
}

/// Label pieces (text, is_subscript, is_main_symbol); None for bare skeletal carbons.
pub fn atom_label(m: &Molecule, i: usize, show_carbons: bool) -> Option<Vec<(String, bool, bool)>> {
    let a = m.atoms[i];
    let nb = m.neighbors(i);
    if a.elem == Elem::C && !show_carbons && !nb.is_empty() {
        return None;
    }
    let h = m.implicit_h(i);
    let mut main = vec![(a.elem.symbol().to_string(), false, true)];
    let mut hs = vec![];
    if h > 0 && a.elem != Elem::H {
        hs.push(("H".to_string(), false, false));
        if h > 1 {
            hs.push((h.to_string(), true, false));
        }
    }
    if hs.is_empty() {
        return Some(main);
    }
    // put H on the left when the bonds come from the right
    let sum_dx: f32 = nb.iter().map(|(j, _)| m.atoms[*j].x - a.x).sum();
    let left = !nb.is_empty() && sum_dx > 0.25 && a.elem != Elem::C;
    if left {
        hs.extend(main);
        Some(hs)
    } else {
        main.extend(hs);
        Some(main)
    }
}

fn label_job(pieces: &[(String, bool, bool)], size: f32, color: Color32) -> LayoutJob {
    let mut job = LayoutJob::default();
    for (t, sub, _) in pieces {
        let fmt = TextFormat {
            font_id: FontId::proportional(if *sub { size * 0.68 } else { size }),
            color,
            valign: if *sub { Align::BOTTOM } else { Align::Center },
            ..Default::default()
        };
        job.append(t, 0.0, fmt);
    }
    job
}

/// Rich text with subscript digits, e.g. a formula "C2H6O".
pub fn formula_job(parts: &[(&str, u32)], size: f32, color: Color32) -> LayoutJob {
    let mut job = LayoutJob::default();
    for (sym, n) in parts {
        job.append(
            sym,
            0.0,
            TextFormat { font_id: FontId::proportional(size), color, valign: Align::Center, ..Default::default() },
        );
        if *n > 1 {
            job.append(
                &n.to_string(),
                0.0,
                TextFormat {
                    font_id: FontId::proportional(size * 0.7),
                    color,
                    valign: Align::BOTTOM,
                    ..Default::default()
                },
            );
        }
    }
    job
}

pub fn draw_molecule(painter: &Painter, m: &Molecule, to_screen: &dyn Fn(f32, f32) -> Pos2, scale: f32, o: &DrawOpts) {
    let n = m.atoms.len();
    let pos: Vec<Pos2> = m.atoms.iter().map(|a| to_screen(a.x, a.y)).collect();
    let rings = molecule_rings(m);
    let stroke_w = (scale * 0.05).clamp(1.4, 2.6);
    let font = (scale * 0.4).clamp(11.0, 22.0);

    // main-chain highlight
    if !o.parent.is_empty() {
        let soft = accent_soft(46);
        let k = o.parent.len();
        let segs = if o.parent_ring { k } else { k.saturating_sub(1) };
        for i in 0..segs {
            let (a, b) = (o.parent[i], o.parent[(i + 1) % k]);
            if a < n && b < n {
                painter.line_segment([pos[a], pos[b]], Stroke::new(scale * 0.26, soft));
            }
        }
        for &a in o.parent {
            if a < n {
                painter.circle_filled(pos[a], scale * 0.13, soft);
            }
        }
    }
    if let Some(bi) = o.hover_bond {
        if let Some(b) = m.bonds.get(bi) {
            painter.line_segment([pos[b.a], pos[b.b]], Stroke::new(scale * 0.3, accent_soft(70)));
        }
    }
    if let Some(ai) = o.hover_atom {
        if ai < n {
            painter.circle_filled(pos[ai], scale * 0.3, accent_soft(70));
        }
    }

    let labels: Vec<Option<Vec<(String, bool, bool)>>> = (0..n).map(|i| atom_label(m, i, o.show_carbons)).collect();

    // bonds
    for b in &m.bonds {
        let (mut p1, mut p2) = (pos[b.a], pos[b.b]);
        let d = p2 - p1;
        let len = d.length();
        if len < 1.0 {
            continue;
        }
        let u = d / len;
        let shrink = scale * 0.3;
        if labels[b.a].is_some() {
            p1 += u * shrink.min(len * 0.4);
        }
        if labels[b.b].is_some() {
            p2 -= u * shrink.min(len * 0.4);
        }
        let normal = vec2(-u.y, u.x);
        let col = if o.error_atoms.contains(&b.a) && o.error_atoms.contains(&b.b) { ERR } else { INK };
        let st = Stroke::new(stroke_w, col);
        let gap = scale * 0.15;
        match b.order {
            2 => {
                let ring = rings.iter().find(|r| r.contains(&b.a) && r.contains(&b.b));
                if let Some(r) = ring {
                    // inner line towards the ring centre
                    let c = r.iter().fold(Vec2::ZERO, |acc, &i| acc + pos[i].to_vec2()) / r.len() as f32;
                    let mid = (p1 + p2.to_vec2()) / 2.0;
                    let side = if (c - mid.to_vec2()).dot(normal) > 0.0 { 1.0 } else { -1.0 };
                    painter.line_segment([p1, p2], st);
                    let off = normal * gap * side;
                    let q1 = p1 + off + u * len * 0.14;
                    let q2 = p2 + off - u * len * 0.14;
                    painter.line_segment([q1, q2], st);
                } else {
                    let off = normal * gap * 0.5;
                    painter.line_segment([p1 + off, p2 + off], st);
                    painter.line_segment([p1 - off, p2 - off], st);
                }
            }
            3 => {
                let off = normal * gap * 0.85;
                painter.line_segment([p1, p2], st);
                painter.line_segment([p1 + off, p2 + off], st);
                painter.line_segment([p1 - off, p2 - off], st);
            }
            _ => {
                painter.line_segment([p1, p2], st);
            }
        }
    }

    // labels
    for i in 0..n {
        let Some(pieces) = &labels[i] else { continue };
        let color = if o.error_atoms.contains(&i) { ERR } else { elem_color(m.atoms[i].elem) };
        let galley = painter.layout_job(label_job(pieces, font, color));
        // centre the main symbol on the atom
        let mut before = 0.0;
        let mut main_w = 0.0;
        for (t, sub, is_main) in pieces {
            let g = painter.layout_no_wrap(
                t.clone(),
                FontId::proportional(if *sub { font * 0.68 } else { font }),
                color,
            );
            if *is_main {
                main_w = g.size().x;
                break;
            }
            before += g.size().x;
        }
        let top_left = pos2(pos[i].x - before - main_w / 2.0, pos[i].y - galley.size().y / 2.0);
        let bg_rect = egui::Rect::from_min_size(top_left, galley.size()).expand(1.5);
        painter.rect_filled(bg_rect, 3.0, Color32::from_rgba_unmultiplied(255, 255, 255, 200));
        painter.galley(top_left, galley, color);
    }

    // error rings
    for &e in o.error_atoms {
        if e < n {
            painter.circle_stroke(pos[e], scale * 0.38, Stroke::new(2.0, ERR));
        }
    }

    // locant bubbles
    if o.show_numbers && o.parent.len() > 1 {
        let ring_centre = if o.parent_ring {
            let v: Vec<Pos2> = o.parent.iter().filter(|&&a| a < n).map(|&a| pos[a]).collect();
            Some(v.iter().fold(Vec2::ZERO, |acc, p| acc + p.to_vec2()) / v.len().max(1) as f32)
        } else {
            None
        };
        for (k, &a) in o.parent.iter().enumerate() {
            if a >= n {
                continue;
            }
            // bisector of the widest empty angle around the atom
            let mut angs: Vec<f32> = m
                .neighbors(a)
                .iter()
                .map(|(j, _)| {
                    let d = pos[*j] - pos[a];
                    d.y.atan2(d.x)
                })
                .collect();
            let dir = if angs.is_empty() {
                vec2(0.0, -1.0)
            } else {
                angs.sort_by(|x, y| x.partial_cmp(y).unwrap());
                let mut best = (0.0f32, -1.0f32);
                for i in 0..angs.len() {
                    let s = angs[i];
                    let e = if i + 1 < angs.len() { angs[i + 1] } else { angs[0] + std::f32::consts::TAU };
                    if e - s > best.1 + 1e-3 {
                        best = (s + (e - s) / 2.0, e - s);
                    }
                }
                vec2(best.0.cos(), best.0.sin())
            };
            let r = if ring_centre.is_some() { (scale * 0.16).clamp(6.0, 9.5) } else { (scale * 0.2).clamp(7.0, 11.0) };
            let at = match ring_centre {
                Some(c) => pos[a] + (c - pos[a].to_vec2()).normalized() * scale * 0.4,
                None => pos[a] + dir * (scale * 0.5 + if labels[a].is_some() { scale * 0.15 } else { 0.0 }),
            };
            painter.circle_filled(at, r, ACCENT);
            painter.text(
                at,
                egui::Align2::CENTER_CENTER,
                (k + 1).to_string(),
                FontId::proportional(r * 1.25),
                Color32::WHITE,
            );
        }
    }
}
