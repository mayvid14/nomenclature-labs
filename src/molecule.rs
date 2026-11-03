//! Shared molecule model: atoms with 2D positions and bonds between them.

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Elem {
    C,
    H,
    O,
    N,
    F,
    Cl,
    Br,
    I,
}

impl Elem {
    pub fn symbol(self) -> &'static str {
        match self {
            Elem::C => "C",
            Elem::H => "H",
            Elem::O => "O",
            Elem::N => "N",
            Elem::F => "F",
            Elem::Cl => "Cl",
            Elem::Br => "Br",
            Elem::I => "I",
        }
    }
    pub fn valence(self) -> u32 {
        match self {
            Elem::C => 4,
            Elem::N => 3,
            Elem::O => 2,
            _ => 1,
        }
    }
    pub fn mass(self) -> f64 {
        match self {
            Elem::C => 12.011,
            Elem::H => 1.008,
            Elem::O => 15.999,
            Elem::N => 14.007,
            Elem::F => 18.998,
            Elem::Cl => 35.45,
            Elem::Br => 79.904,
            Elem::I => 126.904,
        }
    }
    pub fn is_halogen(self) -> bool {
        matches!(self, Elem::F | Elem::Cl | Elem::Br | Elem::I)
    }
    pub fn halo_prefix(self) -> &'static str {
        match self {
            Elem::F => "fluoro",
            Elem::Cl => "chloro",
            Elem::Br => "bromo",
            Elem::I => "iodo",
            _ => "",
        }
    }
    pub const ALL_DRAWABLE: [Elem; 8] = [
        Elem::C,
        Elem::H,
        Elem::O,
        Elem::N,
        Elem::F,
        Elem::Cl,
        Elem::Br,
        Elem::I,
    ];
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Atom {
    pub elem: Elem,
    pub x: f32,
    pub y: f32,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Bond {
    pub a: usize,
    pub b: usize,
    pub order: u8,
}

impl Bond {
    pub fn other(&self, i: usize) -> usize {
        if self.a == i { self.b } else { self.a }
    }
    pub fn has(&self, i: usize) -> bool {
        self.a == i || self.b == i
    }
}

/// Molecule in "bond length" units: a standard bond is 1.0 long.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
}

impl Molecule {
    pub fn add_atom(&mut self, elem: Elem, x: f32, y: f32) -> usize {
        self.atoms.push(Atom { elem, x, y });
        self.atoms.len() - 1
    }

    pub fn bond_between(&self, a: usize, b: usize) -> Option<usize> {
        self.bonds
            .iter()
            .position(|bd| (bd.a == a && bd.b == b) || (bd.a == b && bd.b == a))
    }

    pub fn add_bond(&mut self, a: usize, b: usize, order: u8) -> usize {
        if let Some(i) = self.bond_between(a, b) {
            self.bonds[i].order = order;
            return i;
        }
        self.bonds.push(Bond { a, b, order });
        self.bonds.len() - 1
    }

    pub fn neighbors(&self, i: usize) -> Vec<(usize, u8)> {
        self.bonds
            .iter()
            .filter(|b| b.has(i))
            .map(|b| (b.other(i), b.order))
            .collect()
    }

    pub fn bond_order_sum(&self, i: usize) -> u32 {
        self.bonds
            .iter()
            .filter(|b| b.has(i))
            .map(|b| b.order as u32)
            .sum()
    }

    pub fn implicit_h(&self, i: usize) -> u32 {
        let v = self.atoms[i].elem.valence();
        v.saturating_sub(self.bond_order_sum(i))
    }

    pub fn over_valent(&self, i: usize) -> bool {
        self.bond_order_sum(i) > self.atoms[i].elem.valence()
    }

    pub fn remove_atom(&mut self, i: usize) {
        self.bonds.retain(|b| !b.has(i));
        self.atoms.remove(i);
        for b in &mut self.bonds {
            if b.a > i {
                b.a -= 1;
            }
            if b.b > i {
                b.b -= 1;
            }
        }
    }

    pub fn remove_bond(&mut self, bi: usize) {
        let b = self.bonds.remove(bi);
        // Drop atoms that are now alone, so erasing a line removes its ends too.
        let mut ends = vec![b.a, b.b];
        ends.sort_unstable();
        for &e in ends.iter().rev() {
            if self.neighbors(e).is_empty() && self.atoms.len() > 1 {
                self.remove_atom(e);
            }
        }
    }

    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    /// Hill-order molecular formula pieces, e.g. [("C",2),("H",6),("O",1)].
    pub fn formula(&self) -> Vec<(&'static str, u32)> {
        let mut counts: std::collections::BTreeMap<&'static str, u32> = Default::default();
        let mut c = 0;
        let mut h = 0;
        for (i, a) in self.atoms.iter().enumerate() {
            match a.elem {
                Elem::C => c += 1,
                Elem::H => h += 1,
                e => *counts.entry(e.symbol()).or_default() += 1,
            }
            h += self.implicit_h(i);
        }
        let mut out = vec![];
        if c > 0 {
            out.push(("C", c));
            if h > 0 {
                out.push(("H", h));
            }
            out.extend(counts);
        } else {
            let mut all: Vec<(&'static str, u32)> = counts.into_iter().collect();
            if h > 0 {
                all.push(("H", h));
            }
            all.sort();
            out = all;
        }
        out
    }

    pub fn molar_mass(&self) -> f64 {
        let mut m = 0.0;
        for (i, a) in self.atoms.iter().enumerate() {
            m += a.elem.mass();
            m += self.implicit_h(i) as f64 * Elem::H.mass();
        }
        m
    }

    /// Bounding box in molecule units: (min_x, min_y, max_x, max_y).
    pub fn bounds(&self) -> Option<(f32, f32, f32, f32)> {
        if self.atoms.is_empty() {
            return None;
        }
        let mut b = (f32::MAX, f32::MAX, f32::MIN, f32::MIN);
        for a in &self.atoms {
            b.0 = b.0.min(a.x);
            b.1 = b.1.min(a.y);
            b.2 = b.2.max(a.x);
            b.3 = b.3.max(a.y);
        }
        Some(b)
    }

    /// Connected components as lists of atom indices.
    pub fn components(&self) -> Vec<Vec<usize>> {
        let n = self.atoms.len();
        let mut seen = vec![false; n];
        let mut out = vec![];
        for s in 0..n {
            if seen[s] {
                continue;
            }
            let mut stack = vec![s];
            seen[s] = true;
            let mut comp = vec![];
            while let Some(x) = stack.pop() {
                comp.push(x);
                for (y, _) in self.neighbors(x) {
                    if !seen[y] {
                        seen[y] = true;
                        stack.push(y);
                    }
                }
            }
            comp.sort_unstable();
            out.push(comp);
        }
        out
    }

    /// Copy of the atoms listed (in that order), with bonds among them.
    pub fn submolecule(&self, idx: &[usize]) -> (Molecule, Vec<usize>) {
        let mut map = vec![usize::MAX; self.atoms.len()];
        let mut m = Molecule::default();
        for (k, &i) in idx.iter().enumerate() {
            map[i] = k;
            m.atoms.push(self.atoms[i]);
        }
        for b in &self.bonds {
            if map[b.a] != usize::MAX && map[b.b] != usize::MAX {
                m.bonds.push(Bond { a: map[b.a], b: map[b.b], order: b.order });
            }
        }
        (m, idx.to_vec())
    }
}
