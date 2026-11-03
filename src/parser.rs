//! IUPAC name → structure. Accepts both classic ("2-butene", "2-propanol")
//! and modern ("but-2-ene", "propan-2-ol") styles, plus a list of common names.

use crate::molecule::{Elem, Molecule};
use crate::namer::{self, Suffix};

#[derive(Clone, Copy, Debug, PartialEq)]
enum Loc {
    Num(u32),
    N,
}

#[derive(Clone, Debug)]
enum Sub {
    Halo(Elem),
    Hydroxy,
    Oxo,
    Amino(Vec<Sub>),
    Alkoxy(Box<Sub>),
    Carboxy,
    Cyano,
    Formyl,
    Carbamoyl(Vec<Sub>),
    AlkoxyCarbonyl(Box<Sub>),
    Carbon(Frag),
}

impl Sub {
    fn attach_order(&self) -> u32 {
        match self {
            Sub::Oxo => 2,
            Sub::Carbon(f) => f.attach_order as u32,
            _ => 1,
        }
    }
    fn is_carbon_branch(&self) -> bool {
        matches!(self, Sub::Carbon(_))
    }
}

#[derive(Clone, Debug, Default)]
struct Unsat {
    doubles: Option<Vec<u32>>,
    n_double: usize,
    triples: Option<Vec<u32>>,
    n_triple: usize,
}

#[derive(Clone, Debug)]
struct Frag {
    len: usize,
    ring: bool,
    aromatic: bool,
    unsat: Unsat,
    attach: Option<u32>,
    attach_order: u8,
    prefixes: Vec<Item>,
}

impl Frag {
    fn chain(len: usize) -> Frag {
        Frag {
            len,
            ring: false,
            aromatic: false,
            unsat: Unsat::default(),
            attach: Some(1),
            attach_order: 1,
            prefixes: vec![],
        }
    }
    fn phenyl() -> Frag {
        Frag { len: 6, ring: true, aromatic: true, ..Frag::chain(6) }
    }
}

#[derive(Clone, Debug)]
struct Item {
    locs: Option<Vec<Loc>>,
    count: usize,
    sub: Sub,
}

#[derive(Clone, Debug)]
struct SuffixSpec {
    kind: Suffix,
    count: usize,
    locs: Option<Vec<u32>>,
    carbo: bool,
}

#[derive(Clone, Debug)]
struct Parent {
    len: usize,
    ring: bool,
    aromatic: bool,
    unsat: Unsat,
    suffix: Option<SuffixSpec>,
}

pub struct Built {
    pub mol: Molecule,
    pub backbone: Vec<usize>,
    pub backbone_ring: bool,
    /// The name as it was understood (after common-name lookup).
    pub interpreted: String,
}

// ---------------------------------------------------------------- lexing helpers

fn parse_number_list(s: &str) -> Option<(Vec<Loc>, &str)> {
    let mut locs = vec![];
    let mut rest = s;
    loop {
        let digits: String = rest.chars().take_while(|c| c.is_ascii_digit()).collect();
        if !digits.is_empty() {
            locs.push(Loc::Num(digits.parse().ok()?));
            rest = &rest[digits.len()..];
            // primes like 1' are not supported
        } else if let Some(r) = rest.strip_prefix('n') {
            // "N" locant, must be followed by , or -
            if r.starts_with(',') || r.starts_with('-') {
                locs.push(Loc::N);
                rest = r;
            } else {
                return None;
            }
        } else {
            return None;
        }
        if let Some(r) = rest.strip_prefix(',') {
            rest = r;
        } else {
            break;
        }
    }
    Some((locs, rest))
}

/// "2,3-" at the start of s.
fn parse_locants(s: &str) -> Option<(Vec<Loc>, &str)> {
    let (l, rest) = parse_number_list(s)?;
    let rest = rest.strip_prefix('-')?;
    Some((l, rest))
}

/// "-2,3-" inside a word.
fn parse_inner_locants(s: &str) -> Option<(Vec<u32>, &str)> {
    let r = s.strip_prefix('-')?;
    let (l, rest) = parse_number_list(r)?;
    let rest = rest.strip_prefix('-')?;
    let nums: Option<Vec<u32>> =
        l.iter().map(|x| if let Loc::Num(n) = x { Some(*n) } else { None }).collect();
    Some((nums?, rest))
}

const MULTS: [(&str, usize); 9] = [
    ("di", 2),
    ("tri", 3),
    ("tetra", 4),
    ("penta", 5),
    ("hexa", 6),
    ("bis", 2),
    ("tris", 3),
    ("tetrakis", 4),
    ("pentakis", 5),
];

fn mult_options(s: &str) -> Vec<(usize, &str)> {
    let mut v = vec![(1, s)];
    for (w, k) in MULTS {
        if let Some(r) = s.strip_prefix(w) {
            v.push((k, r));
        }
    }
    v
}

/// Stems that s starts with (longest first).
fn stem_options(s: &str) -> Vec<(usize, &str)> {
    let mut v = vec![];
    for (i, st) in namer::stems().iter().enumerate() {
        if let Some(r) = s.strip_prefix(st) {
            v.push((i + 1, r));
        }
    }
    v.sort_by_key(|x| std::cmp::Reverse(x.0));
    v
}

/// Unsaturation after a stem: "an", "-2-en", "a-1,3-dien", "-1-en-4-yn", "en", "yn", "adien".
/// Returns (unsat, rest) where rest is right after "an"/"en"/"yn".
fn unsat_options(s: &str) -> Vec<(Unsat, &str)> {
    let mut out = vec![];
    if let Some(r) = s.strip_prefix("an") {
        out.push((Unsat::default(), r));
    }
    // en / yn parts
    let starts: Vec<&str> = if let Some(r) = s.strip_prefix('a') { vec![s, r] } else { vec![s] };
    for st in starts {
        // optional locants
        let opts: Vec<(Option<Vec<u32>>, &str)> = {
            let mut o = vec![(None, st)];
            if let Some((l, r)) = parse_inner_locants(st) {
                o.push((Some(l), r));
            }
            o
        };
        for (locs, r) in opts {
            for (k, r2) in mult_options(r) {
                if k > 1 && !r.starts_with(|c: char| c.is_alphabetic()) {
                    continue;
                }
                if let Some(r3) = r2.strip_prefix("en") {
                    let u = Unsat { doubles: locs.clone(), n_double: locs.as_ref().map_or(k, |l| l.len().max(k)), ..Default::default() };
                    out.push((u.clone(), r3));
                    // optional yne after ene
                    let yopts: Vec<(Option<Vec<u32>>, &str)> = {
                        let mut o = vec![(None, r3)];
                        if let Some((l, rr)) = parse_inner_locants(r3) {
                            o.push((Some(l), rr));
                        }
                        o
                    };
                    for (tl, rr) in yopts {
                        for (kt, rr2) in mult_options(rr) {
                            if let Some(rr3) = rr2.strip_prefix("yn") {
                                let mut u2 = u.clone();
                                u2.n_triple = tl.as_ref().map_or(kt, |l| l.len().max(kt));
                                u2.triples = tl.clone();
                                out.push((u2, rr3));
                            }
                        }
                    }
                }
                if let Some(r3) = r2.strip_prefix("yn") {
                    let u = Unsat { triples: locs.clone(), n_triple: locs.as_ref().map_or(k, |l| l.len().max(k)), ..Default::default() };
                    out.push((u, r3));
                }
            }
        }
    }
    out
}

const SUFFIX_WORDS: [(&str, Suffix, bool); 13] = [
    ("carboxylic acid", Suffix::Acid, true),
    ("carboxylate", Suffix::Ester, true),
    ("carboxamide", Suffix::Amide, true),
    ("carbonitrile", Suffix::Nitrile, true),
    ("carbaldehyde", Suffix::Aldehyde, true),
    ("oic acid", Suffix::Acid, false),
    ("oate", Suffix::Ester, false),
    ("amide", Suffix::Amide, false),
    ("amine", Suffix::Amine, false),
    ("nitrile", Suffix::Nitrile, false),
    ("one", Suffix::Ketone, false),
    ("ol", Suffix::Alcohol, false),
    ("al", Suffix::Aldehyde, false),
];

/// Suffix part, must consume all of s. Ok(None) = no suffix.
fn suffix_full(s: &str) -> Option<Option<SuffixSpec>> {
    let s = s.strip_prefix('e').unwrap_or(s);
    if s.is_empty() {
        return Some(None);
    }
    let mut opts: Vec<(Option<Vec<u32>>, &str)> = vec![(None, s)];
    if let Some((l, r)) = parse_inner_locants(s) {
        opts.push((Some(l), r));
    }
    for (locs, r) in opts {
        for (k, r2) in mult_options(r) {
            for (w, kind, carbo) in SUFFIX_WORDS {
                if r2 == w {
                    let count = locs.as_ref().map_or(k, |l| l.len().max(k));
                    return Some(Some(SuffixSpec { kind, count, locs: locs.clone(), carbo }));
                }
            }
        }
    }
    None
}

fn parent_full(s: &str) -> Option<Parent> {
    // benzene family
    let retained: [(&str, Option<Suffix>); 8] = [
        ("benzene", None),
        ("phenol", Some(Suffix::Alcohol)),
        ("benzoic acid", Some(Suffix::Acid)),
        ("benzoate", Some(Suffix::Ester)),
        ("benzaldehyde", Some(Suffix::Aldehyde)),
        ("benzonitrile", Some(Suffix::Nitrile)),
        ("benzamide", Some(Suffix::Amide)),
        ("aniline", Some(Suffix::Amine)),
    ];
    for (w, suf) in retained {
        if s == w {
            return Some(Parent {
                len: 6,
                ring: true,
                aromatic: true,
                unsat: Unsat::default(),
                suffix: suf.map(|kind| SuffixSpec {
                    kind,
                    count: 1,
                    locs: Some(vec![1]),
                    carbo: matches!(kind, Suffix::Acid | Suffix::Ester | Suffix::Aldehyde | Suffix::Nitrile | Suffix::Amide),
                }),
            });
        }
    }
    if let Some(r) = s.strip_prefix("benzen") {
        if let Some(suf) = suffix_full(r) {
            return Some(Parent { len: 6, ring: true, aromatic: true, unsat: Unsat::default(), suffix: suf });
        }
    }
    let (ring, s2) = match s.strip_prefix("cyclo") {
        Some(r) => (true, r),
        None => (false, s),
    };
    for (n, r) in stem_options(s2) {
        for (u, r2) in unsat_options(r) {
            if let Some(suf) = suffix_full(r2) {
                // "methan" alone is not a name
                if suf.is_none() && !r2.starts_with('e') {
                    continue;
                }
                return Some(Parent { len: n, ring, aromatic: false, unsat: u, suffix: suf });
            }
        }
    }
    None
}

// ---------------------------------------------------------------- substituent names

fn simple_words() -> Vec<(&'static str, Sub)> {
    let isopropyl = Sub::Carbon(Frag { attach: Some(2), ..Frag::chain(3) });
    let tbutyl = Sub::Carbon(Frag {
        prefixes: vec![Item { locs: Some(vec![Loc::Num(1), Loc::Num(1)]), count: 2, sub: Sub::Carbon(Frag::chain(1)) }],
        ..Frag::chain(2)
    });
    let sbutyl = Sub::Carbon(Frag { attach: Some(2), ..Frag::chain(4) });
    let ibutyl = Sub::Carbon(Frag {
        prefixes: vec![Item { locs: Some(vec![Loc::Num(2)]), count: 1, sub: Sub::Carbon(Frag::chain(1)) }],
        ..Frag::chain(3)
    });
    let vinyl = Sub::Carbon(Frag {
        unsat: Unsat { doubles: Some(vec![1]), n_double: 1, ..Default::default() },
        ..Frag::chain(2)
    });
    let allyl = Sub::Carbon(Frag {
        unsat: Unsat { doubles: Some(vec![2]), n_double: 1, ..Default::default() },
        ..Frag::chain(3)
    });
    let benzyl = Sub::Carbon(Frag {
        prefixes: vec![Item { locs: Some(vec![Loc::Num(1)]), count: 1, sub: Sub::Carbon(Frag::phenyl()) }],
        ..Frag::chain(1)
    });
    let acetyl = Sub::Carbon(Frag {
        prefixes: vec![Item { locs: Some(vec![Loc::Num(1)]), count: 1, sub: Sub::Oxo }],
        ..Frag::chain(2)
    });
    let ylidene = |n| Sub::Carbon(Frag { attach_order: 2, ..Frag::chain(n) });
    vec![
        ("fluoro", Sub::Halo(Elem::F)),
        ("chloro", Sub::Halo(Elem::Cl)),
        ("bromo", Sub::Halo(Elem::Br)),
        ("iodo", Sub::Halo(Elem::I)),
        ("hydroxy", Sub::Hydroxy),
        ("hydroxyl", Sub::Hydroxy),
        ("oxo", Sub::Oxo),
        ("keto", Sub::Oxo),
        ("amino", Sub::Amino(vec![])),
        ("cyano", Sub::Cyano),
        ("carboxy", Sub::Carboxy),
        ("formyl", Sub::Formyl),
        ("carbamoyl", Sub::Carbamoyl(vec![])),
        ("phenyl", Sub::Carbon(Frag::phenyl())),
        ("phenoxy", Sub::Alkoxy(Box::new(Sub::Carbon(Frag::phenyl())))),
        ("benzyl", benzyl),
        ("vinyl", vinyl.clone()),
        ("ethenyl", vinyl),
        ("allyl", allyl),
        ("isopropyl", isopropyl.clone()),
        ("tert-butyl", tbutyl.clone()),
        ("t-butyl", tbutyl.clone()),
        ("tertbutyl", tbutyl.clone()),
        ("sec-butyl", sbutyl.clone()),
        ("s-butyl", sbutyl),
        ("isobutyl", ibutyl),
        ("isopropoxy", Sub::Alkoxy(Box::new(isopropyl))),
        ("tert-butoxy", Sub::Alkoxy(Box::new(tbutyl))),
        ("acetyl", acetyl),
        ("methylene", ylidene(1)),
        ("ethynyl", Sub::Carbon(Frag {
            unsat: Unsat { triples: Some(vec![1]), n_triple: 1, ..Default::default() },
            ..Frag::chain(2)
        })),
    ]
}

/// Alkyl-type substituent bases without leading locants:
/// "methyl", "cyclohexyl", "propan-2-yl", "prop-2-en-1-yl", "propenyl", "methoxy", "ethylidene".
fn alkyl_options(s: &str) -> Vec<(Sub, &str)> {
    let mut out = vec![];
    let (ring, s2) = match s.strip_prefix("cyclo") {
        Some(r) => (true, r),
        None => (false, s),
    };
    for (n, r) in stem_options(s2) {
        let mk = |u: Unsat, attach: Option<u32>, order: u8| Frag {
            len: n,
            ring,
            aromatic: false,
            unsat: u,
            attach: if ring { Some(1) } else { attach },
            attach_order: order,
            prefixes: vec![],
        };
        // plain: methyl / methylidene
        if let Some(r2) = r.strip_prefix("ylidene") {
            out.push((Sub::Carbon(mk(Unsat::default(), Some(1), 2)), r2));
        }
        if let Some(r2) = r.strip_prefix("yl") {
            out.push((Sub::Carbon(mk(Unsat::default(), Some(1), 1)), r2));
        }
        if let Some(r2) = r.strip_prefix("oxy") {
            out.push((Sub::Alkoxy(Box::new(Sub::Carbon(mk(Unsat::default(), Some(1), 1)))), r2));
        }
        for (u, r2) in unsat_options(r) {
            // propan-2-yl, prop-2-en-1-yl
            if let Some((l, r3)) = parse_inner_locants(r2) {
                if l.len() == 1 {
                    if let Some(r4) = r3.strip_prefix("yl") {
                        let r4b = r4.strip_prefix("idene");
                        let order = if r4b.is_some() { 2 } else { 1 };
                        out.push((Sub::Carbon(mk(u.clone(), Some(l[0]), order)), r4b.unwrap_or(r4)));
                    }
                }
            }
            // propenyl, ethenyl (attach 1)
            if let Some(r3) = r2.strip_prefix("yl") {
                if u.n_double + u.n_triple > 0 {
                    out.push((Sub::Carbon(mk(u.clone(), Some(1), 1)), r3));
                }
            }
        }
    }
    out
}

fn sub_options(s: &str) -> Vec<(Sub, &str)> {
    let mut out = vec![];
    if let Some(inner_start) = s.strip_prefix('(') {
        // matching paren
        let mut depth = 1;
        for (i, ch) in inner_start.char_indices() {
            match ch {
                '(' => depth += 1,
                ')' => {
                    depth -= 1;
                    if depth == 0 {
                        let inner = &inner_start[..i];
                        let rest = &inner_start[i + 1..];
                        if let Some(sub) = complex_sub(inner) {
                            out.push((sub, rest));
                        }
                        break;
                    }
                }
                _ => {}
            }
        }
        return out;
    }
    for (w, sub) in simple_words() {
        if let Some(r) = s.strip_prefix(w) {
            out.push((sub, r));
        }
    }
    for (sub, r) in alkyl_options(s) {
        // alkoxycarbonyl
        if let (Sub::Alkoxy(_), Some(r2)) = (&sub, r.strip_prefix("carbonyl")) {
            if let Sub::Alkoxy(inner) = sub.clone() {
                out.push((Sub::AlkoxyCarbonyl(inner), r2));
            }
        }
        out.push((sub, r));
    }
    // simple words followed by "oxy" / "carbonyl" (e.g. "phenoxycarbonyl")
    let mut extra = vec![];
    for (sub, r) in &out {
        if let Sub::Alkoxy(inner) = sub {
            if let Some(r2) = r.strip_prefix("carbonyl") {
                extra.push((Sub::AlkoxyCarbonyl(inner.clone()), r2));
            }
        }
        if let Sub::Carbon(_) = sub {
            if let Some(r2) = r.strip_prefix("oxy") {
                extra.push((Sub::Alkoxy(Box::new(sub.clone())), r2));
            }
        }
    }
    out.extend(extra);
    out
}

fn mult_sub_options(s: &str) -> Vec<(usize, Sub, &str)> {
    let mut out = vec![];
    for (k, r) in mult_options(s) {
        for (sub, r2) in sub_options(r) {
            out.push((k, sub, r2));
        }
    }
    out
}

/// Items followed by a terminal piece that must consume the rest.
fn parse_seq<T>(s: &str, end: &dyn Fn(&str) -> Option<T>, depth: usize) -> Option<(Vec<Item>, Option<Vec<Loc>>, T)> {
    if depth > 40 {
        return None;
    }
    if let Some(t) = end(s) {
        return Some((vec![], None, t));
    }
    if let Some((locs, rest)) = parse_locants(s) {
        if let Some(t) = end(rest) {
            return Some((vec![], Some(locs), t));
        }
        for (k, sub, r2) in mult_sub_options(rest) {
            let r2 = r2.strip_prefix('-').unwrap_or(r2);
            if let Some((mut items, pend, t)) = parse_seq(r2, end, depth + 1) {
                let count = k.max(locs.len());
                items.insert(0, Item { locs: Some(locs.clone()), count, sub });
                return Some((items, pend, t));
            }
        }
    }
    for (k, sub, r2) in mult_sub_options(s) {
        let r2 = r2.strip_prefix('-').unwrap_or(r2);
        if let Some((mut items, pend, t)) = parse_seq(r2, end, depth + 1) {
            items.insert(0, Item { locs: None, count: k, sub });
            return Some((items, pend, t));
        }
    }
    None
}

enum Base {
    Sub(Sub),
    Amino,
    Carbamoyl,
}

fn base_full(s: &str) -> Option<Base> {
    if s == "amino" {
        return Some(Base::Amino);
    }
    if s == "carbamoyl" {
        return Some(Base::Carbamoyl);
    }
    for (sub, r) in sub_options(s) {
        if r.is_empty() && !s.starts_with('(') {
            return Some(Base::Sub(sub));
        }
    }
    None
}

fn flatten(items: &[Item]) -> Vec<Sub> {
    let mut v = vec![];
    for it in items {
        for _ in 0..it.count {
            v.push(it.sub.clone());
        }
    }
    v
}

fn assign_pending_unsat(u: &mut Unsat, pend: &[u32]) -> bool {
    if u.n_double > 0 && u.doubles.is_none() {
        u.doubles = Some(pend.to_vec());
        u.n_double = u.n_double.max(pend.len());
        true
    } else if u.n_triple > 0 && u.triples.is_none() {
        u.triples = Some(pend.to_vec());
        u.n_triple = u.n_triple.max(pend.len());
        true
    } else {
        false
    }
}

fn nums(l: &[Loc]) -> Option<Vec<u32>> {
    l.iter().map(|x| if let Loc::Num(n) = x { Some(*n) } else { None }).collect()
}

fn complex_sub(inner: &str) -> Option<Sub> {
    let (items, pend, base) = parse_seq(inner, &base_full, 0)?;
    match base {
        Base::Amino => Some(Sub::Amino(flatten(&items))),
        Base::Carbamoyl => Some(Sub::Carbamoyl(flatten(&items))),
        Base::Sub(sub) => {
            let apply = |f: &mut Frag| -> Option<()> {
                f.prefixes.extend(items.iter().cloned());
                if let Some(p) = &pend {
                    let p = nums(p)?;
                    if !assign_pending_unsat(&mut f.unsat, &p) {
                        // "(2-propyl)" style attachment point
                        if p.len() == 1 && !f.ring {
                            f.attach = Some(p[0]);
                        } else {
                            return None;
                        }
                    }
                }
                Some(())
            };
            match sub {
                Sub::Carbon(mut f) => {
                    apply(&mut f)?;
                    Some(Sub::Carbon(f))
                }
                Sub::Alkoxy(b) => match *b {
                    Sub::Carbon(mut f) => {
                        apply(&mut f)?;
                        Some(Sub::Alkoxy(Box::new(Sub::Carbon(f))))
                    }
                    other => if items.is_empty() { Some(Sub::Alkoxy(Box::new(other))) } else { None },
                },
                Sub::AlkoxyCarbonyl(b) => match *b {
                    Sub::Carbon(mut f) => {
                        apply(&mut f)?;
                        Some(Sub::AlkoxyCarbonyl(Box::new(Sub::Carbon(f))))
                    }
                    other => if items.is_empty() { Some(Sub::AlkoxyCarbonyl(Box::new(other))) } else { None },
                },
                Sub::Amino(_) if items.is_empty() => Some(Sub::Amino(vec![])),
                other => if items.is_empty() && pend.is_none() { Some(other) } else { None },
            }
        }
    }
}

// ---------------------------------------------------------------- building

struct Builder {
    mol: Molecule,
}

impl Builder {
    fn add(&mut self, e: Elem) -> usize {
        self.mol.add_atom(e, 0.0, 0.0)
    }
    fn free(&self, a: usize) -> i32 {
        self.mol.atoms[a].elem.valence() as i32 - self.mol.bond_order_sum(a) as i32
    }
    fn bond(&mut self, a: usize, b: usize, o: u8) {
        self.mol.add_bond(a, b, o);
    }

    /// Creates a chain/ring skeleton, returns atom indices (locant order).
    fn skeleton(&mut self, len: usize, ring: bool, aromatic: bool, unsat: &Unsat) -> Result<Vec<usize>, String> {
        if len == 0 {
            return Err("Empty chain.".into());
        }
        if ring && len < 3 {
            return Err("A ring needs at least 3 carbons.".into());
        }
        let atoms: Vec<usize> = (0..len).map(|_| self.add(Elem::C)).collect();
        let nb = if ring { len } else { len - 1 };
        for i in 0..nb {
            self.bond(atoms[i], atoms[(i + 1) % len], 1);
        }
        if aromatic {
            for i in [0, 2, 4] {
                self.bond(atoms[i], atoms[i + 1], 2);
            }
            return Ok(atoms);
        }
        let place = |locs: &Option<Vec<u32>>, count: usize, order: u8, b: &mut Builder| -> Result<(), String> {
            let what = if order == 2 { "double" } else { "triple" };
            match locs {
                Some(ls) => {
                    for &l in ls {
                        if l < 1 || l as usize > nb {
                            return Err(format!(
                                "There's no bond number {l} for a {what} bond in a {len}-carbon {}.",
                                if ring { "ring" } else { "chain" }
                            ));
                        }
                        let (x, y) = (atoms[l as usize - 1], atoms[l as usize % len]);
                        b.bond(x, y, order);
                    }
                }
                None => {
                    // lowest positions, avoiding neighbours of existing multiple bonds first
                    let mut placed = 0;
                    for pass in 0..2 {
                        for i in 0..nb {
                            if placed == count {
                                break;
                            }
                            let (x, y) = (atoms[i], atoms[(i + 1) % len]);
                            if b.mol.bond_between(x, y).map(|bi| b.mol.bonds[bi].order) != Some(1) {
                                continue;
                            }
                            let need = order as i32 - 1;
                            if b.free(x) < need || b.free(y) < need {
                                continue;
                            }
                            let busy = |a: usize, b: &Builder| b.mol.neighbors(a).iter().any(|n| n.1 > 1);
                            if pass == 0 && (busy(x, b) || busy(y, b)) {
                                continue;
                            }
                            b.bond(x, y, order);
                            placed += 1;
                        }
                    }
                    if placed < count {
                        return Err(format!("Not enough room for {count} {what} bonds."));
                    }
                }
            }
            Ok(())
        };
        place(&unsat.doubles, unsat.n_double, 2, self)?;
        place(&unsat.triples, unsat.n_triple, 3, self)?;
        Ok(atoms)
    }

    /// Choose `count` positions for a group needing `need` free valence.
    fn default_positions(&self, atoms: &[usize], count: usize, need: i32, prefer_inner: bool, start: usize) -> Result<Vec<usize>, String> {
        let mut out = vec![];
        let mut used = vec![0i32; atoms.len()];
        let order: Vec<usize> = if prefer_inner && atoms.len() > 2 {
            (1..atoms.len() - 1).chain([0, atoms.len() - 1]).collect()
        } else {
            (start..atoms.len()).chain(0..start).collect()
        };
        for _ in 0..count {
            let mut found = false;
            for &i in &order {
                if self.free(atoms[i]) - used[i] >= need {
                    used[i] += need;
                    out.push(i);
                    found = true;
                    break;
                }
            }
            if !found {
                return Err("There isn't room on the chain for all of these groups.".into());
            }
        }
        Ok(out)
    }

    fn attach(&mut self, at: usize, sub: &Sub, n_atoms: &mut Vec<usize>) -> Result<(), String> {
        match sub {
            Sub::Halo(e) => {
                let x = self.add(*e);
                self.bond(at, x, 1);
            }
            Sub::Hydroxy => {
                let o = self.add(Elem::O);
                self.bond(at, o, 1);
            }
            Sub::Oxo => {
                let o = self.add(Elem::O);
                self.bond(at, o, 2);
            }
            Sub::Amino(subs) => {
                let n = self.add(Elem::N);
                self.bond(at, n, 1);
                for s in subs {
                    self.attach(n, s, n_atoms)?;
                }
            }
            Sub::Alkoxy(inner) => {
                let o = self.add(Elem::O);
                self.bond(at, o, 1);
                self.attach(o, inner, n_atoms)?;
            }
            Sub::Carboxy => {
                let c = self.add(Elem::C);
                self.bond(at, c, 1);
                let o1 = self.add(Elem::O);
                self.bond(c, o1, 2);
                let o2 = self.add(Elem::O);
                self.bond(c, o2, 1);
            }
            Sub::Cyano => {
                let c = self.add(Elem::C);
                self.bond(at, c, 1);
                let n = self.add(Elem::N);
                self.bond(c, n, 3);
            }
            Sub::Formyl => {
                let c = self.add(Elem::C);
                self.bond(at, c, 1);
                let o = self.add(Elem::O);
                self.bond(c, o, 2);
            }
            Sub::Carbamoyl(subs) => {
                let c = self.add(Elem::C);
                self.bond(at, c, 1);
                let o = self.add(Elem::O);
                self.bond(c, o, 2);
                let n = self.add(Elem::N);
                self.bond(c, n, 1);
                for s in subs {
                    self.attach(n, s, n_atoms)?;
                }
            }
            Sub::AlkoxyCarbonyl(inner) => {
                let c = self.add(Elem::C);
                self.bond(at, c, 1);
                let o = self.add(Elem::O);
                self.bond(c, o, 2);
                let o2 = self.add(Elem::O);
                self.bond(c, o2, 1);
                self.attach(o2, inner, n_atoms)?;
            }
            Sub::Carbon(f) => {
                let atoms = self.skeleton(f.len, f.ring, f.aromatic, &f.unsat)?;
                let ap = f.attach.unwrap_or(1);
                if ap < 1 || ap as usize > atoms.len() {
                    return Err(format!("A {}-carbon side group can't attach at position {ap}.", f.len));
                }
                self.bond(at, atoms[ap as usize - 1], f.attach_order);
                let mut local_n = vec![];
                self.apply_items(&atoms, &f.prefixes, &mut local_n, false)?;
            }
        }
        Ok(())
    }

    fn apply_items(&mut self, atoms: &[usize], items: &[Item], n_atoms: &mut Vec<usize>, is_chain_parent: bool) -> Result<(), String> {
        for it in items {
            match &it.locs {
                Some(locs) => {
                    if locs.len() != it.count {
                        return Err(format!(
                            "The numbers ({}) don't match the multiplier for one group (expected {} number{}).",
                            locs.iter().map(|l| match l { Loc::Num(n) => n.to_string(), Loc::N => "N".into() }).collect::<Vec<_>>().join(","),
                            it.count,
                            if it.count == 1 { "" } else { "s" }
                        ));
                    }
                    for l in locs {
                        let at = match l {
                            Loc::Num(n) => {
                                if *n < 1 || *n as usize > atoms.len() {
                                    return Err(format!(
                                        "Position {n} doesn't exist: the main chain/ring only has {} carbon{}.",
                                        atoms.len(),
                                        if atoms.len() == 1 { "" } else { "s" }
                                    ));
                                }
                                atoms[*n as usize - 1]
                            }
                            Loc::N => *n_atoms
                                .iter()
                                .find(|&&na| self.free(na) >= 1)
                                .ok_or("An \"N-\" group needs an amine or amide nitrogen to sit on.")?,
                        };
                        self.attach(at, &it.sub, n_atoms)?;
                    }
                }
                None => {
                    let need = it.sub.attach_order() as i32;
                    let inner = is_chain_parent && it.sub.is_carbon_branch();
                    let pos = self.default_positions(atoms, it.count, need, inner, 0)?;
                    for i in pos {
                        self.attach(atoms[i], &it.sub, n_atoms)?;
                    }
                }
            }
        }
        Ok(())
    }
}

// ---------------------------------------------------------------- public

fn normalize(input: &str) -> String {
    let mut s = input.trim().to_lowercase();
    for (a, b) in [("–", "-"), ("—", "-"), ("‐", "-"), ("[", "("), ("]", ")"), ("{", "("), ("}", ")"), ("′", "'")] {
        s = s.replace(a, b);
    }
    let mut out = String::new();
    let mut last_space = false;
    for ch in s.chars() {
        if ch.is_whitespace() {
            if !last_space {
                out.push(' ');
            }
            last_space = true;
        } else {
            out.push(ch);
            last_space = false;
        }
    }
    out = out.replace(" -", "-").replace("- ", "-").replace(", ", ",").replace(" ,", ",");
    // n-butane, n-propyl
    if !out.contains("amine") && !out.contains("amide") && !out.contains("aniline") {
        for pre in ["n-", "normal-"] {
            if let Some(r) = out.strip_prefix(pre) {
                out = r.to_string();
            }
        }
    }
    out
}

const EXTRA_COMMON: &[(&str, &str)] = &[
    ("isopropanol", "propan-2-ol"),
    ("rubbing alcohol", "propan-2-ol"),
    ("alcohol", "ethanol"),
    ("vinyl chloride", "chloroethene"),
    ("methyl ethyl ketone", "butan-2-one"),
    ("mek", "butan-2-one"),
    ("ether", "ethoxyethane"),
    ("glycol", "ethane-1,2-diol"),
    ("glycerine", "propane-1,2,3-triol"),
    ("glycerin", "propane-1,2,3-triol"),
    ("acetic acid", "ethanoic acid"),
    ("vinegar acid", "ethanoic acid"),
    ("acrylic acid", "prop-2-enoic acid"),
    ("propionaldehyde", "propanal"),
    ("butyraldehyde", "butanal"),
    ("methylamine", "methanamine"),
    ("ethylamine", "ethanamine"),
    ("dimethylamine", "n-methylmethanamine"),
    ("trimethylamine", "n,n-dimethylmethanamine"),
    ("xylene", "1,2-dimethylbenzene"),
    ("cumene", "isopropylbenzene"),
    ("isooctane", "2,2,4-trimethylpentane"),
    ("phenylethene", "ethenylbenzene"),
    ("methane", "methane"),
    ("acetylene", "ethyne"),
    ("ethylene", "ethene"),
    ("propylene", "propene"),
    ("isobutylene", "2-methylpropene"),
    ("catechol", "benzene-1,2-diol"),
    ("hydroquinone", "benzene-1,4-diol"),
    ("resorcinol", "benzene-1,3-diol"),
    ("nitrobenzene", "nitrobenzene"),
];

pub fn parse_name(input: &str) -> Result<Built, String> {
    let mut s = normalize(input);
    if s.is_empty() {
        return Err("Type a name first.".into());
    }
    if let Some(sys) = namer::common_to_systematic(&s) {
        s = sys.to_string();
    } else if let Some((_, sys)) = EXTRA_COMMON.iter().find(|x| x.0 == s) {
        s = sys.to_string();
    }
    if s.contains("nitro") {
        return Err("Nitro groups aren't supported yet.".into());
    }
    let interpreted = s.clone();

    // esters: "<alkyl words> <acid>oate"
    let mut ester_alkyls: Vec<Sub> = vec![];
    let mut body = s.as_str();
    if !body.ends_with(" acid") && body.contains(' ') {
        let idx = body.rfind(' ').unwrap();
        let (words, last) = (&body[..idx], &body[idx + 1..]);
        if !(last.ends_with("oate") || last.ends_with("carboxylate")) {
            return Err(format!(
                "I couldn't read \"{}\". Names are usually written as one word (except acids like \"ethanoic acid\" and esters like \"ethyl ethanoate\").",
                input.trim()
            ));
        }
        for w in words.split(' ') {
            let mut found = false;
            for (k, sub, r) in mult_sub_options(w) {
                if r.is_empty() {
                    for _ in 0..k {
                        ester_alkyls.push(sub.clone());
                    }
                    found = true;
                    break;
                }
            }
            if !found {
                // "1-methylethyl ethanoate"
                if let Some(sub) = complex_sub(w) {
                    ester_alkyls.push(sub);
                    found = true;
                }
            }
            if !found {
                return Err(format!("I couldn't read the ester group \"{w}\"."));
            }
        }
        body = last;
    }

    let (items, pending, mut parent) = parse_seq(body, &parent_full, 0).ok_or_else(|| {
        format!(
            "I couldn't read \"{}\". Check the spelling, e.g. \"2-methylbutane\", \"but-2-ene\" or \"2-butene\", \"propan-2-ol\" or \"2-propanol\".",
            input.trim()
        )
    })?;

    if let Some(p) = &pending {
        let p = nums(p).ok_or("An \"N\" can't number the main chain.")?;
        if !assign_pending_unsat(&mut parent.unsat, &p) {
            match &mut parent.suffix {
                Some(sf) if sf.locs.is_none() || (parent.aromatic && sf.count == 1) => {
                    sf.count = sf.count.max(p.len());
                    sf.locs = Some(p);
                }
                _ => return Err("There are numbers in the name that don't belong to anything.".into()),
            }
        }
    }
    if let Some(sf) = &parent.suffix {
        if sf.kind == Suffix::Ester && ester_alkyls.is_empty() {
            return Err("An ester needs the group on the oxygen as a first word, e.g. \"methyl ethanoate\".".into());
        }
        if sf.kind != Suffix::Ester && !ester_alkyls.is_empty() {
            return Err("Only esters (ending in \"-oate\") are written as two words.".into());
        }
    }

    let mut b = Builder { mol: Molecule::default() };
    let atoms = b.skeleton(parent.len, parent.ring, parent.aromatic, &parent.unsat)?;
    let n = atoms.len();
    let mut n_atoms: Vec<usize> = vec![];
    if let Some(sf) = &parent.suffix {
        let terminal_chain = sf.kind.is_terminal() && !sf.carbo && !parent.ring;
        if sf.kind.is_terminal() && !sf.carbo && parent.ring {
            return Err(format!(
                "On a ring this group is written with \"{}\", e.g. \"cyclohexane{}\".",
                sf.kind.carbo_word(),
                sf.kind.carbo_word()
            ));
        }
        let positions: Vec<usize> = match &sf.locs {
            Some(ls) => {
                let mut v = vec![];
                for &l in ls {
                    if l < 1 || l as usize > n {
                        return Err(format!("Position {l} doesn't exist on a {n}-carbon {}.", if parent.ring { "ring" } else { "chain" }));
                    }
                    v.push(l as usize - 1);
                }
                v
            }
            None => {
                if terminal_chain {
                    match sf.count {
                        1 => vec![0],
                        2 if n > 1 => vec![0, n - 1],
                        _ => return Err("A chain only has two ends for this group.".into()),
                    }
                } else if sf.kind == Suffix::Ketone && !parent.ring {
                    if n < 3 {
                        return Err("A ketone needs at least 3 carbons (the C=O can't be at the end).".into());
                    }
                    b.default_positions(&atoms[1..n - 1], sf.count, 2, false, 0)?.into_iter().map(|i| i + 1).collect()
                } else {
                    let need = if sf.kind == Suffix::Ketone { 2 } else { 1 };
                    b.default_positions(&atoms, sf.count, need, false, 0)?
                }
            }
        };
        if positions.len() != sf.count {
            return Err("The numbers don't match the multiplier on the ending.".into());
        }
        for (k, &pi) in positions.iter().enumerate() {
            let at = atoms[pi];
            let c = if sf.carbo {
                let c = b.add(Elem::C);
                b.bond(at, c, 1);
                c
            } else {
                at
            };
            match sf.kind {
                Suffix::Acid => {
                    let o1 = b.add(Elem::O);
                    b.bond(c, o1, 2);
                    let o2 = b.add(Elem::O);
                    b.bond(c, o2, 1);
                }
                Suffix::Ester => {
                    let o1 = b.add(Elem::O);
                    b.bond(c, o1, 2);
                    let o2 = b.add(Elem::O);
                    b.bond(c, o2, 1);
                    let alk = ester_alkyls.get(k).or(ester_alkyls.first()).cloned().unwrap();
                    b.attach(o2, &alk, &mut n_atoms)?;
                }
                Suffix::Amide => {
                    let o1 = b.add(Elem::O);
                    b.bond(c, o1, 2);
                    let nn = b.add(Elem::N);
                    b.bond(c, nn, 1);
                    n_atoms.push(nn);
                }
                Suffix::Nitrile => {
                    let nn = b.add(Elem::N);
                    b.bond(c, nn, 3);
                }
                Suffix::Aldehyde | Suffix::Ketone => {
                    let o1 = b.add(Elem::O);
                    b.bond(c, o1, 2);
                }
                Suffix::Alcohol => {
                    let o1 = b.add(Elem::O);
                    b.bond(c, o1, 1);
                }
                Suffix::Amine => {
                    let nn = b.add(Elem::N);
                    b.bond(c, nn, 1);
                    n_atoms.push(nn);
                }
            }
        }
    }
    b.apply_items(&atoms, &items, &mut n_atoms, !parent.ring)?;

    // valence check
    for i in 0..b.mol.atoms.len() {
        if b.mol.over_valent(i) {
            let pos = atoms.iter().position(|&a| a == i);
            return Err(match pos {
                Some(p) => format!(
                    "This name puts too many bonds on carbon {} of the main {} (carbon can only make 4).",
                    p + 1,
                    if parent.ring { "ring" } else { "chain" }
                ),
                None => "This name puts too many bonds on one atom (carbon can only make 4).".into(),
            });
        }
    }
    Ok(Built { mol: b.mol, backbone: atoms, backbone_ring: parent.ring, interpreted })
}
