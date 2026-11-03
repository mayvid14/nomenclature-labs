//! Structure → IUPAC name.
//!
//! Covers acyclic and monocyclic carbon skeletons (including benzene) with
//! halogens, alcohols, ethers, aldehydes, ketones, carboxylic acids, esters,
//! amides, nitriles and amines, following the IUPAC 2013 recommendations for
//! choosing and numbering the parent chain.

use crate::molecule::{Elem, Molecule};
use std::cell::RefCell;
use std::collections::HashMap;

/// Naming convention.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Style {
    /// Textbook style: "2-butene", "2-propanol"; the chain with most multiple bonds wins.
    Classic,
    /// IUPAC 2013: "but-2-ene", "propan-2-ol"; the longest chain wins.
    Modern,
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Hash)]
pub enum Suffix {
    Acid,
    Ester,
    Amide,
    Nitrile,
    Aldehyde,
    Ketone,
    Alcohol,
    Amine,
}

impl Suffix {
    pub fn family(self) -> &'static str {
        match self {
            Suffix::Acid => "Carboxylic acid",
            Suffix::Ester => "Ester",
            Suffix::Amide => "Amide",
            Suffix::Nitrile => "Nitrile",
            Suffix::Aldehyde => "Aldehyde",
            Suffix::Ketone => "Ketone",
            Suffix::Alcohol => "Alcohol",
            Suffix::Amine => "Amine",
        }
    }
    pub fn group_desc(self) -> &'static str {
        match self {
            Suffix::Acid => "carboxyl group (–COOH)",
            Suffix::Ester => "ester group (–COO–)",
            Suffix::Amide => "amide group (–CONH₂)",
            Suffix::Nitrile => "nitrile group (–C≡N)",
            Suffix::Aldehyde => "aldehyde group (–CHO)",
            Suffix::Ketone => "carbonyl group (C=O)",
            Suffix::Alcohol => "hydroxyl group (–OH)",
            Suffix::Amine => "amino group (–NH₂)",
        }
    }
    /// Groups whose carbon ends a chain (and is counted in it).
    pub fn is_terminal(self) -> bool {
        matches!(
            self,
            Suffix::Acid | Suffix::Ester | Suffix::Amide | Suffix::Nitrile | Suffix::Aldehyde
        )
    }
    pub fn chain_word(self) -> &'static str {
        match self {
            Suffix::Acid => "oic acid",
            Suffix::Ester => "oate",
            Suffix::Amide => "amide",
            Suffix::Nitrile => "nitrile",
            Suffix::Aldehyde => "al",
            Suffix::Ketone => "one",
            Suffix::Alcohol => "ol",
            Suffix::Amine => "amine",
        }
    }
    pub fn carbo_word(self) -> &'static str {
        match self {
            Suffix::Acid => "carboxylic acid",
            Suffix::Ester => "carboxylate",
            Suffix::Amide => "carboxamide",
            Suffix::Nitrile => "carbonitrile",
            Suffix::Aldehyde => "carbaldehyde",
            other => other.chain_word(),
        }
    }
    fn benzene_retained(self) -> Option<&'static str> {
        match self {
            Suffix::Acid => Some("benzoic acid"),
            Suffix::Ester => Some("benzoate"),
            Suffix::Amide => Some("benzamide"),
            Suffix::Nitrile => Some("benzonitrile"),
            Suffix::Aldehyde => Some("benzaldehyde"),
            Suffix::Alcohol => Some("phenol"),
            Suffix::Amine => Some("aniline"),
            Suffix::Ketone => None,
        }
    }
}

#[derive(Clone, Debug)]
pub struct NameResult {
    pub name: String,
    /// Original atom indices of the parent chain/ring, in locant order (locant = index + 1).
    pub parent_atoms: Vec<usize>,
    pub parent_is_ring: bool,
    pub family: String,
    pub features: Vec<String>,
    pub steps: Vec<String>,
    pub common: Option<&'static str>,
}

#[derive(Clone, Debug)]
pub struct NameError {
    pub message: String,
    pub atoms: Vec<usize>,
}

fn err<T>(msg: impl Into<String>) -> Result<T, NameError> {
    Err(NameError { message: msg.into(), atoms: vec![] })
}

// ---------------------------------------------------------------- words

const STEMS: [&str; 20] = [
    "meth", "eth", "prop", "but", "pent", "hex", "hept", "oct", "non", "dec", "undec", "dodec",
    "tridec", "tetradec", "pentadec", "hexadec", "heptadec", "octadec", "nonadec", "icos",
];

pub fn stem(n: usize) -> String {
    if (1..=20).contains(&n) {
        STEMS[n - 1].to_string()
    } else {
        format!("C{n}-")
    }
}

pub fn stems() -> &'static [&'static str] {
    &STEMS
}

fn mult_simple(k: usize) -> &'static str {
    ["", "", "di", "tri", "tetra", "penta", "hexa", "hepta", "octa", "nona", "deca"]
        .get(k)
        .copied()
        .unwrap_or("poly")
}

fn mult_complex(k: usize) -> &'static str {
    ["", "", "bis", "tris", "tetrakis", "pentakis", "hexakis"]
        .get(k)
        .copied()
        .unwrap_or("polykis")
}

fn join_nums(v: &[u32]) -> String {
    v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")
}

/// "propane", "prop-1-ene", "buta-1,3-diene", "pent-1-en-4-yne".
pub fn hydride_name(stem: &str, doubles: &[u32], triples: &[u32], show_locs: bool) -> String {
    if doubles.is_empty() && triples.is_empty() {
        return format!("{stem}ane");
    }
    let mut s = stem.to_string();
    let first_mult = if !doubles.is_empty() { doubles.len() } else { triples.len() };
    if first_mult > 1 {
        s.push('a');
    }
    if !doubles.is_empty() {
        if show_locs {
            s += &format!("-{}-", join_nums(doubles));
        }
        s += mult_simple(doubles.len());
        s += "ene";
    }
    if !triples.is_empty() {
        if !doubles.is_empty() {
            s.pop(); // "ene" -> "en" before "yne"
        }
        if show_locs {
            s += &format!("-{}-", join_nums(triples));
        }
        s += mult_simple(triples.len());
        s += "yne";
    }
    s
}

fn starts_with_vowel(s: &str) -> bool {
    matches!(s.chars().next(), Some('a' | 'e' | 'i' | 'o' | 'u' | 'y'))
}

fn attach_suffix(hydride: &str, locs: Option<String>, count: usize, word: &str) -> String {
    let mword = format!("{}{}", mult_simple(count), word);
    let mut base = hydride.to_string();
    if base.ends_with('e') && starts_with_vowel(&mword) {
        base.pop();
    }
    match locs {
        Some(l) => format!("{base}-{l}-{mword}"),
        None => format!("{base}{mword}"),
    }
}

// ---------------------------------------------------------------- substituents

#[derive(Clone, Debug, PartialEq)]
struct Sub {
    text: String,
    /// Needs parentheses when cited (has locants or inner prefixes).
    paren: bool,
    /// Is itself substituted (uses bis/tris when multiplied).
    substituted: bool,
}

fn simple(t: &str) -> Sub {
    Sub { text: t.to_string(), paren: false, substituted: false }
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
enum Loc {
    N,
    Num(u32),
}

#[derive(Clone, Debug)]
struct Prefix {
    loc: Loc,
    sub: Sub,
}

fn alpha_key(s: &Sub) -> String {
    let t = s.text.trim_start_matches("tert-").trim_start_matches("sec-");
    t.chars().filter(|c| c.is_ascii_alphabetic()).collect::<String>().to_lowercase()
}

/// Groups identical prefixes, multiplies them and orders alphabetically.
fn assemble_prefixes(prefixes: &[Prefix], show_locs: bool) -> String {
    let mut groups: Vec<(Sub, Vec<Loc>)> = vec![];
    for p in prefixes {
        if let Some(g) = groups.iter_mut().find(|g| g.0.text == p.sub.text) {
            g.1.push(p.loc);
        } else {
            groups.push((p.sub.clone(), vec![p.loc]));
        }
    }
    for g in &mut groups {
        g.1.sort();
    }
    groups.sort_by(|a, b| alpha_key(&a.0).cmp(&alpha_key(&b.0)).then(a.1.cmp(&b.1)));
    let mut parts: Vec<String> = vec![];
    for (sub, locs) in &groups {
        let k = locs.len();
        let mult = if k > 1 {
            if sub.substituted { mult_complex(k) } else { mult_simple(k) }
        } else {
            ""
        };
        let name = if sub.paren || (k > 1 && sub.substituted) {
            format!("({})", sub.text)
        } else {
            sub.text.clone()
        };
        let any_n = locs.iter().any(|l| *l == Loc::N);
        if show_locs || any_n {
            let ls: Vec<String> = locs
                .iter()
                .map(|l| match l {
                    Loc::N => "N".to_string(),
                    Loc::Num(n) => n.to_string(),
                })
                .collect();
            parts.push(format!("{}-{}{}", ls.join(","), mult, name));
        } else {
            parts.push(format!("{mult}{name}"));
        }
    }
    let mut out = String::new();
    for p in parts {
        if !out.is_empty() && p.starts_with(|c: char| c.is_ascii_digit() || c == 'N') {
            out.push('-');
        }
        out += &p;
    }
    out
}

fn alkoxy(inner: &Sub) -> Sub {
    let t = &inner.text;
    let short = match t.as_str() {
        "methyl" => Some("methoxy"),
        "ethyl" => Some("ethoxy"),
        "propyl" => Some("propoxy"),
        "butyl" => Some("butoxy"),
        "phenyl" => Some("phenoxy"),
        "tert-butyl" => Some("tert-butoxy"),
        _ => None,
    };
    if let Some(s) = short {
        return simple(s);
    }
    for st in ["methyl", "ethyl", "propyl", "butyl"] {
        if t.ends_with(st) && inner.substituted {
            let text = format!("{}oxy", &t[..t.len() - 2]);
            return Sub { text, paren: true, substituted: true };
        }
    }
    let paren = inner.paren || inner.substituted;
    let text = if inner.paren { format!("{t}oxy") } else { format!("{t}oxy") };
    Sub { text, paren, substituted: inner.substituted }
}

// ---------------------------------------------------------------- graph

#[derive(Clone, Copy, PartialEq, Debug)]
enum CType {
    Plain,
    Group(Suffix),
}

struct G {
    el: Vec<Elem>,
    adj: Vec<Vec<(usize, u8)>>,
    orig: Vec<usize>,
}

impl G {
    fn order(&self, a: usize, b: usize) -> u8 {
        self.adj[a].iter().find(|x| x.0 == b).map(|x| x.1).unwrap_or(0)
    }
    fn deg(&self, a: usize) -> usize {
        self.adj[a].len()
    }
    fn carbons_of(&self, a: usize) -> Vec<usize> {
        self.adj[a].iter().filter(|x| self.el[x.0] == Elem::C).map(|x| x.0).collect()
    }
}

/// Simple rings of a graph given as adjacency. Err(atoms) for fused/spiro systems.
fn rings_of(adj: &[Vec<(usize, u8)>]) -> Result<Vec<Vec<usize>>, Vec<usize>> {
    let n = adj.len();
    let mut visited = vec![false; n];
    let mut parent = vec![usize::MAX; n];
    let mut depth = vec![0usize; n];
    let mut rings: Vec<Vec<usize>> = vec![];
    fn dfs(
        u: usize,
        adj: &[Vec<(usize, u8)>],
        visited: &mut Vec<bool>,
        parent: &mut Vec<usize>,
        depth: &mut Vec<usize>,
        rings: &mut Vec<Vec<usize>>,
    ) {
        visited[u] = true;
        for &(v, _) in &adj[u] {
            if v == parent[u] {
                continue;
            }
            if !visited[v] {
                parent[v] = u;
                depth[v] = depth[u] + 1;
                dfs(v, adj, visited, parent, depth, rings);
            } else if depth[v] < depth[u] {
                let mut cyc = vec![u];
                let mut x = u;
                while x != v {
                    x = parent[x];
                    cyc.push(x);
                }
                rings.push(cyc);
            }
        }
    }
    for s in 0..n {
        if !visited[s] {
            dfs(s, adj, &mut visited, &mut parent, &mut depth, &mut rings);
        }
    }
    let mut count = vec![0; n];
    for r in &rings {
        for &a in r {
            count[a] += 1;
        }
    }
    let shared: Vec<usize> = (0..n).filter(|&i| count[i] > 1).collect();
    if !shared.is_empty() {
        return Err(shared);
    }
    Ok(rings)
}

/// Rings of a whole molecule (atom indices in cyclic order). Fused systems give partial results.
pub fn molecule_rings(m: &Molecule) -> Vec<Vec<usize>> {
    let mut adj = vec![vec![]; m.atoms.len()];
    for b in &m.bonds {
        adj[b.a].push((b.b, b.order));
        adj[b.b].push((b.a, b.order));
    }
    match rings_of(&adj) {
        Ok(r) => r,
        Err(_) => {
            // still return the fundamental cycles for drawing purposes
            let mut visited = vec![false; adj.len()];
            let mut out = vec![];
            let _ = &mut visited;
            out.clear();
            out
        }
    }
}

// ---------------------------------------------------------------- namer

#[derive(PartialEq, Eq, PartialOrd, Ord, Debug, Clone)]
struct Key {
    neg_suffix: i32,
    k1: i32,
    k2: i32,
    neg_dbl: i32,
    suffix_locs: Vec<u32>,
    mult_locs: Vec<u32>,
    dbl_locs: Vec<u32>,
    neg_nprefix: i32,
    prefix_locs: Vec<u32>,
    alpha_locs: Vec<u32>,
    text: String,
}

struct Eval {
    key: Key,
    atoms: Vec<usize>,
    is_ring: bool,
    aromatic: bool,
    suffix_locs: Vec<u32>,
    carbo: bool,
    doubles: Vec<u32>,
    triples: Vec<u32>,
    prefixes: Vec<Prefix>,
    ester_alkyls: Vec<Sub>,
}

struct Namer<'a> {
    g: &'a G,
    style: Style,
    ring_id: Vec<Option<usize>>,
    rings: Vec<Vec<usize>>,
    aromatic: Vec<bool>,
    ctype: Vec<CType>,
    principal: Option<Suffix>,
    memo: RefCell<HashMap<(usize, usize), Result<Sub, String>>>,
}

impl<'a> Namer<'a> {
    fn new(g: &'a G, style: Style) -> Result<Self, NameError> {
        let n = g.el.len();
        let rings = match rings_of(&g.adj) {
            Ok(r) => r,
            Err(atoms) => {
                return Err(NameError {
                    message: "Fused, bridged or spiro ring systems (rings sharing atoms) aren't supported yet. Try separate rings joined by a bond.".into(),
                    atoms: atoms.iter().map(|&a| g.orig[a]).collect(),
                });
            }
        };
        let mut ring_id = vec![None; n];
        for (ri, r) in rings.iter().enumerate() {
            for &a in r {
                if g.el[a] != Elem::C {
                    return Err(NameError {
                        message: "Rings containing O or N (heterocycles) aren't supported yet.".into(),
                        atoms: vec![g.orig[a]],
                    });
                }
                ring_id[a] = Some(ri);
            }
        }
        let aromatic: Vec<bool> = rings
            .iter()
            .map(|r| {
                if r.len() != 6 {
                    return false;
                }
                let o: Vec<u8> = (0..6).map(|i| g.order(r[i], r[(i + 1) % 6])).collect();
                (o[0] == 2 && o[1] == 1 && o[2] == 2 && o[3] == 1 && o[4] == 2 && o[5] == 1)
                    || (o[0] == 1 && o[1] == 2 && o[2] == 1 && o[3] == 2 && o[4] == 1 && o[5] == 2)
            })
            .collect();

        // validate heteroatoms
        for a in 0..n {
            let e = g.el[a];
            let bad = |msg: &str| -> Result<Self, NameError> {
                Err(NameError { message: msg.to_string(), atoms: vec![g.orig[a]] })
            };
            if e == Elem::C {
                continue;
            }
            for &(b, o) in &g.adj[a] {
                if g.el[b] != Elem::C {
                    return bad(&format!(
                        "A {}–{} bond isn't part of the groups this app can name (only bonds to carbon are supported).",
                        e.symbol(),
                        g.el[b].symbol()
                    ));
                }
                if e.is_halogen() && o != 1 {
                    return bad("Halogens can only make one single bond.");
                }
            }
            if e == Elem::O {
                let d = g.deg(a);
                if d == 0 {
                    return bad("A lone oxygen is water, not an organic compound.");
                }
                if d == 2 && g.adj[a].iter().any(|x| x.1 != 1) {
                    return bad("This oxygen has too many bonds.");
                }
            }
            if e == Elem::N {
                let orders: Vec<u8> = g.adj[a].iter().map(|x| x.1).collect();
                if orders.contains(&2) {
                    return bad("C=N double bonds (imines) aren't supported yet.");
                }
                if orders.contains(&3) && orders.len() > 1 {
                    return bad("This nitrogen has too many bonds.");
                }
            }
        }

        // classify carbons
        let mut ctype = vec![CType::Plain; n];
        for c in 0..n {
            if g.el[c] != Elem::C {
                continue;
            }
            let mut dbl_o = false;
            let mut oh = false;
            let mut o_link = false;
            let mut n_single = false;
            let mut n_triple = false;
            let mut halo = false;
            for &(b, o) in &g.adj[c] {
                match g.el[b] {
                    Elem::O if o == 2 => dbl_o = true,
                    Elem::O if g.deg(b) == 1 => oh = true,
                    Elem::O => o_link = true,
                    Elem::N if o == 3 => n_triple = true,
                    Elem::N => n_single = true,
                    e if e.is_halogen() => halo = true,
                    _ => {}
                }
            }
            let nc = g.carbons_of(c).len();
            ctype[c] = if n_triple {
                CType::Group(Suffix::Nitrile)
            } else if dbl_o {
                if oh {
                    CType::Group(Suffix::Acid)
                } else if o_link {
                    CType::Group(Suffix::Ester)
                } else if n_single {
                    CType::Group(Suffix::Amide)
                } else if halo {
                    return Err(NameError {
                        message: "Acid halides (–COCl and similar) aren't supported yet.".into(),
                        atoms: vec![g.orig[c]],
                    });
                } else if nc <= 1 {
                    CType::Group(Suffix::Aldehyde)
                } else {
                    CType::Group(Suffix::Ketone)
                }
            } else {
                CType::Plain
            };
        }
        // anhydrides / ester checks
        for a in 0..n {
            if g.el[a] == Elem::O && g.deg(a) == 2 {
                let carbonyls = g.adj[a]
                    .iter()
                    .filter(|x| ctype[x.0] == CType::Group(Suffix::Ester))
                    .count();
                if carbonyls > 1 {
                    return Err(NameError {
                        message: "Acid anhydrides aren't supported yet.".into(),
                        atoms: vec![g.orig[a]],
                    });
                }
            }
        }

        let mut namer = Namer {
            g,
            style,
            ring_id,
            rings,
            aromatic,
            ctype,
            principal: None,
            memo: RefCell::new(HashMap::new()),
        };
        // principal group
        let mut best: Option<Suffix> = None;
        for a in 0..n {
            let mut cand = vec![];
            if let CType::Group(s) = namer.ctype[a] {
                cand.push(s);
            }
            if g.el[a] == Elem::O && g.deg(a) == 1 && g.adj[a][0].1 == 1 {
                let c = g.adj[a][0].0;
                if namer.ctype[c] == CType::Plain || namer.ctype[c] == CType::Group(Suffix::Ketone) {
                    cand.push(Suffix::Alcohol);
                }
            }
            if namer.is_amine_n(a) {
                cand.push(Suffix::Amine);
            }
            for s in cand {
                best = Some(best.map_or(s, |b: Suffix| b.min(s)));
            }
        }
        namer.principal = best;
        Ok(namer)
    }

    fn is_amine_n(&self, a: usize) -> bool {
        let g = self.g;
        g.el[a] == Elem::N
            && g.adj[a].iter().all(|x| x.1 == 1)
            && g.adj[a]
                .iter()
                .all(|x| self.ctype[x.0] != CType::Group(Suffix::Amide))
            && !g.adj[a].is_empty()
    }

    fn ester_link(&self, c: usize) -> Option<usize> {
        self.g.adj[c]
            .iter()
            .find(|x| self.g.el[x.0] == Elem::O && x.1 == 1 && self.g.deg(x.0) == 2)
            .map(|x| x.0)
    }

    /// Heteroatoms that belong to the functional group on carbon c.
    fn group_atoms(&self, c: usize) -> Vec<usize> {
        let g = self.g;
        let CType::Group(s) = self.ctype[c] else { return vec![] };
        g.adj[c]
            .iter()
            .filter(|&&(b, o)| match s {
                Suffix::Acid => g.el[b] == Elem::O && g.deg(b) == 1,
                Suffix::Ester => {
                    g.el[b] == Elem::O && (o == 2 || Some(b) == self.ester_link(c))
                }
                Suffix::Amide => (g.el[b] == Elem::O && o == 2) || g.el[b] == Elem::N,
                Suffix::Nitrile => g.el[b] == Elem::N,
                Suffix::Aldehyde | Suffix::Ketone => g.el[b] == Elem::O && o == 2,
                _ => false,
            })
            .map(|x| x.0)
            .collect()
    }

    // ------------------------------------------------ substituent names

    fn sub(&self, a: usize, b: usize) -> Result<Sub, String> {
        if let Some(r) = self.memo.borrow().get(&(a, b)) {
            return r.clone();
        }
        let r = self.sub_inner(a, b);
        self.memo.borrow_mut().insert((a, b), r.clone());
        r
    }

    fn n_subs(&self, n_atom: usize, except: usize) -> Result<Vec<Sub>, String> {
        let mut v = vec![];
        for &(x, _) in &self.g.adj[n_atom] {
            if x == except {
                continue;
            }
            if let CType::Group(Suffix::Amide | Suffix::Ester | Suffix::Acid) = self.ctype[x] {
                if self.group_atoms(x).contains(&n_atom) && self.ctype[x] == CType::Group(Suffix::Amide) {
                    return Err("Amides joined through nitrogen (acylamino groups) aren't supported yet.".into());
                }
            }
            v.push(self.sub(n_atom, x)?);
        }
        Ok(v)
    }

    fn grouped_text(subs: &[Sub]) -> String {
        let prefixes: Vec<Prefix> =
            subs.iter().map(|s| Prefix { loc: Loc::Num(1), sub: s.clone() }).collect();
        // multiply without locants; wrap complex names
        let mut groups: Vec<(Sub, usize)> = vec![];
        for p in &prefixes {
            if let Some(g) = groups.iter_mut().find(|g| g.0.text == p.sub.text) {
                g.1 += 1;
            } else {
                groups.push((p.sub.clone(), 1));
            }
        }
        groups.sort_by(|a, b| alpha_key(&a.0).cmp(&alpha_key(&b.0)));
        let n_groups = groups.len();
        let mut out = String::new();
        for (i, (s, k)) in groups.iter().enumerate() {
            let m = if *k > 1 { mult_simple(*k) } else { "" };
            let wrap = s.paren || (n_groups > 1 && i + 1 < n_groups);
            if wrap {
                out += &format!("{m}({})", s.text);
            } else {
                out += &format!("{m}{}", s.text);
            }
        }
        out
    }

    fn sub_inner(&self, a: usize, b: usize) -> Result<Sub, String> {
        let g = self.g;
        let e = g.el[b];
        let ord = g.order(a, b);
        if e.is_halogen() {
            return Ok(simple(e.halo_prefix()));
        }
        match e {
            Elem::O => {
                if ord == 2 {
                    return Ok(simple("oxo"));
                }
                if g.deg(b) == 1 {
                    return Ok(simple("hydroxy"));
                }
                let c2 = g.adj[b].iter().find(|x| x.0 != a).map(|x| x.0).unwrap();
                if self.ctype[c2] == CType::Group(Suffix::Ester) && self.ester_link(c2) == Some(b) {
                    return Err("Acyloxy side groups (–O–CO–R) aren't supported yet.".into());
                }
                let inner = self.sub(b, c2)?;
                Ok(alkoxy(&inner))
            }
            Elem::N => {
                if ord == 3 {
                    return Err("Unexpected C≡N arrangement.".into());
                }
                let subs = self.n_subs(b, a)?;
                if subs.is_empty() {
                    Ok(simple("amino"))
                } else {
                    let t = Self::grouped_text(&subs);
                    Ok(Sub { text: format!("{t}amino"), paren: true, substituted: true })
                }
            }
            Elem::C => {
                match self.ctype[b] {
                    CType::Group(Suffix::Acid) => return Ok(simple("carboxy")),
                    CType::Group(Suffix::Nitrile) => return Ok(simple("cyano")),
                    CType::Group(Suffix::Amide) => {
                        let n_atom = g.adj[b].iter().find(|x| g.el[x.0] == Elem::N).unwrap().0;
                        let subs = self.n_subs(n_atom, b)?;
                        if subs.is_empty() {
                            return Ok(simple("carbamoyl"));
                        }
                        let t = Self::grouped_text(&subs);
                        return Ok(Sub {
                            text: format!("{t}carbamoyl"),
                            paren: true,
                            substituted: true,
                        });
                    }
                    CType::Group(Suffix::Ester) => {
                        if let Some(o) = self.ester_link(b) {
                            if o != a {
                                let c2 = g.adj[o].iter().find(|x| x.0 != b).unwrap().0;
                                let inner = self.sub(o, c2)?;
                                let ax = alkoxy(&inner);
                                return Ok(Sub {
                                    text: format!("{}carbonyl", ax.text),
                                    paren: true,
                                    substituted: true,
                                });
                            }
                        }
                        return Err("This ester arrangement isn't supported yet.".into());
                    }
                    CType::Group(Suffix::Aldehyde) if g.carbons_of(b).iter().all(|&x| x == a) => {
                        return Ok(simple("formyl"));
                    }
                    _ => {}
                }
                if self.ring_id[b].is_some() {
                    self.ring_sub(a, b)
                } else {
                    self.alkyl_sub(a, b)
                }
            }
            _ => Err(format!("Element {} isn't supported.", e.symbol())),
        }
    }

    fn chain_eligible_sub(&self, x: usize) -> bool {
        self.g.el[x] == Elem::C
            && self.ring_id[x].is_none()
            && matches!(
                self.ctype[x],
                CType::Plain | CType::Group(Suffix::Ketone) | CType::Group(Suffix::Aldehyde)
            )
    }

    fn alkyl_sub(&self, a: usize, b: usize) -> Result<Sub, String> {
        let g = self.g;
        // all simple paths starting at b, not going back through a
        let mut paths: Vec<Vec<usize>> = vec![];
        let mut stack: Vec<Vec<usize>> = vec![vec![b]];
        while let Some(p) = stack.pop() {
            let last = *p.last().unwrap();
            for &(x, _) in &g.adj[last] {
                if x == a || p.contains(&x) || !self.chain_eligible_sub(x) {
                    continue;
                }
                let mut q = p.clone();
                q.push(x);
                stack.push(q);
            }
            paths.push(p);
        }
        let mut best: Option<(Key, Sub)> = None;
        for p in paths {
            let mut prefixes = vec![];
            for (i, &at) in p.iter().enumerate() {
                for &(x, _) in &g.adj[at] {
                    if p.contains(&x) || (i == 0 && x == a) {
                        continue;
                    }
                    prefixes.push(Prefix { loc: Loc::Num(i as u32 + 1), sub: self.sub(at, x)? });
                }
            }
            let (doubles, triples) = self.chain_unsat(&p, false);
            let key = self.make_key(0, true, p.len(), &[], &doubles, &triples, &prefixes);
            let text = self.alkyl_text(p.len(), &doubles, &triples, &prefixes, g.order(a, b));
            let key = Key { text: text.text.clone(), ..key };
            if best.as_ref().is_none_or(|bk| key < bk.0) {
                best = Some((key, text));
            }
        }
        Ok(best.unwrap().1)
    }

    fn alkyl_text(
        &self,
        len: usize,
        doubles: &[u32],
        triples: &[u32],
        prefixes: &[Prefix],
        attach_order: u8,
    ) -> Sub {
        let classic = self.style == Style::Classic;
        let ending = match attach_order {
            2 => "ylidene",
            3 => "ylidyne",
            _ => "yl",
        };
        let base = if doubles.is_empty() && triples.is_empty() {
            format!("{}{ending}", stem(len))
        } else if len == 2 && doubles == [1] && attach_order == 1 {
            "ethenyl".into()
        } else if len == 2 && triples == [1] {
            "ethynyl".into()
        } else if classic {
            // "2-propenyl", "1,3-butadienyl"
            let mut h = hydride_name(&stem(len), doubles, triples, false);
            h.pop();
            let mut locs: Vec<u32> = doubles.iter().chain(triples.iter()).copied().collect();
            locs.sort();
            format!("{}-{h}{ending}", join_nums(&locs))
        } else {
            let mut h = hydride_name(&stem(len), doubles, triples, true);
            h.pop();
            format!("{h}-1-{ending}")
        };
        let pre = assemble_prefixes(prefixes, true);
        let text = format!("{pre}{base}");
        let mapped: Option<(&str, bool)> = match (text.as_str(), classic) {
            ("1-methylethyl", true) => Some(("isopropyl", false)),
            ("1-methylethyl", false) => Some(("propan-2-yl", true)),
            ("1,1-dimethylethyl", _) => Some(("tert-butyl", false)),
            ("1-methylpropyl", true) => Some(("sec-butyl", false)),
            ("1-methylpropyl", false) => Some(("butan-2-yl", true)),
            ("2-methylpropyl", true) => Some(("isobutyl", false)),
            ("1-oxoethyl", _) => Some(("acetyl", false)),
            ("phenylmethyl", _) => Some(("benzyl", false)),
            _ => None,
        };
        if let Some((t, paren)) = mapped {
            return Sub { text: t.into(), paren, substituted: false };
        }
        let has_digit = base.chars().any(|c| c.is_ascii_digit());
        Sub { text, paren: !prefixes.is_empty() || has_digit, substituted: !prefixes.is_empty() }
    }

    fn ring_sub(&self, a: usize, b: usize) -> Result<Sub, String> {
        let g = self.g;
        let ri = self.ring_id[b].unwrap();
        let ring = &self.rings[ri];
        let n = ring.len();
        let k = ring.iter().position(|&x| x == b).unwrap();
        let arom = self.aromatic[ri];
        let mut best: Option<(Key, Sub)> = None;
        for dir in [1isize, -1] {
            let atoms: Vec<usize> = (0..n)
                .map(|i| ring[((k as isize + dir * i as isize).rem_euclid(n as isize)) as usize])
                .collect();
            let mut prefixes = vec![];
            for (i, &at) in atoms.iter().enumerate() {
                for &(x, _) in &g.adj[at] {
                    if atoms.contains(&x) || (i == 0 && x == a) {
                        continue;
                    }
                    prefixes.push(Prefix { loc: Loc::Num(i as u32 + 1), sub: self.sub(at, x)? });
                }
            }
            let (doubles, triples) = if arom { (vec![], vec![]) } else { self.chain_unsat(&atoms, true) };
            let key = self.make_key(0, true, n, &[], &doubles, &triples, &prefixes);
            let ending = if g.order(a, b) == 2 { "ylidene" } else { "yl" };
            let base = if arom {
                "phenyl".to_string()
            } else if doubles.is_empty() && triples.is_empty() {
                format!("cyclo{}{ending}", stem(n))
            } else if self.style == Style::Classic {
                let mut h = hydride_name(&stem(n), &doubles, &triples, false);
                h.pop();
                let mut locs: Vec<u32> = doubles.iter().chain(triples.iter()).copied().collect();
                locs.sort();
                format!("{}-cyclo{h}{ending}", join_nums(&locs))
            } else {
                let mut h = hydride_name(&stem(n), &doubles, &triples, true);
                h.pop();
                format!("cyclo{h}-1-{ending}")
            };
            let pre = assemble_prefixes(&prefixes, true);
            let has_digit = base.chars().any(|c| c.is_ascii_digit());
            let sub = Sub {
                text: format!("{pre}{base}"),
                paren: !prefixes.is_empty() || has_digit,
                substituted: !prefixes.is_empty(),
            };
            let key = Key { text: sub.text.clone(), ..key };
            if best.as_ref().is_none_or(|bk| key < bk.0) {
                best = Some((key, sub));
            }
        }
        Ok(best.unwrap().1)
    }

    fn chain_unsat(&self, atoms: &[usize], ring: bool) -> (Vec<u32>, Vec<u32>) {
        let n = atoms.len();
        let mut d = vec![];
        let mut t = vec![];
        let m = if ring { n } else { n.saturating_sub(1) };
        for i in 0..m {
            let o = self.g.order(atoms[i], atoms[(i + 1) % n]);
            if o == 2 {
                d.push(i as u32 + 1);
            } else if o == 3 {
                t.push(i as u32 + 1);
            }
        }
        (d, t)
    }

    #[allow(clippy::too_many_arguments)]
    fn make_key(
        &self,
        n_suffix: usize,
        len_first: bool,
        len: usize,
        suffix_locs: &[u32],
        doubles: &[u32],
        triples: &[u32],
        prefixes: &[Prefix],
    ) -> Key {
        let mut mult: Vec<u32> = doubles.iter().chain(triples.iter()).copied().collect();
        mult.sort();
        let mut sl = suffix_locs.to_vec();
        sl.sort();
        let mut pl: Vec<u32> = prefixes
            .iter()
            .filter_map(|p| if let Loc::Num(n) = p.loc { Some(n) } else { None })
            .collect();
        pl.sort();
        let mut alpha: Vec<(String, u32)> = prefixes
            .iter()
            .filter_map(|p| if let Loc::Num(n) = p.loc { Some((alpha_key(&p.sub), n)) } else { None })
            .collect();
        alpha.sort();
        let (nl, nm) = (-(len as i32), -(mult.len() as i32));
        Key {
            neg_suffix: -(n_suffix as i32),
            k1: if len_first { nl } else { nm },
            k2: if len_first { nm } else { nl },
            neg_dbl: -(doubles.len() as i32),
            suffix_locs: sl,
            mult_locs: mult,
            dbl_locs: doubles.to_vec(),
            neg_nprefix: -(prefixes.len() as i32),
            prefix_locs: pl,
            alpha_locs: alpha.into_iter().map(|x| x.1).collect(),
            text: String::new(),
        }
    }

    // ------------------------------------------------ parent evaluation

    fn eval_parent(&self, atoms: &[usize], ring: Option<usize>) -> Result<Eval, String> {
        let g = self.g;
        let is_ring = ring.is_some();
        let aromatic = ring.map(|r| self.aromatic[r]).unwrap_or(false);
        let p = self.principal;
        let mut suffix_locs = vec![];
        let mut prefixes = vec![];
        let mut ester_alkyls = vec![];
        let mut carbo = false;

        let extras = |c: usize, prefixes: &mut Vec<Prefix>, ester_alkyls: &mut Vec<Sub>| -> Result<(), String> {
            match self.ctype[c] {
                CType::Group(Suffix::Ester) => {
                    let o = self.ester_link(c).unwrap();
                    let c2 = g.adj[o].iter().find(|x| x.0 != c).unwrap().0;
                    ester_alkyls.push(self.sub(o, c2)?);
                }
                CType::Group(Suffix::Amide) => {
                    let n_atom = g.adj[c].iter().find(|x| g.el[x.0] == Elem::N).unwrap().0;
                    for s in self.n_subs(n_atom, c)? {
                        prefixes.push(Prefix { loc: Loc::N, sub: s });
                    }
                }
                _ => {}
            }
            Ok(())
        };

        for (i, &a) in atoms.iter().enumerate() {
            let loc = i as u32 + 1;
            let mut consumed: Vec<usize> = vec![];
            if let CType::Group(s) = self.ctype[a] {
                if Some(s) == p && (s.is_terminal() || s == Suffix::Ketone) {
                    suffix_locs.push(loc);
                    consumed = self.group_atoms(a);
                    extras(a, &mut prefixes, &mut ester_alkyls)?;
                }
            }
            for &(b, ord) in &g.adj[a] {
                if atoms.contains(&b) || consumed.contains(&b) {
                    continue;
                }
                if p == Some(Suffix::Alcohol) && g.el[b] == Elem::O && ord == 1 && g.deg(b) == 1 {
                    suffix_locs.push(loc);
                    continue;
                }
                if p == Some(Suffix::Amine) && self.is_amine_n(b) {
                    suffix_locs.push(loc);
                    for s in self.n_subs(b, a)? {
                        prefixes.push(Prefix { loc: Loc::N, sub: s });
                    }
                    continue;
                }
                if is_ring && g.el[b] == Elem::C {
                    if let CType::Group(s) = self.ctype[b] {
                        if Some(s) == p && s.is_terminal() {
                            suffix_locs.push(loc);
                            carbo = true;
                            extras(b, &mut prefixes, &mut ester_alkyls)?;
                            continue;
                        }
                    }
                }
                prefixes.push(Prefix { loc: Loc::Num(loc), sub: self.sub(a, b)? });
            }
        }
        let (doubles, triples) = if aromatic {
            (vec![], vec![])
        } else {
            self.chain_unsat(atoms, is_ring)
        };
        let n_mult_for_key = if aromatic { 3 } else { doubles.len() + triples.len() };
        let mut key = self.make_key(
            suffix_locs.len(),
            self.style == Style::Modern || is_ring,
            atoms.len(),
            &suffix_locs,
            &doubles,
            &triples,
            &prefixes,
        );
        if self.style == Style::Modern || is_ring {
            key.k2 = -(n_mult_for_key as i32);
        } else {
            key.k1 = -(n_mult_for_key as i32);
        }
        if aromatic {
            key.neg_dbl = -3;
        }
        suffix_locs.sort();
        Ok(Eval {
            key,
            atoms: atoms.to_vec(),
            is_ring,
            aromatic,
            suffix_locs,
            carbo,
            doubles,
            triples,
            prefixes,
            ester_alkyls,
        })
    }

    fn chain_eligible_parent(&self, x: usize) -> bool {
        if self.g.el[x] != Elem::C || self.ring_id[x].is_some() {
            return false;
        }
        match self.ctype[x] {
            CType::Plain | CType::Group(Suffix::Ketone) | CType::Group(Suffix::Aldehyde) => true,
            CType::Group(s) => Some(s) == self.principal,
        }
    }

    fn assemble(&self, ev: &Eval) -> String {
        let mut body = match self.style {
            Style::Modern => self.assemble_modern(ev),
            Style::Classic => self.assemble_classic(ev),
        };
        if !ev.ester_alkyls.is_empty() {
            // "ethyl ethanoate", "diethyl butanedioate", "ethyl methyl butanedioate"
            let mut groups: Vec<(String, usize)> = vec![];
            for s in &ev.ester_alkyls {
                if let Some(g) = groups.iter_mut().find(|g| g.0 == s.text) {
                    g.1 += 1;
                } else {
                    groups.push((s.text.clone(), 1));
                }
            }
            groups.sort_by(|a, b| alpha_key(&simple(&a.0)).cmp(&alpha_key(&simple(&b.0))));
            let words: Vec<String> = groups
                .iter()
                .map(|(t, k)| {
                    if *k > 1 {
                        if t.chars().any(|c| c.is_ascii_digit()) {
                            format!("{}({t})", mult_simple(*k))
                        } else {
                            format!("{}{t}", mult_simple(*k))
                        }
                    } else {
                        t.clone()
                    }
                })
                .collect();
            body = format!("{} {}", words.join(" "), body);
        }
        body
    }

    fn assemble_classic(&self, ev: &Eval) -> String {
        let p = self.principal;
        let n = ev.atoms.len();
        let n_suffix = ev.suffix_locs.len();
        let n_mult = ev.doubles.len() + ev.triples.len();
        let num_prefix = ev.prefixes.iter().filter(|x| matches!(x.loc, Loc::Num(_))).count();
        let terminal = p.is_some_and(|s| s.is_terminal()) && !ev.carbo;
        let ring = ev.is_ring;

        let show_num = if ev.aromatic {
            num_prefix + n_suffix > 1
        } else if ring {
            num_prefix + n_suffix + n_mult > 1
        } else if n == 1 {
            false
        } else if n == 2 {
            !terminal && num_prefix + n_suffix > 1
        } else {
            true
        };
        let show_unsat = n_mult > 0
            && if ring { !(n_mult == 1 && n_suffix == 0) } else { !(n == 2 || (n == 3 && n_mult == 1)) };
        let ketone_implied =
            p == Some(Suffix::Ketone) && n_suffix == 1 && !ring && (n == 3 || (n == 4 && num_prefix == 0));
        let show_suffix = n_suffix > 0
            && !terminal
            && !ketone_implied
            && !(ring && n_suffix == 1 && n_mult == 0)
            && (ring || n > 2 || n_suffix > 1)
            && (show_num || ring);

        let prefix_str = assemble_prefixes(&ev.prefixes, show_num);
        let mut lead: Option<String> = None;
        let core: String;
        if ev.aromatic {
            match p {
                Some(s) if n_suffix == 1 && s.benzene_retained().is_some() => {
                    core = s.benzene_retained().unwrap().to_string();
                }
                Some(s) if n_suffix > 0 => {
                    let word = if ev.carbo { s.carbo_word() } else { s.chain_word() };
                    if show_num {
                        lead = Some(join_nums(&ev.suffix_locs));
                    }
                    core = attach_suffix("benzene", None, n_suffix, word);
                }
                _ => core = "benzene".to_string(),
            }
        } else {
            let st = if ring { format!("cyclo{}", stem(n)) } else { stem(n) };
            let mut hyd;
            if !ev.doubles.is_empty() && !ev.triples.is_empty() && show_unsat {
                // "1-penten-4-yne"
                hyd = st.clone();
                if ev.doubles.len() > 1 {
                    hyd.push('a');
                }
                hyd += mult_simple(ev.doubles.len());
                hyd += "en";
                hyd += &format!("-{}-{}yne", join_nums(&ev.triples), mult_simple(ev.triples.len()));
                lead = Some(join_nums(&ev.doubles));
            } else {
                hyd = hydride_name(&st, &ev.doubles, &ev.triples, false);
                if show_unsat {
                    let l = if ev.doubles.is_empty() { &ev.triples } else { &ev.doubles };
                    lead = Some(join_nums(l));
                }
            }
            core = match p {
                Some(s) if n_suffix > 0 => {
                    let word = if ev.carbo { s.carbo_word() } else { s.chain_word() };
                    if lead.is_some() {
                        let inner = if show_suffix { Some(join_nums(&ev.suffix_locs)) } else { None };
                        attach_suffix(&hyd, inner, n_suffix, word)
                    } else {
                        if show_suffix {
                            lead = Some(join_nums(&ev.suffix_locs));
                        }
                        attach_suffix(&hyd, None, n_suffix, word)
                    }
                }
                _ => hyd,
            };
        }
        match lead {
            Some(l) if prefix_str.is_empty() => format!("{l}-{core}"),
            Some(l) => format!("{prefix_str}-{l}-{core}"),
            None => format!("{prefix_str}{core}"),
        }
    }

    fn assemble_modern(&self, ev: &Eval) -> String {
        let p = self.principal;
        let n = ev.atoms.len();
        let n_suffix = ev.suffix_locs.len();
        let num_prefix = ev.prefixes.iter().filter(|x| matches!(x.loc, Loc::Num(_))).count();
        let terminal = p.is_some_and(|s| s.is_terminal()) && !ev.carbo;

        let (show_num, show_unsat) = if ev.aromatic {
            (num_prefix + n_suffix > 1, false)
        } else if ev.is_ring {
            let items = num_prefix + n_suffix + ev.doubles.len() + ev.triples.len();
            (items > 1, items > 1)
        } else if n == 1 {
            (false, false)
        } else if n == 2 {
            if terminal { (false, false) } else { (num_prefix + n_suffix > 1, false) }
        } else {
            (true, true)
        };

        let prefix_str = assemble_prefixes(&ev.prefixes, show_num);
        let suffix_locs_str = if show_num && !terminal && n_suffix > 0 {
            Some(join_nums(&ev.suffix_locs))
        } else {
            None
        };

        let parent_str = if ev.aromatic {
            match p {
                Some(s) if n_suffix == 1 && s.benzene_retained().is_some() => {
                    s.benzene_retained().unwrap().to_string()
                }
                Some(s) if n_suffix > 0 => {
                    let word = if ev.carbo { s.carbo_word() } else { s.chain_word() };
                    let locs = if show_num { Some(join_nums(&ev.suffix_locs)) } else { None };
                    attach_suffix("benzene", locs, n_suffix, word)
                }
                _ => "benzene".to_string(),
            }
        } else {
            let st = stem(n);
            let hyd = hydride_name(&st, &ev.doubles, &ev.triples, show_unsat);
            let hyd = if ev.is_ring { format!("cyclo{hyd}") } else { hyd };
            match p {
                Some(s) if n_suffix > 0 => {
                    let word = if ev.carbo { s.carbo_word() } else { s.chain_word() };
                    attach_suffix(&hyd, suffix_locs_str, n_suffix, word)
                }
                _ => hyd,
            }
        };
        format!("{prefix_str}{parent_str}")
    }

    fn best_parent(&self) -> Result<Eval, String> {
        let g = self.g;
        let consider = |ev: Eval, best: &mut Option<Eval>| {
            if best.as_ref().is_none_or(|b| ev.key < b.key) {
                *best = Some(ev);
            }
        };
        let mut first_err: Option<String> = None;
        let mut best_ring: Option<Eval> = None;
        for (ri, ring) in self.rings.iter().enumerate() {
            let n = ring.len();
            for k in 0..n {
                for dir in [1isize, -1] {
                    let atoms: Vec<usize> = (0..n)
                        .map(|i| ring[((k as isize + dir * i as isize).rem_euclid(n as isize)) as usize])
                        .collect();
                    let mut ev = match self.eval_parent(&atoms, Some(ri)) {
                        Ok(ev) => ev,
                        Err(e) => {
                            first_err.get_or_insert(e);
                            continue;
                        }
                    };
                    ev.key.text = self.assemble(&ev);
                    consider(ev, &mut best_ring);
                }
            }
        }
        // chain candidates: every directed simple path through eligible carbons
        let mut best_chain: Option<Eval> = None;
        let n = g.el.len();
        for s in 0..n {
            if !self.chain_eligible_parent(s) {
                continue;
            }
            let mut stack: Vec<Vec<usize>> = vec![vec![s]];
            while let Some(p) = stack.pop() {
                let last = *p.last().unwrap();
                for &(x, _) in &g.adj[last] {
                    if p.contains(&x) || !self.chain_eligible_parent(x) {
                        continue;
                    }
                    let mut q = p.clone();
                    q.push(x);
                    stack.push(q);
                }
                let mut ev = match self.eval_parent(&p, None) {
                    Ok(ev) => ev,
                    Err(e) => {
                        first_err.get_or_insert(e);
                        continue;
                    }
                };
                if let Some(b) = &best_chain {
                    if (ev.key.neg_suffix, ev.key.k1, ev.key.k2) > (b.key.neg_suffix, b.key.k1, b.key.k2) {
                        continue;
                    }
                }
                ev.key.text = self.assemble(&ev);
                consider(ev, &mut best_chain);
            }
        }
        match (best_ring, best_chain) {
            (Some(r), Some(c)) => {
                if r.key.neg_suffix != c.key.neg_suffix {
                    Ok(if r.key.neg_suffix < c.key.neg_suffix { r } else { c })
                } else if self.style == Style::Modern || r.atoms.len() >= c.atoms.len() {
                    Ok(r)
                } else {
                    Ok(c)
                }
            }
            (Some(r), None) => Ok(r),
            (None, Some(c)) => Ok(c),
            (None, None) => Err(first_err.unwrap_or_else(|| "Nothing to name.".to_string())),
        }
    }
}

// ---------------------------------------------------------------- public API

const COMMON: &[(&str, &str)] = &[
    ("methanal", "formaldehyde"),
    ("ethanal", "acetaldehyde"),
    ("propan-2-one", "acetone"),
    ("methanoic acid", "formic acid"),
    ("ethanoic acid", "acetic acid"),
    ("propanoic acid", "propionic acid"),
    ("butanoic acid", "butyric acid"),
    ("methylbenzene", "toluene"),
    ("ethene", "ethylene"),
    ("ethyne", "acetylene"),
    ("prop-1-ene", "propylene"),
    ("methanol", "methyl alcohol (wood alcohol)"),
    ("ethanol", "ethyl alcohol"),
    ("propan-2-ol", "isopropyl alcohol"),
    ("trichloromethane", "chloroform"),
    ("tetrachloromethane", "carbon tetrachloride"),
    ("dichloromethane", "methylene chloride"),
    ("2-methylpropane", "isobutane"),
    ("2-methylbutane", "isopentane"),
    ("2,2-dimethylpropane", "neopentane"),
    ("ethyl ethanoate", "ethyl acetate"),
    ("methyl ethanoate", "methyl acetate"),
    ("ethoxyethane", "diethyl ether"),
    ("methoxymethane", "dimethyl ether"),
    ("1,2-dimethylbenzene", "o-xylene"),
    ("1,3-dimethylbenzene", "m-xylene"),
    ("1,4-dimethylbenzene", "p-xylene"),
    ("ethenylbenzene", "styrene"),
    ("phenylmethanol", "benzyl alcohol"),
    ("ethane-1,2-diol", "ethylene glycol"),
    ("propane-1,2,3-triol", "glycerol"),
    ("ethanenitrile", "acetonitrile"),
    ("ethanamide", "acetamide"),
    ("ethanedioic acid", "oxalic acid"),
    ("2-methylbuta-1,3-diene", "isoprene"),
    ("1-phenylethan-1-one", "acetophenone"),
    ("benzene-1,4-dicarboxylic acid", "terephthalic acid"),
    ("2-hydroxybenzoic acid", "salicylic acid"),
    ("methylphenol", "cresol"),
    ("4-methylphenol", "p-cresol"),
];

pub fn common_name(name: &str) -> Option<&'static str> {
    COMMON.iter().find(|c| c.0 == name).map(|c| c.1)
}

pub fn common_to_systematic(common: &str) -> Option<&'static str> {
    let c = common.trim().to_lowercase();
    COMMON
        .iter()
        .find(|x| {
            x.1 == c
                || x.1.split(" (").next() == Some(c.as_str())
        })
        .map(|x| x.0)
}

/// Name every separate molecule in the drawing.
pub fn name_all(m: &Molecule, style: Style) -> Vec<(Vec<usize>, Result<NameResult, NameError>)> {
    m.components()
        .into_iter()
        .map(|comp| {
            let (sub, map) = m.submolecule(&comp);
            let r = name_molecule(&sub, style).map(|mut r| {
                r.parent_atoms = r.parent_atoms.iter().map(|&i| map[i]).collect();
                r
            });
            let r = r.map_err(|mut e| {
                e.atoms = e.atoms.iter().map(|&i| map[i]).collect();
                e
            });
            (comp, r)
        })
        .collect()
}

/// Name one connected molecule.
pub fn name_molecule(m: &Molecule, style: Style) -> Result<NameResult, NameError> {
    // valence
    let over: Vec<usize> = (0..m.atoms.len()).filter(|&i| m.over_valent(i)).collect();
    if !over.is_empty() {
        let a = m.atoms[over[0]];
        return Err(NameError {
            message: format!(
                "{} can make at most {} bond{}, but this one has {}. Remove a bond or lower the bond order.",
                elem_word(a.elem),
                a.elem.valence(),
                if a.elem.valence() == 1 { "" } else { "s" },
                m.bond_order_sum(over[0])
            ),
            atoms: over,
        });
    }
    // strip explicit hydrogens
    let heavy: Vec<usize> = (0..m.atoms.len()).filter(|&i| m.atoms[i].elem != Elem::H).collect();
    if heavy.is_empty() {
        return err("Only hydrogen atoms here. Draw a carbon skeleton.");
    }
    let mut idx = vec![usize::MAX; m.atoms.len()];
    for (k, &i) in heavy.iter().enumerate() {
        idx[i] = k;
    }
    let mut g = G { el: heavy.iter().map(|&i| m.atoms[i].elem).collect(), adj: vec![vec![]; heavy.len()], orig: heavy.clone() };
    for b in &m.bonds {
        let (x, y) = (idx[b.a], idx[b.b]);
        if x != usize::MAX && y != usize::MAX {
            g.adj[x].push((y, b.order));
            g.adj[y].push((x, b.order));
        }
    }
    if !g.el.contains(&Elem::C) {
        return err("There's no carbon here, so it isn't an organic compound. Add some carbon atoms.");
    }
    let namer = Namer::new(&g, style)?;
    let ev = namer.best_parent().map_err(|m| NameError { message: m, atoms: vec![] })?;
    let name = namer.assemble(&ev);
    let common = common_name(&name);

    // features & explanation
    let mut features = vec![];
    let has = |pred: &dyn Fn(usize) -> bool| (0..g.el.len()).any(pred);
    let any_double = (0..g.el.len()).any(|a| g.adj[a].iter().any(|&(b, o)| o == 2 && g.el[b] == Elem::C && g.el[a] == Elem::C && !(namer.ring_id[a].is_some() && namer.ring_id[a] == namer.ring_id[b] && namer.aromatic[namer.ring_id[a].unwrap()])));
    let any_triple = (0..g.el.len()).any(|a| g.adj[a].iter().any(|&(b, o)| o == 3 && g.el[b] == Elem::C && g.el[a] == Elem::C));
    let any_arom = namer.aromatic.iter().any(|&x| x);
    let any_halo = has(&|a| g.el[a].is_halogen());
    let any_ether = has(&|a| {
        g.el[a] == Elem::O
            && g.deg(a) == 2
            && g.adj[a].iter().all(|x| namer.ctype[x.0] != CType::Group(Suffix::Ester))
    });
    let mut groups_present: Vec<Suffix> = vec![];
    for a in 0..g.el.len() {
        if let CType::Group(s) = namer.ctype[a] {
            groups_present.push(s);
        }
        if g.el[a] == Elem::O && g.deg(a) == 1 && g.adj[a][0].1 == 1 && namer.ctype[g.adj[a][0].0] == CType::Plain {
            groups_present.push(Suffix::Alcohol);
        }
        if namer.is_amine_n(a) {
            groups_present.push(Suffix::Amine);
        }
    }
    groups_present.sort();
    groups_present.dedup();
    for s in &groups_present {
        features.push(format!("{} ({})", s.family(), s.group_desc()));
    }
    if any_double {
        features.push("Alkene (C=C double bond)".into());
    }
    if any_triple {
        features.push("Alkyne (C≡C triple bond)".into());
    }
    if any_arom {
        features.push("Benzene ring (aromatic)".into());
    }
    if any_halo {
        features.push("Halogen substituent".into());
    }
    if any_ether {
        features.push("Ether (C–O–C)".into());
    }

    let family = if let Some(p) = namer.principal {
        p.family().to_string()
    } else if any_halo {
        "Haloalkane (halogen compound)".into()
    } else if any_ether {
        "Ether".into()
    } else if any_arom {
        "Aromatic hydrocarbon (arene)".into()
    } else if any_triple {
        "Alkyne".into()
    } else if any_double {
        "Alkene".into()
    } else if !namer.rings.is_empty() {
        "Cycloalkane".into()
    } else {
        "Alkane".into()
    };

    let mut steps = vec![];
    let n = ev.atoms.len();
    if let Some(p) = namer.principal {
        steps.push(format!(
            "Main functional group: {} — it decides the ending \"{}\".",
            p.group_desc(),
            if ev.carbo { p.carbo_word() } else { p.chain_word() }
        ));
    }
    if ev.aromatic {
        steps.push("Parent: a benzene ring (six carbons, alternating double bonds).".into());
    } else if ev.is_ring {
        steps.push(format!("Parent: a ring of {n} carbons → \"cyclo{}\".", stem(n)));
    } else {
        let why = if namer.principal.is_some() { "longest chain holding the main group" } else { "longest carbon chain" };
        steps.push(format!("Parent: the {why} has {n} carbon{} → \"{}\".", if n == 1 { "" } else { "s" }, stem(n)));
    }
    if !ev.doubles.is_empty() || !ev.triples.is_empty() {
        let mut parts = vec![];
        if !ev.doubles.is_empty() {
            parts.push(format!("double bond{} at {}", if ev.doubles.len() > 1 { "s" } else { "" }, join_nums(&ev.doubles)));
        }
        if !ev.triples.is_empty() {
            parts.push(format!("triple bond{} at {}", if ev.triples.len() > 1 { "s" } else { "" }, join_nums(&ev.triples)));
        }
        steps.push(format!("Multiple bonds: {} → \"-ene\"/\"-yne\".", parts.join(", ")));
    }
    if n > 1 {
        steps.push(if namer.principal.is_some() {
            "Numbering starts from the end that gives the main group the lowest number, then multiple bonds, then side groups.".into()
        } else if !ev.doubles.is_empty() || !ev.triples.is_empty() {
            "Numbering gives multiple bonds the lowest numbers first, then side groups.".into()
        } else {
            "Numbering gives the side groups the lowest possible numbers.".into()
        });
    }
    if !ev.prefixes.is_empty() {
        let mut list: Vec<String> = ev
            .prefixes
            .iter()
            .map(|p| {
                let l = match p.loc {
                    Loc::N => "N".to_string(),
                    Loc::Num(k) => k.to_string(),
                };
                format!("{l}-{}", p.sub.text)
            })
            .collect();
        list.sort();
        steps.push(format!("Side groups (listed alphabetically, ignoring di/tri): {}.", list.join(", ")));
    }
    if !ev.ester_alkyls.is_empty() {
        steps.push("Esters are named as two words: the group on the oxygen first, then the acid part ending in \"-oate\".".into());
    }

    Ok(NameResult {
        name,
        parent_atoms: ev.atoms.iter().map(|&a| g.orig[a]).collect(),
        parent_is_ring: ev.is_ring,
        family,
        features,
        steps,
        common,
    })
}

fn elem_word(e: Elem) -> &'static str {
    match e {
        Elem::C => "Carbon",
        Elem::H => "Hydrogen",
        Elem::O => "Oxygen",
        Elem::N => "Nitrogen",
        Elem::F => "Fluorine",
        Elem::Cl => "Chlorine",
        Elem::Br => "Bromine",
        Elem::I => "Iodine",
    }
}
