#![cfg(test)]
use crate::molecule::{Elem, Molecule};
use crate::namer::{molecule_rings, name_molecule, Style};
use crate::parser::parse_name;

/// Minimal SMILES reader for tests (C c O N F Cl Br I, = #, branches, ring digits).
fn smiles(s: &str) -> Molecule {
    let mut m = Molecule::default();
    let ch: Vec<char> = s.chars().collect();
    let mut i = 0;
    let mut prev: Option<usize> = None;
    let mut stack: Vec<Option<usize>> = vec![];
    let mut order = 1u8;
    let mut rings: std::collections::HashMap<char, (usize, u8)> = Default::default();
    let mut arom = vec![];
    while i < ch.len() {
        let c = ch[i];
        let mut atom: Option<(Elem, bool)> = None;
        match c {
            '(' => stack.push(prev),
            ')' => prev = stack.pop().unwrap(),
            '=' => order = 2,
            '#' => order = 3,
            '0'..='9' => {
                let p = prev.unwrap();
                if let Some((q, o)) = rings.remove(&c) {
                    m.add_bond(p, q, o.max(order));
                } else {
                    rings.insert(c, (p, order));
                }
                order = 1;
            }
            'C' if ch.get(i + 1) == Some(&'l') => {
                atom = Some((Elem::Cl, false));
                i += 1;
            }
            'B' if ch.get(i + 1) == Some(&'r') => {
                atom = Some((Elem::Br, false));
                i += 1;
            }
            'C' => atom = Some((Elem::C, false)),
            'c' => atom = Some((Elem::C, true)),
            'O' => atom = Some((Elem::O, false)),
            'N' => atom = Some((Elem::N, false)),
            'F' => atom = Some((Elem::F, false)),
            'I' => atom = Some((Elem::I, false)),
            _ => panic!("bad smiles {s}"),
        }
        if let Some((e, a)) = atom {
            let k = m.add_atom(e, 0.0, 0.0);
            arom.push(a);
            if let Some(p) = prev {
                m.add_bond(p, k, order);
            }
            order = 1;
            prev = Some(k);
        }
        i += 1;
    }
    for r in molecule_rings(&m) {
        if r.len() == 6 && r.iter().all(|&a| arom[a]) {
            for k in [0, 2, 4] {
                let bi = m.bond_between(r[k], r[k + 1]).unwrap();
                m.bonds[bi].order = 2;
            }
        }
    }
    m
}

fn nm(s: &str, st: Style) -> String {
    match name_molecule(&smiles(s), st) {
        Ok(r) => r.name,
        Err(e) => format!("ERR: {}", e.message),
    }
}

const CASES: &[(&str, &str, &str)] = &[
    ("C", "methane", "methane"),
    ("CCCC", "butane", "butane"),
    ("CC(C)C", "2-methylpropane", "2-methylpropane"),
    ("CC(C)(C)C", "2,2-dimethylpropane", "2,2-dimethylpropane"),
    ("CCC(C)CC", "3-methylpentane", "3-methylpentane"),
    ("CCC(CC)C(C)C", "3-ethyl-2-methylpentane", "3-ethyl-2-methylpentane"),
    ("CCC(C(C)C)CCC", "3-ethyl-2-methylhexane", "3-ethyl-2-methylhexane"),
    ("CCCC(C)C", "2-methylpentane", "2-methylpentane"),
    ("C=CC", "propene", "prop-1-ene"),
    ("CC=CC", "2-butene", "but-2-ene"),
    ("C=CC=C", "1,3-butadiene", "buta-1,3-diene"),
    ("C#CC", "propyne", "prop-1-yne"),
    ("C=CC#C", "1-buten-3-yne", "but-1-en-3-yne"),
    ("C=C(CC)CCC", "2-ethyl-1-pentene", "3-methylidenehexane"),
    ("CC(C)=C", "2-methylpropene", "2-methylprop-1-ene"),
    ("CCO", "ethanol", "ethanol"),
    ("CC(O)C", "2-propanol", "propan-2-ol"),
    ("OCCO", "1,2-ethanediol", "ethane-1,2-diol"),
    ("CC(C)CO", "2-methyl-1-propanol", "2-methylpropan-1-ol"),
    ("CC(C)(C)O", "2-methyl-2-propanol", "2-methylpropan-2-ol"),
    ("CC=CC(O)C", "3-penten-2-ol", "pent-3-en-2-ol"),
    ("NCCO", "2-aminoethanol", "2-aminoethan-1-ol"),
    ("CC(=O)C", "propanone", "propan-2-one"),
    ("CCC(=O)C", "butanone", "butan-2-one"),
    ("CC(=O)CC(=O)C", "2,4-pentanedione", "pentane-2,4-dione"),
    ("CCCC=O", "butanal", "butanal"),
    ("C=CC=O", "propenal", "prop-2-enal"),
    ("C=O", "methanal", "methanal"),
    ("CC(=O)O", "ethanoic acid", "ethanoic acid"),
    ("CC(C)C(=O)O", "2-methylpropanoic acid", "2-methylpropanoic acid"),
    ("OC(=O)CC(=O)O", "propanedioic acid", "propanedioic acid"),
    ("ClCC(=O)O", "chloroethanoic acid", "chloroethanoic acid"),
    ("CCOC(=O)C", "ethyl ethanoate", "ethyl ethanoate"),
    ("COC(=O)CCC", "methyl butanoate", "methyl butanoate"),
    ("CCN", "ethanamine", "ethanamine"),
    ("CC(N)C", "2-propanamine", "propan-2-amine"),
    ("CCNC", "N-methylethanamine", "N-methylethanamine"),
    ("CC#N", "ethanenitrile", "ethanenitrile"),
    ("CC(N)=O", "ethanamide", "ethanamide"),
    ("CCCl", "chloroethane", "chloroethane"),
    ("CC(Cl)C", "2-chloropropane", "2-chloropropane"),
    ("ClC(Cl)Cl", "trichloromethane", "trichloromethane"),
    ("BrCC(C)CCl", "1-bromo-3-chloro-2-methylpropane", "1-bromo-3-chloro-2-methylpropane"),
    ("COC", "methoxymethane", "methoxymethane"),
    ("CCOCC", "ethoxyethane", "ethoxyethane"),
    ("CCCOC", "1-methoxypropane", "1-methoxypropane"),
    ("C1CCCCC1", "cyclohexane", "cyclohexane"),
    ("CC1CCCCC1", "methylcyclohexane", "methylcyclohexane"),
    ("OC1CCCCC1", "cyclohexanol", "cyclohexanol"),
    ("CC1CCCCC1O", "2-methylcyclohexanol", "2-methylcyclohexan-1-ol"),
    ("C1=CCCCC1", "cyclohexene", "cyclohexene"),
    ("CC1=CCCCC1", "1-methylcyclohexene", "1-methylcyclohex-1-ene"),
    ("CC1C=CCCC1", "3-methylcyclohexene", "3-methylcyclohex-1-ene"),
    ("OC1CCC(=O)CC1", "4-hydroxycyclohexanone", "4-hydroxycyclohexan-1-one"),
    ("OC(=O)C1CCCCC1", "cyclohexanecarboxylic acid", "cyclohexanecarboxylic acid"),
    ("CC(C)C1CCCCC1", "isopropylcyclohexane", "(propan-2-yl)cyclohexane"),
    ("c1ccccc1", "benzene", "benzene"),
    ("Cc1ccccc1", "methylbenzene", "methylbenzene"),
    ("Oc1ccccc1", "phenol", "phenol"),
    ("Oc1ccc(Cl)cc1", "4-chlorophenol", "4-chlorophenol"),
    ("OC(=O)c1ccccc1", "benzoic acid", "benzoic acid"),
    ("Cc1ccccc1C", "1,2-dimethylbenzene", "1,2-dimethylbenzene"),
    ("Nc1ccccc1", "aniline", "aniline"),
    ("OCc1ccccc1", "phenylmethanol", "phenylmethanol"),
    ("CCOC(=O)c1ccccc1", "ethyl benzoate", "ethyl benzoate"),
    ("C=Cc1ccccc1", "ethenylbenzene", "ethenylbenzene"),
    ("CCCCCCCc1ccccc1", "1-phenylheptane", "heptylbenzene"),
    ("OCCc1ccccc1", "2-phenylethanol", "2-phenylethan-1-ol"),
    ("CC(C)(C)CC(C)C", "2,2,4-trimethylpentane", "2,2,4-trimethylpentane"),
    ("CCC(C)(C)CC", "3,3-dimethylpentane", "3,3-dimethylpentane"),
    ("CC(Br)C(C)Cl", "2-bromo-3-chlorobutane", "2-bromo-3-chlorobutane"),
];

#[test]
fn names_from_structures() {
    let mut bad = vec![];
    for (smi, classic, modern) in CASES {
        let c = nm(smi, Style::Classic);
        let m = nm(smi, Style::Modern);
        if &c != classic {
            bad.push(format!("{smi}: classic got {c:?}, want {classic:?}"));
        }
        if &m != modern {
            bad.push(format!("{smi}: modern got {m:?}, want {modern:?}"));
        }
    }
    assert!(bad.is_empty(), "\n{}", bad.join("\n"));
}

#[test]
fn round_trip_names() {
    let mut bad = vec![];
    for (smi, classic, modern) in CASES {
        for (name, st) in [(classic, Style::Classic), (modern, Style::Modern)] {
            match parse_name(name) {
                Ok(b) => {
                    let got = name_molecule(&b.mol, st).map(|r| r.name).unwrap_or_else(|e| format!("ERR {}", e.message));
                    if &got != name {
                        bad.push(format!("{name} ({smi}) -> {got}"));
                    }
                }
                Err(e) => bad.push(format!("{name}: parse error {e}")),
            }
        }
    }
    assert!(bad.is_empty(), "\n{}", bad.join("\n"));
}

#[test]
fn alternative_spellings() {
    let pairs = [
        ("2-ethylbutane", "3-methylpentane"),
        ("4-methylpentane", "2-methylpentane"),
        ("acetone", "propanone"),
        ("toluene", "methylbenzene"),
        ("isopropyl alcohol", "2-propanol"),
        ("methylpropane", "2-methylpropane"),
        ("dimethylpropane", "2,2-dimethylpropane"),
        ("2-methyl-2-butene", "2-methyl-2-butene"),
        ("2-methylbut-2-ene", "2-methyl-2-butene"),
        ("3-methyl-1-butanol", "3-methyl-1-butanol"),
        ("2-(1-methylethyl)pentane", "2,3-dimethylhexane"),
        ("ethyl acetate", "ethyl ethanoate"),
        ("1,3-butadiene", "1,3-butadiene"),
        ("2,4-hexadien-1-ol", "2,4-hexadien-1-ol"),
        ("cyclohex-2-en-1-ol", "2-cyclohexen-1-ol"),
        ("tert-butylbenzene", "tert-butylbenzene"),
        ("3-(prop-2-en-1-yl)hexane", "4-ethyl-1-heptene"),
        ("4-(dimethylamino)butan-2-one", "4-(dimethylamino)-2-butanone"),
        ("N,N-dimethylethanamine", "N,N-dimethylethanamine"),
        ("methyl 2-methylpropanoate", "methyl 2-methylpropanoate"),
        ("1,4-benzenediol", "1,4-benzenediol"),
        ("benzene-1,4-dicarboxylic acid", "1,4-benzenedicarboxylic acid"),
        ("2-methoxypropane", "2-methoxypropane"),
        ("3-ethyl-2,4-dimethylhexane", "3-ethyl-2,4-dimethylhexane"),
        ("pent-1-en-4-yne", "1-penten-4-yne"),
    ];
    let mut bad = vec![];
    for (inp, want) in pairs {
        match parse_name(inp) {
            Ok(b) => {
                let got = name_molecule(&b.mol, Style::Classic).map(|r| r.name).unwrap_or_else(|e| format!("ERR {}", e.message));
                if got != want {
                    bad.push(format!("{inp} -> {got} (want {want})"));
                }
            }
            Err(e) => bad.push(format!("{inp}: {e}")),
        }
    }
    assert!(bad.is_empty(), "\n{}", bad.join("\n"));
}
