use crate::layout;
use crate::molecule::{Elem, Molecule};
use crate::namer::{self, NameError, NameResult, Style};
use crate::parser;
use crate::render::{self, ACCENT, ACCENT_DARK, BG, BORDER, CARD, DrawOpts, ERR, INK, MUTED, accent_soft};
use egui::{Color32, CornerRadius, FontId, Pos2, Rect, RichText, Sense, Stroke, StrokeKind, Vec2, pos2, vec2};
use std::f32::consts::PI;
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

#[derive(Clone, Copy, PartialEq, Debug)]
enum Tool {
    Bond(u8),
    Chain,
    Element(Elem),
    Ring(usize, bool),
    Erase,
    Move,
}

#[derive(Clone, Copy, PartialEq)]
enum Mode {
    DrawToName,
    NameToDraw,
}

struct Drag {
    from_atom: Option<usize>,
    start: (f32, f32),
}

struct Typed {
    input: String,
    mol: Molecule,
    parent: Vec<usize>,
    parent_ring: bool,
    names: Option<Result<NameResult, NameError>>,
    other_name: Option<String>,
}

pub struct App {
    mol: Molecule,
    undo: Vec<Molecule>,
    redo: Vec<Molecule>,
    tool: Tool,
    drag: Option<Drag>,
    pan: Vec2,
    scale: f32,
    style: Style,
    show_numbers: bool,
    show_carbons: bool,
    mode: Mode,
    name_input: String,
    typed: Option<Result<Typed, String>>,
    names: Vec<(Vec<usize>, Result<NameResult, NameError>)>,
    other_names: Vec<Option<String>>,
    dirty: bool,
    toast: Option<(String, f64)>,
    frame_count: u64,
    shot_path: Option<PathBuf>,
    shot_requested: bool,
    shot_done: Arc<AtomicBool>,
}

pub struct Options {
    pub screenshot: Option<PathBuf>,
    pub name: Option<String>,
    pub canvas: Option<String>,
    pub modern: bool,
}

fn other(s: Style) -> Style {
    match s {
        Style::Classic => Style::Modern,
        Style::Modern => Style::Classic,
    }
}

fn build_from_name(name: &str, style: Style) -> Result<Typed, String> {
    let b = parser::parse_name(name)?;
    let mut mol = b.mol;
    let named = namer::name_molecule(&mol, style);
    let (parent, ring) = match &named {
        Ok(r) => (r.parent_atoms.clone(), r.parent_is_ring),
        Err(_) => (b.backbone.clone(), b.backbone_ring),
    };
    layout::layout(&mut mol, &parent, ring);
    let other_name = namer::name_molecule(&mol, other(style)).ok().map(|r| r.name);
    Ok(Typed { input: name.to_string(), mol, parent, parent_ring: ring, names: Some(named), other_name })
}

impl App {
    pub fn new(cc: &eframe::CreationContext<'_>, opts: Options) -> Self {
        setup_style(&cc.egui_ctx);
        let style = if opts.modern { Style::Modern } else { Style::Classic };
        let start = opts.canvas.clone().unwrap_or_else(|| "4-methyl-3-penten-2-ol".to_string());
        let mol = build_from_name(&start, style).map(|t| t.mol).unwrap_or_default();
        let mut app = App {
            mol,
            undo: vec![],
            redo: vec![],
            tool: Tool::Bond(1),
            drag: None,
            pan: Vec2::ZERO,
            scale: 58.0,
            style,
            show_numbers: true,
            show_carbons: false,
            mode: Mode::DrawToName,
            name_input: "2-methylbutanoic acid".into(),
            typed: None,
            names: vec![],
            other_names: vec![],
            dirty: true,
            toast: None,
            frame_count: 0,
            shot_path: opts.screenshot,
            shot_requested: false,
            shot_done: Arc::new(AtomicBool::new(false)),
        };
        if let Some(n) = opts.name {
            app.name_input = n;
            app.mode = Mode::NameToDraw;
        }
        app.typed = Some(build_from_name(&app.name_input, app.style));
        app
    }

    fn refresh(&mut self) {
        self.names = namer::name_all(&self.mol, self.style);
        self.other_names = namer::name_all(&self.mol, other(self.style))
            .into_iter()
            .map(|(_, r)| r.ok().map(|x| x.name))
            .collect();
        self.dirty = false;
    }

    fn push_undo(&mut self) {
        self.undo.push(self.mol.clone());
        if self.undo.len() > 200 {
            self.undo.remove(0);
        }
        self.redo.clear();
        self.dirty = true;
    }

    fn do_undo(&mut self) {
        if let Some(m) = self.undo.pop() {
            self.redo.push(std::mem::replace(&mut self.mol, m));
            self.dirty = true;
        }
    }

    fn do_redo(&mut self) {
        if let Some(m) = self.redo.pop() {
            self.undo.push(std::mem::replace(&mut self.mol, m));
            self.dirty = true;
        }
    }

    fn set_toast(&mut self, ctx: &egui::Context, msg: &str) {
        let t = ctx.input(|i| i.time);
        self.toast = Some((msg.to_string(), t + 3.5));
    }

    // ------------------------------------------------------------ geometry helpers

    fn hit_atom(&self, p: (f32, f32), skip: Option<usize>) -> Option<usize> {
        let mut best = None;
        let mut bd = 0.32f32;
        for (i, a) in self.mol.atoms.iter().enumerate() {
            if Some(i) == skip {
                continue;
            }
            let d = ((a.x - p.0).powi(2) + (a.y - p.1).powi(2)).sqrt();
            if d < bd {
                bd = d;
                best = Some(i);
            }
        }
        best
    }

    fn hit_bond(&self, p: (f32, f32)) -> Option<usize> {
        let mut best = None;
        let mut bd = 0.2f32;
        for (i, b) in self.mol.bonds.iter().enumerate() {
            let (a1, a2) = (self.mol.atoms[b.a], self.mol.atoms[b.b]);
            let (dx, dy) = (a2.x - a1.x, a2.y - a1.y);
            let l2 = dx * dx + dy * dy;
            if l2 < 1e-6 {
                continue;
            }
            let t = (((p.0 - a1.x) * dx + (p.1 - a1.y) * dy) / l2).clamp(0.0, 1.0);
            if !(0.12..=0.88).contains(&t) {
                continue;
            }
            let (cx, cy) = (a1.x + t * dx, a1.y + t * dy);
            let d = ((p.0 - cx).powi(2) + (p.1 - cy).powi(2)).sqrt();
            if d < bd {
                bd = d;
                best = Some(i);
            }
        }
        best
    }

    /// Best direction (radians) for a new bond from atom a.
    fn free_direction(&self, a: usize) -> f32 {
        let at = self.mol.atoms[a];
        let mut angs: Vec<f32> = self
            .mol
            .neighbors(a)
            .iter()
            .map(|(j, _)| {
                let b = self.mol.atoms[*j];
                (b.y - at.y).atan2(b.x - at.x)
            })
            .collect();
        match angs.len() {
            0 => -PI / 6.0,
            1 => {
                let back = angs[0];
                let linear = self.mol.neighbors(a).iter().any(|x| x.1 == 3);
                if linear {
                    return back + PI;
                }
                let c1 = back + 2.0 * PI / 3.0;
                let c2 = back - 2.0 * PI / 3.0;
                // zig-zag: keep the same vertical sense as the bond we came from
                if back.sin().abs() < 0.05 {
                    if c1.sin() < 0.0 { c1 } else { c2 }
                } else if c1.sin() * back.sin() > 0.0 {
                    c1
                } else {
                    c2
                }
            }
            _ => {
                angs.sort_by(|x, y| x.partial_cmp(y).unwrap());
                let mut best = (0.0, -1.0);
                for i in 0..angs.len() {
                    let s = angs[i];
                    let e = if i + 1 < angs.len() { angs[i + 1] } else { angs[0] + 2.0 * PI };
                    if e - s > best.1 {
                        best = (s + (e - s) / 2.0, e - s);
                    }
                }
                best.0
            }
        }
    }

    fn snap_end(&self, start: (f32, f32), cur: (f32, f32), from: Option<usize>) -> ((f32, f32), Option<usize>) {
        if let Some(hit) = self.hit_atom(cur, from) {
            let a = self.mol.atoms[hit];
            return ((a.x, a.y), Some(hit));
        }
        let ang = (cur.1 - start.1).atan2(cur.0 - start.0);
        let step = PI / 6.0;
        let snapped = (ang / step).round() * step;
        let end = (start.0 + snapped.cos(), start.1 + snapped.sin());
        if let Some(hit) = self.hit_atom(end, from) {
            let a = self.mol.atoms[hit];
            return ((a.x, a.y), Some(hit));
        }
        (end, None)
    }

    fn atom_at_or_new(&mut self, p: (f32, f32), existing: Option<usize>, elem: Elem) -> usize {
        match existing {
            Some(i) => i,
            None => self.mol.add_atom(elem, p.0, p.1),
        }
    }

    fn chain_points(start: (f32, f32), cur: (f32, f32)) -> Vec<(f32, f32)> {
        let (dx, dy) = (cur.0 - start.0, cur.1 - start.1);
        let dist = (dx * dx + dy * dy).sqrt();
        let n = ((dist / 0.866).round() as usize).clamp(1, 20);
        let ang = dy.atan2(dx);
        let step = PI / 6.0;
        let base = (ang / step).round() * step;
        let mut pts = vec![start];
        let mut p = start;
        for i in 0..n {
            let a = if i % 2 == 0 { base - PI / 6.0 } else { base + PI / 6.0 };
            p = (p.0 + a.cos(), p.1 + a.sin());
            pts.push(p);
        }
        pts
    }

    fn place_ring(&mut self, center: (f32, f32), size: usize, aromatic: bool, first_angle: f32) -> Vec<usize> {
        let r = 0.5 / (PI / size as f32).sin();
        let mut idx = vec![];
        for i in 0..size {
            let a = first_angle + 2.0 * PI * i as f32 / size as f32;
            let p = (center.0 + r * a.cos(), center.1 + r * a.sin());
            let existing = self.hit_atom(p, None).filter(|_| false);
            let k = self.atom_at_or_new(p, existing, Elem::C);
            idx.push(k);
        }
        for i in 0..size {
            let order = if aromatic && i % 2 == 0 { 2 } else { 1 };
            self.mol.add_bond(idx[i], idx[(i + 1) % size], order);
        }
        idx
    }

    fn erase_at(&mut self, p: (f32, f32)) -> bool {
        if let Some(a) = self.hit_atom(p, None) {
            let nbrs: Vec<usize> = self.mol.neighbors(a).iter().map(|x| x.0).collect();
            self.mol.remove_atom(a);
            // drop carbons left completely alone
            let mut lonely: Vec<usize> = nbrs
                .into_iter()
                .map(|j| if j > a { j - 1 } else { j })
                .filter(|&j| self.mol.neighbors(j).is_empty() && self.mol.atoms[j].elem == Elem::C)
                .collect();
            lonely.sort_unstable();
            for j in lonely.into_iter().rev() {
                self.mol.remove_atom(j);
            }
            return true;
        }
        if let Some(b) = self.hit_bond(p) {
            self.mol.remove_bond(b);
            return true;
        }
        false
    }

    fn cleanup_layout(&mut self) {
        if self.mol.is_empty() {
            return;
        }
        if self.names.len() != 1 {
            return;
        }
        let (parent, ring) = match &self.names[0].1 {
            Ok(r) => (r.parent_atoms.clone(), r.parent_is_ring),
            Err(_) => (vec![0], false),
        };
        self.push_undo();
        layout::layout(&mut self.mol, &parent, ring);
        self.pan = Vec2::ZERO;
    }

    // ------------------------------------------------------------ UI

    fn header(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            // logo mark: a small hexagon
            let (r, _) = ui.allocate_exact_size(vec2(34.0, 34.0), Sense::hover());
            let c = r.center();
            let pts: Vec<Pos2> = (0..6)
                .map(|i| {
                    let a = PI / 6.0 + i as f32 * PI / 3.0;
                    c + vec2(a.cos(), a.sin()) * 14.0
                })
                .collect();
            ui.painter().add(egui::Shape::convex_polygon(pts.clone(), ACCENT, Stroke::NONE));
            ui.painter().line_segment([pts[0] + (c - pts[0]) * 0.3, pts[1] + (c - pts[1]) * 0.3], Stroke::new(2.0, Color32::WHITE));
            ui.painter().line_segment([pts[2] + (c - pts[2]) * 0.3, pts[3] + (c - pts[3]) * 0.3], Stroke::new(2.0, Color32::WHITE));
            ui.painter().line_segment([pts[4] + (c - pts[4]) * 0.3, pts[5] + (c - pts[5]) * 0.3], Stroke::new(2.0, Color32::WHITE));
            ui.add_space(4.0);
            ui.vertical(|ui| {
                ui.label(RichText::new("Nomenclature Lab").size(20.0).strong().color(INK));
                ui.label(RichText::new("Draw a structure to get its name, or type a name to see the structure").size(12.5).color(MUTED));
            });
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.toggle_value(&mut self.show_carbons, "Show C and H labels");
                ui.toggle_value(&mut self.show_numbers, "Show main chain numbers");
                ui.add_space(10.0);
                let before = self.style;
                ui.selectable_value(&mut self.style, Style::Modern, "IUPAC 2013 (but-2-ene)");
                ui.selectable_value(&mut self.style, Style::Classic, "Textbook (2-butene)");
                ui.label(RichText::new("Naming style:").color(MUTED));
                if before != self.style {
                    self.dirty = true;
                    if let Some(Ok(t)) = &self.typed {
                        let input = t.input.clone();
                        self.typed = Some(build_from_name(&input, self.style));
                    }
                }
            });
        });
    }

    fn tool_button(&mut self, ui: &mut egui::Ui, tool: Tool, tip: &str) {
        let selected = self.tool == tool;
        let (rect, resp) = ui.allocate_exact_size(vec2(38.0, 38.0), Sense::click());
        let p = ui.painter();
        let fill = if selected {
            accent_soft(40)
        } else if resp.hovered() {
            Color32::from_rgb(238, 241, 245)
        } else {
            CARD
        };
        p.rect_filled(rect, 8.0, fill);
        p.rect_stroke(rect, 8.0, Stroke::new(if selected { 1.6 } else { 1.0 }, if selected { ACCENT } else { BORDER }), StrokeKind::Inside);
        let c = rect.center();
        let col = if selected { ACCENT_DARK } else { INK };
        let st = Stroke::new(2.0, col);
        match tool {
            Tool::Bond(o) => {
                let a = c + vec2(-9.0, 7.0);
                let b = c + vec2(9.0, -7.0);
                let n = (vec2(7.0, 9.0)).normalized();
                match o {
                    1 => {
                        p.line_segment([a, b], st);
                    }
                    2 => {
                        p.line_segment([a + n * 2.5, b + n * 2.5], st);
                        p.line_segment([a - n * 2.5, b - n * 2.5], st);
                    }
                    _ => {
                        p.line_segment([a, b], st);
                        p.line_segment([a + n * 4.0, b + n * 4.0], st);
                        p.line_segment([a - n * 4.0, b - n * 4.0], st);
                    }
                }
            }
            Tool::Chain => {
                let pts = vec![c + vec2(-12.0, 4.0), c + vec2(-6.0, -4.0), c + vec2(0.0, 4.0), c + vec2(6.0, -4.0), c + vec2(12.0, 4.0)];
                p.line(pts, st);
            }
            Tool::Element(e) => {
                p.text(c, egui::Align2::CENTER_CENTER, e.symbol(), FontId::proportional(17.0), render::elem_color(e));
            }
            Tool::Ring(n, arom) => {
                let r = 11.0;
                let pts: Vec<Pos2> = (0..n)
                    .map(|i| {
                        let a = -PI / 2.0 + 2.0 * PI * i as f32 / n as f32;
                        c + vec2(a.cos(), a.sin()) * r
                    })
                    .collect();
                p.add(egui::Shape::closed_line(pts.clone(), st));
                if arom {
                    for i in [0, 2, 4] {
                        let a = pts[i] + (c - pts[i]) * 0.3;
                        let b = pts[(i + 1) % n] + (c - pts[(i + 1) % n]) * 0.3;
                        p.line_segment([a, b], Stroke::new(1.6, col));
                    }
                }
            }
            Tool::Erase => {
                let r = Rect::from_center_size(c, vec2(16.0, 10.0));
                let pts = [r.left_bottom() + vec2(3.0, 0.0), r.left_top() + vec2(7.0, 0.0), r.right_top(), r.right_bottom() - vec2(4.0, 0.0)];
                p.add(egui::Shape::closed_line(pts.to_vec(), st));
                p.line_segment([pts[1] + vec2(5.0, 0.0), pts[0] + vec2(5.0, 0.0)], st);
            }
            Tool::Move => {
                for d in [vec2(1.0, 0.0), vec2(-1.0, 0.0), vec2(0.0, 1.0), vec2(0.0, -1.0)] {
                    p.arrow(c, d * 12.0, Stroke::new(1.8, col));
                }
            }
        }
        if resp.clicked() {
            self.tool = tool;
        }
        resp.on_hover_text(tip);
    }

    fn toolbar(&mut self, ui: &mut egui::Ui) {
        ui.horizontal_wrapped(|ui| {
            ui.spacing_mut().item_spacing = vec2(5.0, 6.0);
            self.tool_button(ui, Tool::Bond(1), "Single bond — drag to draw a line, click a bond to cycle single/double/triple (1)");
            self.tool_button(ui, Tool::Bond(2), "Double bond (2)");
            self.tool_button(ui, Tool::Bond(3), "Triple bond (3)");
            self.tool_button(ui, Tool::Chain, "Chain — drag to draw a zig-zag carbon chain");
            ui.add_space(6.0);
            for e in [Elem::C, Elem::O, Elem::N, Elem::Cl, Elem::Br, Elem::F, Elem::I, Elem::H] {
                self.tool_button(ui, Tool::Element(e), &format!("{} — click an atom to change it, or drag from an atom to add one", e.symbol()));
            }
            ui.add_space(6.0);
            for n in [3, 4, 5, 6] {
                self.tool_button(ui, Tool::Ring(n, false), &format!("{n}-carbon ring — click empty space, or click an atom to attach it"));
            }
            self.tool_button(ui, Tool::Ring(6, true), "Benzene ring");
            ui.add_space(6.0);
            self.tool_button(ui, Tool::Erase, "Eraser — click or drag over atoms and bonds");
            self.tool_button(ui, Tool::Move, "Move atoms — drag an atom (right-drag anywhere to pan, scroll to zoom)");
        });
        ui.add_space(4.0);
        ui.horizontal(|ui| {
            if ui.add_enabled(!self.undo.is_empty(), egui::Button::new("↶ Undo")).clicked() {
                self.do_undo();
            }
            if ui.add_enabled(!self.redo.is_empty(), egui::Button::new("↷ Redo")).clicked() {
                self.do_redo();
            }
            if ui.button("Tidy drawing").on_hover_text("Re-draw neatly with the main chain as a zig-zag").clicked() {
                self.cleanup_layout();
            }
            if ui.button("Clear").clicked() && !self.mol.is_empty() {
                self.push_undo();
                self.mol = Molecule::default();
                self.pan = Vec2::ZERO;
            }
        });
    }

    fn canvas(&mut self, ui: &mut egui::Ui) {
        let size = ui.available_size() - vec2(0.0, 26.0);
        let (rect, resp) = ui.allocate_exact_size(size.max(vec2(200.0, 200.0)), Sense::click_and_drag());
        let painter = ui.painter_at(rect);
        painter.rect_filled(rect, 10.0, CARD);
        painter.rect_stroke(rect, 10.0, Stroke::new(1.0, BORDER), StrokeKind::Inside);
        // dot grid
        let sp = self.scale * 0.5;
        let origin = rect.center() + self.pan;
        let gx0 = ((rect.left() - origin.x) / sp).ceil() as i32;
        let gx1 = ((rect.right() - origin.x) / sp).floor() as i32;
        let gy0 = ((rect.top() - origin.y) / sp).ceil() as i32;
        let gy1 = ((rect.bottom() - origin.y) / sp).floor() as i32;
        if (gx1 - gx0) * (gy1 - gy0) < 20000 {
            for gx in gx0..=gx1 {
                for gy in gy0..=gy1 {
                    painter.circle_filled(pos2(origin.x + gx as f32 * sp, origin.y + gy as f32 * sp), 1.0, Color32::from_rgb(214, 219, 226));
                }
            }
        }

        let scale = self.scale;
        let to_screen = move |x: f32, y: f32| pos2(origin.x + x * scale, origin.y + y * scale);
        let to_mol = |p: Pos2| ((p.x - origin.x) / scale, (p.y - origin.y) / scale);

        // zoom & pan
        if resp.hovered() {
            let (scroll, zoom) = ui.input(|i| (i.smooth_scroll_delta.y, i.zoom_delta()));
            let f = zoom * (1.0 + scroll * 0.0015);
            if (f - 1.0).abs() > 1e-4 {
                self.scale = (self.scale * f).clamp(18.0, 140.0);
            }
        }
        if resp.dragged_by(egui::PointerButton::Secondary) || resp.dragged_by(egui::PointerButton::Middle) {
            self.pan += resp.drag_delta();
        }

        let pointer = resp.interact_pointer_pos().or(resp.hover_pos());
        let mp = pointer.map(to_mol);
        let mut hover_atom = None;
        let mut hover_bond = None;
        if let Some(p) = mp {
            if resp.hovered() || resp.dragged() {
                hover_atom = self.hit_atom(p, None);
                if hover_atom.is_none() {
                    hover_bond = self.hit_bond(p);
                }
            }
        }

        // primary drag
        if resp.drag_started_by(egui::PointerButton::Primary) {
            let origin_p = ui.input(|i| i.pointer.press_origin()).unwrap_or(pointer.unwrap_or(rect.center()));
            let start_m = to_mol(origin_p);
            let from = self.hit_atom(start_m, None);
            let start = from.map(|a| (self.mol.atoms[a].x, self.mol.atoms[a].y)).unwrap_or(start_m);
            match self.tool {
                Tool::Move => {
                    if from.is_some() {
                        self.push_undo();
                    }
                }
                Tool::Erase => {
                    self.push_undo();
                    self.erase_at(start_m);
                }
                _ => {}
            }
            self.drag = Some(Drag { from_atom: from, start });
        }
        if resp.dragged_by(egui::PointerButton::Primary) {
            if let (Some(d), Some(p)) = (&self.drag, mp) {
                match self.tool {
                    Tool::Move => {
                        if let Some(a) = d.from_atom {
                            if a < self.mol.atoms.len() {
                                self.mol.atoms[a].x = p.0;
                                self.mol.atoms[a].y = p.1;
                                self.dirty = true;
                            }
                        }
                    }
                    Tool::Erase => {
                        if self.erase_at(p) {
                            self.dirty = true;
                        }
                    }
                    _ => {}
                }
            }
        }
        if resp.drag_stopped_by(egui::PointerButton::Primary) {
            if let (Some(d), Some(p)) = (self.drag.take(), mp) {
                match self.tool {
                    Tool::Bond(_) | Tool::Element(_) => {
                        let (o, elem) = match self.tool {
                            Tool::Bond(o) => (o, Elem::C),
                            Tool::Element(e) => (1, e),
                            _ => (1, Elem::C),
                        };
                        let dist = ((p.0 - d.start.0).powi(2) + (p.1 - d.start.1).powi(2)).sqrt();
                        if dist > 0.3 {
                            let (end, end_atom) = self.snap_end(d.start, p, d.from_atom);
                            if end_atom.is_none() || end_atom != d.from_atom {
                                self.push_undo();
                                let a = self.atom_at_or_new(d.start, d.from_atom, Elem::C);
                                let b = self.atom_at_or_new(end, end_atom, elem);
                                if a != b {
                                    self.mol.add_bond(a, b, o);
                                }
                            }
                        }
                    }
                    Tool::Chain => {
                        let pts = Self::chain_points(d.start, p);
                        if pts.len() > 1 {
                            self.push_undo();
                            let mut prev = self.atom_at_or_new(d.start, d.from_atom, Elem::C);
                            for q in &pts[1..] {
                                let ex = self.hit_atom(*q, Some(prev));
                                let k = self.atom_at_or_new(*q, ex, Elem::C);
                                if k != prev {
                                    self.mol.add_bond(prev, k, 1);
                                }
                                prev = k;
                            }
                        }
                    }
                    _ => {}
                }
            }
            self.drag = None;
        }
        if resp.clicked() {
            if let Some(p) = mp {
                self.click(ui.ctx(), p);
            }
        }

        // keyboard shortcuts (when no text box has focus)
        let typing = ui.ctx().egui_wants_keyboard_input();
        if !typing {
            let (k1, k2, k3, undo, redo) = ui.input(|i| {
                (
                    i.key_pressed(egui::Key::Num1),
                    i.key_pressed(egui::Key::Num2),
                    i.key_pressed(egui::Key::Num3),
                    i.modifiers.command && !i.modifiers.shift && i.key_pressed(egui::Key::Z),
                    i.modifiers.command && i.modifiers.shift && i.key_pressed(egui::Key::Z),
                )
            });
            if k1 {
                self.tool = Tool::Bond(1);
            }
            if k2 {
                self.tool = Tool::Bond(2);
            }
            if k3 {
                self.tool = Tool::Bond(3);
            }
            if undo {
                self.do_undo();
            }
            if redo {
                self.do_redo();
            }
        }

        if self.dirty {
            self.refresh();
        }

        // what to highlight
        let mut parent: Vec<usize> = vec![];
        let mut parent_ring = false;
        let mut errs: Vec<usize> = vec![];
        for (_, r) in &self.names {
            match r {
                Ok(n) => {
                    if self.names.len() == 1 {
                        parent = n.parent_atoms.clone();
                        parent_ring = n.parent_is_ring;
                    }
                }
                Err(e) => errs.extend(e.atoms.iter().copied()),
            }
        }
        let opts = DrawOpts {
            show_carbons: self.show_carbons,
            parent: if self.show_numbers { &parent } else { &[] },
            parent_ring,
            show_numbers: self.show_numbers,
            error_atoms: &errs,
            hover_atom: if self.drag.is_none() { hover_atom } else { None },
            hover_bond: if self.drag.is_none() { hover_bond } else { None },
        };
        render::draw_molecule(&painter, &self.mol, &to_screen, scale, &opts);

        // drag preview
        if let (Some(d), Some(p)) = (&self.drag, mp) {
            let preview = Stroke::new(2.0, accent_soft(200));
            match self.tool {
                Tool::Bond(_) | Tool::Element(_) => {
                    let (end, _) = self.snap_end(d.start, p, d.from_atom);
                    painter.line_segment([to_screen(d.start.0, d.start.1), to_screen(end.0, end.1)], preview);
                    painter.circle_filled(to_screen(end.0, end.1), 4.0, ACCENT);
                }
                Tool::Chain => {
                    let pts: Vec<Pos2> = Self::chain_points(d.start, p).iter().map(|q| to_screen(q.0, q.1)).collect();
                    let n = pts.len() - 1;
                    let last = *pts.last().unwrap();
                    painter.line(pts, preview);
                    painter.text(last + vec2(10.0, -10.0), egui::Align2::LEFT_BOTTOM, format!("{n} bond{}", if n == 1 { "" } else { "s" }), FontId::proportional(13.0), ACCENT_DARK);
                }
                _ => {}
            }
        }

        if self.mol.is_empty() {
            painter.text(
                rect.center(),
                egui::Align2::CENTER_CENTER,
                "Drag here to draw a bond — each line is a carbon–carbon bond",
                FontId::proportional(16.0),
                MUTED,
            );
        }

        // toast
        let now = ui.input(|i| i.time);
        if let Some((msg, until)) = &self.toast {
            if now < *until {
                let g = painter.layout_no_wrap(msg.clone(), FontId::proportional(14.0), Color32::WHITE);
                let r = Rect::from_center_size(pos2(rect.center().x, rect.bottom() - 30.0), g.size() + vec2(24.0, 14.0));
                painter.rect_filled(r, 8.0, Color32::from_rgba_unmultiplied(31, 41, 51, 230));
                painter.galley(r.min + vec2(12.0, 7.0), g, Color32::WHITE);
                ui.ctx().request_repaint();
            } else {
                self.toast = None;
            }
        }

        let hint = match self.tool {
            Tool::Bond(_) => "Drag to draw a bond · click an atom to grow the chain · click a bond to make it single → double → triple",
            Tool::Chain => "Drag to draw a zig-zag chain; the number of bonds is shown while you drag",
            Tool::Element(_) => "Click an atom to change its element · drag from an atom to add this element · click empty space for a lone atom",
            Tool::Ring(..) => "Click empty space to place a ring · click an atom to attach a ring to it",
            Tool::Erase => "Click or drag over atoms and bonds to remove them",
            Tool::Move => "Drag atoms to move them · right-drag to pan · scroll to zoom",
        };
        ui.add_space(4.0);
        ui.label(RichText::new(hint).size(12.5).color(MUTED));
    }

    fn click(&mut self, ctx: &egui::Context, p: (f32, f32)) {
        let atom = self.hit_atom(p, None);
        let bond = if atom.is_none() { self.hit_bond(p) } else { None };
        match self.tool {
            Tool::Bond(o) => {
                if let Some(a) = atom {
                    self.push_undo();
                    let ang = self.free_direction(a);
                    let at = self.mol.atoms[a];
                    let q = (at.x + ang.cos(), at.y + ang.sin());
                    let ex = self.hit_atom(q, Some(a));
                    let b = self.atom_at_or_new(q, ex, Elem::C);
                    self.mol.add_bond(a, b, o);
                } else if let Some(bi) = bond {
                    self.push_undo();
                    let cur = self.mol.bonds[bi].order;
                    self.mol.bonds[bi].order = if o == 1 { cur % 3 + 1 } else if cur == o { 1 } else { o };
                } else {
                    self.push_undo();
                    let a = self.mol.add_atom(Elem::C, p.0, p.1);
                    let ang = -PI / 6.0;
                    let b = self.mol.add_atom(Elem::C, p.0 + ang.cos(), p.1 + ang.sin());
                    self.mol.add_bond(a, b, o);
                }
            }
            Tool::Chain => {
                if atom.is_none() && bond.is_none() {
                    self.set_toast(ctx, "Drag to draw a chain");
                }
            }
            Tool::Element(e) => {
                if let Some(a) = atom {
                    self.push_undo();
                    let cur = self.mol.atoms[a].elem;
                    self.mol.atoms[a].elem = if cur == e && e != Elem::C { Elem::C } else { e };
                } else if bond.is_none() {
                    self.push_undo();
                    self.mol.add_atom(e, p.0, p.1);
                }
            }
            Tool::Ring(n, arom) => {
                if bond.is_some() {
                    self.set_toast(ctx, "Rings sharing a bond (fused rings) aren't supported — click an atom instead");
                    return;
                }
                self.push_undo();
                let r = 0.5 / (PI / n as f32).sin();
                if let Some(a) = atom {
                    let ang = self.free_direction(a);
                    let at = self.mol.atoms[a];
                    let first = (at.x + ang.cos(), at.y + ang.sin());
                    let center = (first.0 + ang.cos() * r, first.1 + ang.sin() * r);
                    let idx = self.place_ring(center, n, arom, ang + PI);
                    self.mol.add_bond(a, idx[0], 1);
                } else {
                    self.place_ring(p, n, arom, -PI / 2.0);
                }
            }
            Tool::Erase => {
                self.push_undo();
                if !self.erase_at(p) {
                    self.undo.pop();
                }
            }
            Tool::Move => {}
        }
    }

    fn card<R>(ui: &mut egui::Ui, add: impl FnOnce(&mut egui::Ui) -> R) -> R {
        egui::Frame::new()
            .fill(CARD)
            .stroke(Stroke::new(1.0, BORDER))
            .corner_radius(10.0)
            .inner_margin(14.0)
            .show(ui, |ui| {
                ui.set_width(ui.available_width());
                add(ui)
            })
            .inner
    }

    fn mode_switch(&mut self, ui: &mut egui::Ui) {
        let w = ui.available_width();
        let (rect, _) = ui.allocate_exact_size(vec2(w, 42.0), Sense::hover());
        ui.painter().rect_filled(rect, 10.0, Color32::from_rgb(229, 233, 239));
        let halves = [
            (Rect::from_min_max(rect.min, pos2(rect.center().x, rect.max.y)), Mode::DrawToName, "Drawing  →  Name"),
            (Rect::from_min_max(pos2(rect.center().x, rect.min.y), rect.max), Mode::NameToDraw, "Name  →  Drawing"),
        ];
        for (r, m, label) in halves {
            let r = r.shrink(4.0);
            let resp = ui.interact(r, ui.scope_id().with(label), Sense::click());
            let sel = self.mode == m;
            if sel {
                ui.painter().rect_filled(r, 8.0, CARD);
                ui.painter().rect_stroke(r, 8.0, Stroke::new(1.0, BORDER), StrokeKind::Inside);
            }
            ui.painter().text(
                r.center(),
                egui::Align2::CENTER_CENTER,
                label,
                FontId::proportional(15.0),
                if sel { ACCENT_DARK } else if resp.hovered() { INK } else { MUTED },
            );
            if resp.clicked() {
                self.mode = m;
            }
        }
    }

    fn info_rows(ui: &mut egui::Ui, mol: &Molecule, family: Option<&str>) {
        egui::Grid::new(ui.scope_id().with("info")).num_columns(2).spacing(vec2(16.0, 6.0)).show(ui, |ui| {
            ui.label(RichText::new("Formula").color(MUTED));
            let job = render::formula_job(&mol.formula(), 16.0, INK);
            ui.label(job);
            ui.end_row();
            ui.label(RichText::new("Molar mass").color(MUTED));
            ui.label(RichText::new(format!("{:.2} g/mol", mol.molar_mass())).color(INK));
            ui.end_row();
            if let Some(f) = family {
                ui.label(RichText::new("Family").color(MUTED));
                ui.label(RichText::new(f).color(INK));
                ui.end_row();
            }
        });
    }

    fn name_card(&self, ui: &mut egui::Ui, idx: usize, total: usize, r: &NameResult, other_name: Option<&String>) {
        Self::card(ui, |ui| {
            let title = if total > 1 { format!("MOLECULE {} — NAME", idx + 1) } else { "NAME".to_string() };
            ui.label(RichText::new(title).size(11.5).color(MUTED).strong());
            ui.add_space(2.0);
            ui.horizontal_wrapped(|ui| {
                ui.add(egui::Label::new(RichText::new(&r.name).size(27.0).strong().color(INK)).selectable(true).wrap());
            });
            ui.horizontal(|ui| {
                if ui.small_button("Copy name").clicked() {
                    ui.ctx().copy_text(r.name.clone());
                }
            });
            if let Some(o) = other_name {
                if o != &r.name {
                    let label = match self.style {
                        Style::Classic => "IUPAC 2013 style:",
                        Style::Modern => "Textbook style:",
                    };
                    ui.horizontal_wrapped(|ui| {
                        ui.label(RichText::new(label).color(MUTED));
                        ui.label(RichText::new(o).color(INK));
                    });
                }
            }
            if let Some(c) = r.common {
                ui.horizontal_wrapped(|ui| {
                    ui.label(RichText::new("Common name:").color(MUTED));
                    ui.label(RichText::new(c).color(INK));
                });
            }
        });
    }

    fn explain_card(ui: &mut egui::Ui, r: &NameResult) {
        Self::card(ui, |ui| {
            ui.label(RichText::new("HOW IT'S NAMED").size(11.5).color(MUTED).strong());
            ui.add_space(4.0);
            for (i, s) in r.steps.iter().enumerate() {
                ui.horizontal_wrapped(|ui| {
                    ui.label(RichText::new(format!("{}.", i + 1)).color(ACCENT_DARK).strong());
                    ui.label(RichText::new(s).color(INK));
                });
            }
            if !r.features.is_empty() {
                ui.add_space(6.0);
                ui.label(RichText::new("GROUPS FOUND").size(11.5).color(MUTED).strong());
                ui.horizontal_wrapped(|ui| {
                    for f in &r.features {
                        egui::Frame::new()
                            .fill(accent_soft(28))
                            .corner_radius(CornerRadius::same(12))
                            .inner_margin(egui::Margin::symmetric(9, 3))
                            .show(ui, |ui| {
                                ui.label(RichText::new(f).size(12.5).color(ACCENT_DARK));
                            });
                    }
                });
            }
        });
    }

    fn error_card(ui: &mut egui::Ui, msg: &str) {
        egui::Frame::new()
            .fill(Color32::from_rgb(254, 242, 242))
            .stroke(Stroke::new(1.0, Color32::from_rgb(252, 202, 202)))
            .corner_radius(10.0)
            .inner_margin(14.0)
            .show(ui, |ui| {
                ui.set_width(ui.available_width());
                ui.label(RichText::new("CAN'T NAME THIS YET").size(11.5).color(ERR).strong());
                ui.add_space(2.0);
                ui.label(RichText::new(msg).color(INK));
            });
    }

    fn draw_to_name(&mut self, ui: &mut egui::Ui) {
        if self.dirty {
            self.refresh();
        }
        if self.mol.is_empty() {
            Self::card(ui, |ui| {
                ui.label(RichText::new("START DRAWING").size(11.5).color(MUTED).strong());
                ui.add_space(4.0);
                for s in [
                    "Drag on the canvas to draw a line. Every corner and line end is a carbon atom.",
                    "Click a line to turn it into a double, then a triple bond.",
                    "Pick O, N or a halogen and click an atom to change it.",
                    "The name appears here as you draw.",
                ] {
                    ui.label(RichText::new(format!("•  {s}")).color(INK));
                }
            });
            return;
        }
        let names = self.names.clone();
        let total = names.len();
        for (i, (comp, r)) in names.iter().enumerate() {
            match r {
                Ok(res) => {
                    self.name_card(ui, i, total, res, self.other_names.get(i).and_then(|x| x.as_ref()));
                    ui.add_space(8.0);
                    if total == 1 {
                        Self::card(ui, |ui| Self::info_rows(ui, &self.mol, Some(&res.family)));
                        ui.add_space(8.0);
                        Self::explain_card(ui, res);
                    }
                }
                Err(e) => {
                    let msg = if total > 1 { format!("Molecule {}: {}", i + 1, e.message) } else { e.message.clone() };
                    Self::error_card(ui, &msg);
                    let _ = comp;
                }
            }
            ui.add_space(8.0);
        }
        if total > 1 {
            ui.label(RichText::new("Tip: the drawing has separate pieces, so each one is named on its own.").size(12.5).color(MUTED));
        }
    }

    fn name_to_draw(&mut self, ui: &mut egui::Ui) {
        Self::card(ui, |ui| {
            ui.label(RichText::new("TYPE A NAME").size(11.5).color(MUTED).strong());
            ui.add_space(4.0);
            let mut go = false;
            ui.horizontal(|ui| {
                let w = ui.available_width() - 74.0;
                let resp = ui.add(
                    egui::TextEdit::singleline(&mut self.name_input)
                        .hint_text("e.g. 3-methyl-1-butene")
                        .font(FontId::proportional(17.0))
                        .desired_width(w),
                );
                if resp.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)) {
                    go = true;
                }
                if ui.add_sized(vec2(64.0, 30.0), egui::Button::new(RichText::new("Draw").strong().color(Color32::WHITE)).fill(ACCENT)).clicked() {
                    go = true;
                }
            });
            ui.add_space(6.0);
            ui.label(RichText::new("Try one:").size(12.5).color(MUTED));
            ui.horizontal_wrapped(|ui| {
                for ex in [
                    "2-methylbutane",
                    "2-butene",
                    "2-propanol",
                    "propanone",
                    "ethyl ethanoate",
                    "cyclohexanol",
                    "4-chlorophenol",
                    "3-ethyl-2,4-dimethylhexane",
                    "1-penten-4-yne",
                    "N-methylethanamine",
                    "acetone",
                ] {
                    if ui.small_button(ex).clicked() {
                        self.name_input = ex.to_string();
                        go = true;
                    }
                }
            });
            if go {
                self.typed = Some(build_from_name(&self.name_input.clone(), self.style));
            }
        });
        ui.add_space(8.0);

        let mut send_to_canvas = false;
        match &self.typed {
            None => {}
            Some(Err(e)) => Self::error_card(ui, e),
            Some(Ok(t)) => {
                Self::card(ui, |ui| {
                    ui.label(RichText::new("STRUCTURE").size(11.5).color(MUTED).strong());
                    let w = ui.available_width();
                    let (rect, _) = ui.allocate_exact_size(vec2(w, 280.0), Sense::hover());
                    let painter = ui.painter_at(rect);
                    painter.rect_filled(rect, 8.0, Color32::from_rgb(250, 251, 252));
                    if let Some((x0, y0, x1, y1)) = t.mol.bounds() {
                        let bw = (x1 - x0).max(0.5) + 1.6;
                        let bh = (y1 - y0).max(0.5) + 1.6;
                        let sc = (rect.width() / bw).min(rect.height() / bh).min(56.0);
                        let c = rect.center();
                        let (mx, my) = ((x0 + x1) / 2.0, (y0 + y1) / 2.0);
                        let ts = move |x: f32, y: f32| pos2(c.x + (x - mx) * sc, c.y + (y - my) * sc);
                        render::draw_molecule(
                            &painter,
                            &t.mol,
                            &ts,
                            sc,
                            &DrawOpts {
                                show_carbons: self.show_carbons,
                                parent: if self.show_numbers { &t.parent } else { &[] },
                                parent_ring: t.parent_ring,
                                show_numbers: self.show_numbers,
                                error_atoms: &[],
                                hover_atom: None,
                                hover_bond: None,
                            },
                        );
                    }
                    ui.add_space(6.0);
                    if ui.button("Open this on the canvas to edit  →").clicked() {
                        send_to_canvas = true;
                    }
                });
                ui.add_space(8.0);
                if let Some(Ok(r)) = &t.names {
                    let typed_norm: String = t.input.to_lowercase().chars().filter(|c| !c.is_whitespace()).collect();
                    let matches = |s: &str| s.to_lowercase().chars().filter(|c| !c.is_whitespace()).collect::<String>() == typed_norm;
                    let ok = matches(&r.name) || t.other_name.as_deref().is_some_and(matches);
                    let is_common = !ok && namer::common_to_systematic(&t.input).is_some();
                    Self::card(ui, |ui| {
                        ui.label(RichText::new("STANDARD NAME").size(11.5).color(MUTED).strong());
                        ui.add(egui::Label::new(RichText::new(&r.name).size(24.0).strong().color(INK)).selectable(true).wrap());
                        if ok {
                            ui.label(RichText::new("✔ The name you typed is the standard name.").color(ACCENT_DARK));
                        } else if is_common {
                            ui.label(RichText::new(format!("\"{}\" is a common name; above is its systematic name.", t.input.trim())).color(MUTED));
                        } else {
                            ui.label(
                                RichText::new(format!(
                                    "\"{}\" describes this structure, but the standard name is the one above (it uses the correct main chain and lowest numbers).",
                                    t.input.trim()
                                ))
                                .color(Color32::from_rgb(180, 83, 9)),
                            );
                        }
                        if let Some(o) = &t.other_name {
                            if o != &r.name {
                                let label = match self.style {
                                    Style::Classic => "IUPAC 2013 style:",
                                    Style::Modern => "Textbook style:",
                                };
                                ui.horizontal_wrapped(|ui| {
                                    ui.label(RichText::new(label).color(MUTED));
                                    ui.label(RichText::new(o).color(INK));
                                });
                            }
                        }
                        if let Some(c) = r.common {
                            ui.horizontal_wrapped(|ui| {
                                ui.label(RichText::new("Common name:").color(MUTED));
                                ui.label(RichText::new(c).color(INK));
                            });
                        }
                        ui.add_space(6.0);
                        Self::info_rows(ui, &t.mol, Some(&r.family));
                    });
                    ui.add_space(8.0);
                    Self::explain_card(ui, r);
                } else if let Some(Err(e)) = &t.names {
                    Self::error_card(ui, &e.message);
                }
            }
        }
        if send_to_canvas {
            if let Some(Ok(t)) = &self.typed {
                let m = t.mol.clone();
                self.push_undo();
                self.mol = m;
                self.pan = Vec2::ZERO;
                self.mode = Mode::DrawToName;
            }
        }
    }

    fn output_panel(&mut self, ui: &mut egui::Ui) {
        self.mode_switch(ui);
        ui.add_space(10.0);
        egui::ScrollArea::vertical().auto_shrink([false, false]).show(ui, |ui| match self.mode {
            Mode::DrawToName => self.draw_to_name(ui),
            Mode::NameToDraw => self.name_to_draw(ui),
        });
    }

    fn handle_screenshot(&mut self, ctx: &egui::Context) {
        let Some(path) = self.shot_path.clone() else { return };
        self.frame_count += 1;
        ctx.request_repaint();
        if self.frame_count >= 120 && !self.shot_requested {
            self.shot_requested = true;
            let done = self.shot_done.clone();
            ctx.request_screenshot(move |img| {
                let [w, h] = img.size;
                let mut buf = Vec::with_capacity(w * h * 4);
                for p in &img.pixels {
                    buf.extend_from_slice(&[p.r(), p.g(), p.b(), 255]);
                }
                if let Some(parent) = path.parent() {
                    let _ = std::fs::create_dir_all(parent);
                }
                match image::RgbaImage::from_raw(w as u32, h as u32, buf) {
                    Some(im) => {
                        if let Err(e) = im.save(&path) {
                            eprintln!("could not save screenshot: {e}");
                        }
                    }
                    None => eprintln!("bad screenshot buffer"),
                }
                done.store(true, Ordering::SeqCst);
            });
        }
        if self.shot_done.load(Ordering::SeqCst) || self.frame_count > 900 {
            ctx.send_viewport_cmd(egui::ViewportCommand::Close);
        }
    }
}

fn setup_style(ctx: &egui::Context) {
    let mut visuals = egui::Visuals::light();
    visuals.panel_fill = BG;
    visuals.window_fill = CARD;
    visuals.extreme_bg_color = Color32::WHITE;
    visuals.selection.bg_fill = accent_soft(60);
    visuals.selection.stroke = Stroke::new(1.0, ACCENT_DARK);
    visuals.widgets.inactive.corner_radius = CornerRadius::same(6);
    visuals.widgets.hovered.corner_radius = CornerRadius::same(6);
    visuals.widgets.active.corner_radius = CornerRadius::same(6);
    ctx.set_visuals(visuals);
    ctx.global_style_mut(|s| {
        use egui::TextStyle::*;
        s.text_styles.insert(Body, FontId::proportional(14.5));
        s.text_styles.insert(Button, FontId::proportional(14.0));
        s.text_styles.insert(Small, FontId::proportional(12.0));
        s.text_styles.insert(Heading, FontId::proportional(22.0));
        s.spacing.button_padding = vec2(9.0, 5.0);
        s.spacing.item_spacing = vec2(8.0, 6.0);
    });
}

impl eframe::App for App {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        let ctx = ui.ctx().clone();
        self.handle_screenshot(&ctx);
        if self.dirty {
            self.refresh();
        }
        egui::Panel::top("header")
            .frame(egui::Frame::new().fill(CARD).inner_margin(egui::Margin::symmetric(18, 10)).stroke(Stroke::new(1.0, BORDER)))
            .show(ui, |ui| self.header(ui));
        egui::Panel::right("output")
            .resizable(true)
            .default_size(470.0)
            .min_size(360.0)
            .frame(egui::Frame::new().fill(BG).inner_margin(egui::Margin { left: 8, right: 16, top: 14, bottom: 14 }))
            .show(ui, |ui| self.output_panel(ui));
        egui::CentralPanel::default()
            .frame(egui::Frame::new().fill(BG).inner_margin(egui::Margin { left: 16, right: 8, top: 14, bottom: 10 }))
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.label(RichText::new("CANVAS").size(11.5).color(MUTED).strong());
                    ui.label(RichText::new("— draw your structure here").size(12.5).color(MUTED));
                });
                ui.add_space(4.0);
                self.toolbar(ui);
                ui.add_space(6.0);
                self.canvas(ui);
            });
    }
}
