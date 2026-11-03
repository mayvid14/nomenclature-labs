mod app;
mod layout;
mod molecule;
mod namer;
mod parser;
mod render;
mod tests;

use std::path::PathBuf;

fn main() -> eframe::Result {
    let args: Vec<String> = std::env::args().collect();
    let mut opts = app::Options { screenshot: None, name: None, canvas: None, modern: false };
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--screenshot" if i + 1 < args.len() => {
                opts.screenshot = Some(PathBuf::from(&args[i + 1]));
                i += 1;
            }
            "--name" if i + 1 < args.len() => {
                opts.name = Some(args[i + 1].clone());
                i += 1;
            }
            "--canvas" if i + 1 < args.len() => {
                opts.canvas = Some(args[i + 1].clone());
                i += 1;
            }
            "--modern" => opts.modern = true,
            _ => {}
        }
        i += 1;
    }
    let native = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1360.0, 860.0])
            .with_min_inner_size([900.0, 600.0])
            .with_title("Nomenclature Lab"),
        ..Default::default()
    };
    eframe::run_native("Nomenclature Lab", native, Box::new(move |cc| Ok(Box::new(app::App::new(cc, opts)))))
}
