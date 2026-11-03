mod app;
mod layout;
mod molecule;
mod namer;
mod parser;
mod render;
mod tests;

#[cfg(not(target_arch = "wasm32"))]
use std::path::PathBuf;

#[cfg(not(target_arch = "wasm32"))]
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

#[cfg(target_arch = "wasm32")]
fn main() {
    use eframe::wasm_bindgen::JsCast as _;

    let params = web_sys::window()
        .and_then(|w| w.location().search().ok())
        .map(|q| web_sys::UrlSearchParams::new_with_str(&q).unwrap())
        .unwrap();
    let opts = app::Options {
        screenshot: None,
        name: params.get("name"),
        canvas: params.get("canvas"),
        modern: params.has("modern"),
    };
    wasm_bindgen_futures::spawn_local(async move {
        let canvas = web_sys::window()
            .and_then(|w| w.document())
            .and_then(|d| d.get_element_by_id("the_canvas_id"))
            .and_then(|e| e.dyn_into::<web_sys::HtmlCanvasElement>().ok())
            .expect("canvas #the_canvas_id not found");
        let result = eframe::WebRunner::new()
            .start(canvas, eframe::WebOptions::default(), Box::new(move |cc| Ok(Box::new(app::App::new(cc, opts)))))
            .await;
        if let Some(loading) = web_sys::window().and_then(|w| w.document()).and_then(|d| d.get_element_by_id("loading")) {
            match result {
                Ok(_) => loading.remove(),
                Err(e) => loading.set_inner_html(&format!("The app failed to start: {e:?}")),
            }
        }
    });
}
