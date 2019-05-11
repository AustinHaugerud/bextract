pub fn load_lines(path: &str) -> Result<Vec<String>, String> {
    use std::fs::File;
    use std::io::{BufRead, BufReader};

    let file = File::open(&path).map_err(|e| e.to_string())?;

    let buffered_lines = BufReader::new(file).lines();

    let mut lines = vec![];

    for rline in buffered_lines {
        let line = rline.map_err(|e| e.to_string())?;
        lines.push(line);
    }

    Ok(lines)
}
