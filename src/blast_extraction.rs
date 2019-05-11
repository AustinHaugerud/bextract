use std::sync::Arc;

pub type Lines = Arc<Vec<String>>;

pub struct BlastHits {
    hit_zones: Vec<HitsZone>,
}

impl BlastHits {
    pub fn get_hit_zones(&self) -> &Vec<HitsZone> {
        &self.hit_zones
    }
}

pub struct HitsZone {
    query: String,
    hits: Vec<BlastHit>,
}

impl HitsZone {
    pub fn get_query(&self) -> &str {
        &self.query
    }

    pub fn get_hits(&self) -> &Vec<BlastHit> {
        &self.hits
    }
}

pub struct BlastHit {
    record_ref: String,
    evalue: f64,
    subject_bounds: (usize, usize),
}

impl BlastHit {
    pub fn get_record_ref(&self) -> &str {
        &self.record_ref
    }

    pub fn get_evalue(&self) -> f64 {
        self.evalue
    }

    pub fn get_subject_bounds(&self) -> (usize, usize) {
        self.subject_bounds
    }
}

pub fn load_blast_file_hits(path: &str) -> Result<BlastHits, String> {
    use crate::file_util::load_lines;
    use rayon::prelude::*;

    let lines = Arc::new(load_lines(path)?);
    let zone_line_indices = find_zone_line_indices(&lines);

    let zones = partition_zones(zone_line_indices, lines.len());

    let hit_zones: Vec<HitsZone> = zones
        .par_iter()
        .map(|&zone| {
            let lines = lines.clone();
            process_zone(lines, zone)
        })
        .collect();

    Ok(BlastHits { hit_zones })
}

fn find_zone_line_indices(lines: &Lines) -> Vec<usize> {
    let mut indices = vec![];

    for (index, line) in lines.iter().enumerate() {
        if line.contains("Query=") {
            indices.push(index);
        }
    }

    indices
}

fn partition_zones(indices: Vec<usize>, num_lines: usize) -> Vec<(usize, usize)> {
    let mut zones = vec![];

    zones.push((indices[0], indices[1]));

    for i in 1..indices.len() - 1 {
        zones.push((indices[i], indices[i + 1] - 1));
    }

    zones.push((*indices.last().unwrap(), num_lines - 1));

    zones
}

fn process_zone(lines: Lines, zone: (usize, usize)) -> HitsZone {
    let hit_starts = find_hit_starts(&lines, zone);
    let hit_partitions = partition_hit_starts(hit_starts, lines.len());

    let mut hits_zone = HitsZone {
        query: lines[zone.0].to_string(),
        hits: vec![],
    };

    for partition in hit_partitions {
        let (record_ref, evalue, subject_bounds) = parse_hit(&lines, partition);
        hits_zone.hits.push(BlastHit {
            record_ref,
            evalue,
            subject_bounds,
        });
    }

    hits_zone
}

fn find_hit_starts(lines: &Lines, (beg, end): (usize, usize)) -> Vec<usize> {
    let mut hit_starts = vec![];

    for i in beg..=end {
        if lines.get(i).unwrap().contains("> NODE") {
            hit_starts.push(i);
        }
    }

    hit_starts
}

fn partition_hit_starts(hit_starts: Vec<usize>, num_lines: usize) -> Vec<(usize, usize)> {
    if hit_starts.is_empty() {
        return vec![];
    }

    let mut hit_partitions = vec![];
    hit_partitions.push((hit_starts[0], hit_starts[1]));

    for i in 1..hit_starts.len() - 1 {
        hit_partitions.push((hit_starts[i], hit_starts[i + 1] - 1));
    }

    hit_partitions.push((*hit_starts.last().unwrap(), num_lines));

    hit_partitions
}

fn parse_hit(lines: &Lines, partition: (usize, usize)) -> (String, f64, (usize, usize)) {
    (
        parse_hit_record_ref(&lines, partition),
        parse_hit_evalue(&lines, partition),
        parse_hit_subject_bounds(&lines, partition),
    )
}

fn parse_hit_record_ref(lines: &Lines, (beg, _): (usize, usize)) -> String {
    lines[beg].replace(">", "").trim().to_string()
}

fn parse_hit_evalue(lines: &Lines, (beg, end): (usize, usize)) -> f64 {
    let expr = retrieve_eval_expression(&lines, (beg, end));
    let parts = expr.split('=').collect::<Vec<&str>>();
    parts[1].trim().parse::<f64>().unwrap_or_else(|_| {
        panic!("Failed to parse evalue in partition: {:?}", (beg, end));
    })
}

fn parse_hit_subject_bounds(lines: &Lines, (beg, end): (usize, usize)) -> (usize, usize) {
    let subject_line = retrieve_subject(&lines, (beg, end));
    let parts = subject_line.split_ascii_whitespace().collect::<Vec<&str>>();

    let lval = parts[1].parse::<usize>().unwrap_or_else(|_| {
        panic!(
            "Failed to get subject left bound in partition: {:?}",
            (beg, end)
        );
    });

    let rval = parts[3].parse::<usize>().unwrap_or_else(|_| {
        panic!(
            "Failed to get subject right bound in partition: {:?}",
            (beg, end)
        );
    });

    (lval, rval)
}

fn retrieve_eval_expression(lines: &Lines, (beg, end): (usize, usize)) -> String {
    let mut result = None;

    for i in beg..=end {
        let line = lines.get(i).unwrap();
        if line.contains("Expect") {
            let parts = line.split(',').collect::<Vec<&str>>();
            result = Some(parts[1].to_string());
            break;
        }
    }

    result.unwrap_or_else(|| {
        panic!(
            "Failed to retrieve eval expression in partition: {:?}",
            (beg, end)
        );
    })
}

fn retrieve_subject(lines: &Lines, (beg, end): (usize, usize)) -> String {
    let mut result = None;

    for i in beg..=end {
        let line = lines.get(i).unwrap();
        if line.contains("Sbjct") {
            result = Some(line.to_string());
            break;
        }
    }

    result.unwrap_or_else(|| {
        panic!("Failed to retrieve subject in partition: {:?}", (beg, end));
    })
}
