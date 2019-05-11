use std::collections::HashMap;

pub struct Fasta {
    sequences: HashMap<String, Vec<u8>>,
}

impl Fasta {
    pub fn extract_sequence(
        &self,
        record_id: &str,
        (beg, end): (usize, usize),
        extension: usize,
    ) -> Option<Vec<u8>> {
        let seq = self.sequences.get(record_id)?;

        let nbeg = if beg < extension { 0 } else { beg - extension };

        let nend = {
            let mut index = end + extension;
            if index > seq.len() - 1 {
                index = seq.len() - 1;
            }
            index
        };

        let new_slice = &seq[nbeg..nend];
        Some(new_slice.to_vec())
    }
}

pub fn load_fasta(path: &str) -> Result<Fasta, String> {
    use bio::io::fasta;

    let reader = fasta::Reader::from_file(path).map_err(|e| e.to_string())?;

    let mut sequences = HashMap::new();

    for result in reader.records() {
        let record = result.map_err(|e| e.to_string())?;
        let id = record.id().to_string();
        let seq = record.seq().to_vec();
        sequences.insert(id, seq);
    }

    Ok(Fasta { sequences })
}
