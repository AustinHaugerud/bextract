extern crate bio;
extern crate clap;
extern crate rayon;

mod blast_extraction;
mod fasta;
mod file_util;

use clap::{App, Arg, ArgMatches};

struct AppArguments {
    blast_input_path: String,
    sequence_input_path: String,
    emax: f64,
    output_path: String,
    extension: usize,
}

impl AppArguments {
    fn parse_arguments(matches: &ArgMatches) -> Result<AppArguments, Vec<String>> {
        let mut errors = vec![];

        let blast_input_path = matches.value_of("blast_input").unwrap().to_string();
        let sequence_input_path = matches.value_of("sequence_input").unwrap().to_string();
        let output_path = matches.value_of("output_path").unwrap().to_string();
        let emax_str = matches.value_of("emax").unwrap();
        let extension_str = matches.value_of("extension").unwrap();

        let valid_blast_input = std::path::Path::new(&blast_input_path).is_file();

        if !valid_blast_input {
            errors.push(String::from("Could not obtain blast input file."));
        }

        let valid_sequence_input = std::path::Path::new(&sequence_input_path).is_file();

        if !valid_sequence_input {
            errors.push(String::from("Could not obtain sequence input file."));
        }

        let emax_result = emax_str.parse::<f64>();

        if emax_result.is_err() {
            errors.push(String::from("Invalid emax argument."));
        }

        let extension_result = extension_str.parse::<usize>();

        if extension_result.is_err() {
            errors.push(String::from("Invalid extension argument."));
        }

        if valid_blast_input
            && valid_sequence_input
            && emax_result.is_ok()
            && extension_result.is_ok()
        {
            Ok(AppArguments {
                blast_input_path,
                sequence_input_path,
                emax: emax_result.unwrap(),
                output_path,
                extension: extension_result.unwrap(),
            })
        } else {
            Err(errors)
        }
    }
}

fn main() {
    let app = App::new("bextract")
        .version("1.0")
        .author("Austin Jenkins")
        .about("Extract sequences referencing blast hits and criteria.")
        .arg(
            Arg::with_name("blast_input")
                .short("b")
                .long("binput")
                .help("The blast input file.")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("sequence_input")
                .short("c")
                .long("sinput")
                .help("The fasta file to extract from.")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("emax")
                .short("e")
                .long("emax")
                .help("The maximum evalue allowed before discarding.")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("output_path")
                .short("o")
                .long("output")
                .help("The output file path.")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("extension")
                .short("p")
                .long("extension")
                .help("How much to extend when extracting sequence based on subject bounds.")
                .takes_value(true)
                .required(true),
        )
        .get_matches();

    match AppArguments::parse_arguments(&app) {
        Ok(arguments) => process(&arguments),
        Err(errors) => display_errors(&errors),
    }
}

fn process(args: &AppArguments) {
    use crate::fasta::{load_fasta, Fasta};
    use bio::io::fasta;
    use std::fs::File;
    use std::io::BufWriter;
    use std::io::Write;

    let fasta_data = load_fasta(&args.sequence_input_path).expect("Failed to get fasta data.");

    let output = File::create(&args.output_path)
        .expect(&format!("Failed to create file at: {}", &args.output_path));
    let writer = BufWriter::new(output);
    let mut fasta_writer = fasta::Writer::new(writer);

    use blast_extraction::load_blast_file_hits;
    let blast_hits =
        load_blast_file_hits(&args.blast_input_path).expect("Failed to load blast hits.");

    let mut wanted_hits = vec![];

    for zone in blast_hits.get_hit_zones() {
        for hit in zone.get_hits() {
            if hit.get_evalue() <= args.emax {
                wanted_hits.push(hit);
            }
        }
    }

    for hit in wanted_hits {
        if let Some(sequence_data) = fasta_data.extract_sequence(
            hit.get_record_ref(),
            hit.get_subject_bounds(),
            args.extension,
        ) {
            let text_slice = bio::utils::TextSlice::from(&sequence_data);
            let mut record = fasta::Record::with_attrs(hit.get_record_ref(), None, text_slice);
            fasta_writer
                .write_record(&record)
                .expect("Failed to write fasta record.");
        }
    }
}

fn display_errors(errors: &[String]) {
    for error in errors.iter() {
        println!("Error: {}", error);
    }
}
