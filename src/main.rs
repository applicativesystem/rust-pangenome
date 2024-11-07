mod args;
use args::AlignmentArgs;
use clap::Parser;
use std::fs::File;
use std::io::{Write, BufReader, BufRead};
use std::collections::HashSet;
#[allow(dead_code)]

/*
 *Author Gaurav Sablok
 *Universitat Potsdam
 *Date 2024-11-7

rust implementation of metagenome annotations. given an metagenome and
the alignment files, will write all the annotations.
 * */

fn main() {

    let args = AlignmentArgs::parse();
    metagenome_annotate(&args.alignment_arg, &args.fasta_arg);
}

fn metagenome_annotate(path: &str, fasta: &str) {

        let f = File::open(&path).expect("file not present");
        let read = BufReader::new(f);

  // struct for the genome alignment types
        #[derive(Debug)]
        struct AlignmentGFF {
            id: String,
            genomefeature : String,
            start: u32,
            end: u32,
            strand: String,
        }

  // puting the struct alignment to the use
        let mut vectorhold = Vec::new();
        for i in read.lines() {
            let line = i
                       .expect("line not present");
            if line.starts_with("#") {
            continue
            } else {
            let mut linehold = line.to_string().split("\t").collect::<Vec<&str>>();
            let mut idhold = &linehold[0];
            let mut genomefeaturehold = &linehold[2];
            let mut starthold = &linehold[3].to_string().parse::<u32>();
            let mut endhold = &linehold[4].to_string().parse::<u32>();
            let mut strandhold = &linehold[5].to_string();
            vectorhold.push(AlignmentGFF{
            id:idhold.to_string(),
            genomefeature: genomefeaturehold.to_string(),
            start: starthold.parse::<u32>().unwrap(),
            end: endhold.parse::<u32>().unwrap(),
            strand: strandhold.to_string()
            })
            }
            }

    // struct for the sequence types
        #[derive(Debug)]
        struct Sequence {
            id: String,
            sequence:String
            }

        let mut id_hold = Vec::new();
        let mut sequence_hold = Vec::new();
        let mut fasta_open = File::open(&args.fasta_arg).expect("file not present");
        let mut fasta_read = BufReader::new(fasta_open);
        for i in fasta_read.lines() {
            if i.starts_with(">") {
            let mut line = i.expect("line not present").split(" ").collect::<Vec<&str>>();
            let mut header = line[0];
                id_hold.push(header);
            }
            if ! i.starts_with(">"){
            let line = i.expect("file not present");
            sequence_hold.push(line);
            }
        }

        let mut final_seq = Vec::new();
        for i in 0..id_hold.len() {
        final_seq.push(Sequence{
        id: id_hold[i].to_string(),
        sequence: sequence_hold[i].to_string(),
        })
        }

    // struct for the positive strand annotations
     #[derive(Debug)]
     struct Positive {
     id: String,
     genomefeature: String,
     start:u32,
     end: u32,
     strand: String
     }

    // struct for the negative strand annotations
     #[derive(Debug)]
     struct Negative {
     id:String,
     genomefeature:String,
     start:u32,
     end:u32,
     strand:String,
     }

     // putting the struct to the use

     let mut positive = Vec::new();
     let mut negative = Vec::new();

     for add in vectorhold.iter() {
       if i.strand == "+" {
        positive.push(Positive{
        id: add.id,
        genomefeature: add.genomefeature,
        start: add.start,
        end: add.end,
        strand: add.strand,
       })
       }
       if i.strand == "-" {
        negative.push(Negative{
        id: i.id,
        genomefeature: i.genomefeature,
        start: i.start,
        end: i.end,
        strand: i.strand})
       }
       }

    // struct for the sequence classification
     struct CaptureSeq {
     id: String,
     seq: String}

     let mut mrna_capture = Vec::new();
     let mut cds_capture = Vec::new();

     for i in sequence_hold.iter() {
      for j in vectorhold.iter() {
      if j.hold == "mRNA" {
        mrna_capture.push( CaptureSeq{
        id: i.id,
        seq: j.seq[i.start:i.end]}),
      }
      }
      }
      if j.hold == "CDS" {
      cds_capture.push(CaptureSeq{
        id: i.id,
        seq: j.seq[i.start: i.end],
      })
      }

      let mut mrna_file = File::create("mRNA.fasta").expect("file not present");
      for i in mrna_capture.iter() {
      writeln!(mrna_file, ">{:?}\n{:?}\n", i.id, i.seq)
      }

      let cds_file = File::create("cds.fasta").expect("file not present");
      for i in cds_capture.iter() {
      writeln!(cds_file, ">{:?}\n{:?}\n", i.id, i.seq)
      }
}
