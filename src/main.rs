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
    metagenome_annotate(&args.alignment_arg, args.fasta_arg);
}

fn metagenome_annotate(path: &str, fasta: &str) {

        let f = File::open(&path).expect("file not present");
        let read = BufReader::new(f);

  // struct for the genome alignment types
        struct AlignmentHold {
            id: String,
            genomefeature : String,
            start: String,
            end: String,
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
            idhold = &linehold[0];
            genomefeaturehold = &linehold[2];
            starthold = &linehold[3];
            endhold = &linehold[4].to_string();
            strandhold = &linehold[5].to_string();
            vectorhold.push(Alignmenthold{
            id:idhold,
            genomefeature: genomefeaturehold,
            start: starthold,
            end: endhold,
            strand: strandhold
            })
            }
            }

    // struct for the sequence types
        struct Sequence {
            id: String,
            sequence:String
            }

        let mut id_hold = Vec::new();
        let mut sequence_hold = Vec::new();
        let mut fasta_open = File::open(&args.fasta_arg);
        let mut fasta_read = BufReader::new(fasta_open);
        for i in fasta_read.lines() {
            if i.starts_with(">") {
            let mut line = i.expect("line not present").split(" ").collect::<Vec<&str>>();
            let mut header = line[0];
                id_hold.push(header);
            }
            if ! i.starts_with(">"){
            let line = i.expect("file not present");
            sequence_hold.push(sequence);
            }
        }

        let mut final_seq = Vec::new();
        for i in 0..id_hold.len() {
        fastaseq.push(Sequence{
        id: id_hold[i],
        sequence: sequence_hold[i],
        })
        }

    // struct for the positive strand annotations
     struct Positive {
     id: String,
     genomefeature: String,
     start:String,
     end: String,
     strand: String
     }

    // struct for the negative strand annotations
     struct Negative {
     id:String,
     genomefeature:String,
     start:String,
     end:String,
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
}
