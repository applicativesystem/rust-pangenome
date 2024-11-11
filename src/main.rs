mod args;
use args::AlignmentArgs;
use clap::Parser;
use std::fs::File;
use std::io::{Write, BufReader, BufRead};
#[allow(dead_code)]

/*
 *Author Gaurav Sablok
 *Universitat Potsdam
 *Date 2024-11-11

rust implementation of metagenome annotations. given an metagenome and
the alignment files, will write all the annotations.
 * */

fn main() {

    let args = AlignmentArgs::parse();
    metagenome_annotate(&args.alignment_arg, &args.fasta_arg);
}

fn metagenome_annotate(path: &str, fasta: &str) {

    #[derive(Debug, Clone)]
    struct AlignmentGFF {
        id: String,
        genomefeature : String,
        start: u32,
        end: u32,
        strand: String,
    }
    let mut vectorhold = Vec::new();
    let mut vectorstring = Vec::new();
    let f = File::open(&path).expect("file not present");
    let read = BufReader::new(f); 
     for gffreadline in read.lines(){
     let gffline = gffreadline
                       .expect("line not present");
            if gffline.starts_with("#") {
            continue
            } else {
                vectorhold.push(gffline)    
            }
     }
     for i in vectorhold.iter() {
        let addline = i.split("\t").collect::<Vec<&str>>();
        let mut idhold = addline[0];
        let mut genomefeaturehold = addline[2];
        let mut starthold = addline[3].to_string().parse::<u32>().expect("number not present");
        let mut endhold = addline[4].to_string().parse::<u32>().expect("number not present");
        let mut strandhold = addline[6].to_string();
        vectorstring.push(AlignmentGFF{
        id:idhold.to_string(),
        genomefeature: genomefeaturehold.to_string(),
        start: starthold,
        end: endhold,
        strand: strandhold.to_string()
     })
     }
     
     #[derive(Debug, Clone)]
     struct Positive {
     id: String,
     genomefeature: String,
     start:u32,
     end: u32,
     strand: String
     }
    
     #[derive(Debug, Clone)]
     struct Negative {
     id:String,
     genomefeature:String,
     start:u32,
     end:u32,
     strand:String,
     }

     #[derive(Debug, Clone)]
    struct Sequence {
        id: String,
        sequence:String
    }

    #[derive(Debug, Clone)]
     struct CaptureSeq {
     id: String,
     seq: String,
     strand: String,
     }

    let mut positive = Vec::new();
    let mut negative = Vec::new();
    let new_positive = vectorstring.clone();
    let new_negative = vectorstring.clone();
    for i in new_positive.into_iter() { 
        if i.strand == "+" {
            positive.push(Positive{
            id: i.id,
            genomefeature: i.genomefeature,
            start: i.start,
            end: i.end,
            strand: i.strand,
           })
           }
    }
    for i in new_negative.into_iter() {
        if i.strand == "-" {
         negative.push(Negative{
         id: i.id,
         genomefeature: i.genomefeature,
         start: i.start,
         end: i.end,
         strand: i.strand,
        })
        }
    }
    let mut header = vec![];
    let mut sequence = vec![];
    let f = File::open(&fasta).expect("file not present");
     let read = BufReader::new(f);
     for i in read.lines() {
     let line = i.expect("line not present");
     if line.starts_with(">") {
         header.push(line)
     } else {
         sequence.push(line)
     }
   }
    let mut final_seq = Vec::new();
    for i in 0..header.len() {
    final_seq.push(Sequence{
        id: header[i].to_string(),
        sequence: sequence[i].to_string(),
        })
    }

    // implementing the borrowing and the copy traits here. 


    let mut mrna_capture = Vec::new();
    for vectiter in vectorstring.iter_mut() {
        for seqiter in &mut final_seq.iter_mut() {
            if vectiter.genomefeature == "mRNA" {
                let seqhold = seqiter.sequence[vectiter.start..vectiter.end].to_string();
                mrna_capture.push( CaptureSeq{
                id: vectiter.id,
                seq: seqhold,
                strand : vectiter.strand,
      })
      }
    }
   }
   //  borrowing trait error cannot move out of `vectiter.id` which is behind a mutable reference

    let mut mrna_length: Vec<usize> = Vec::new();
    let mut cds_length: Vec<usize> = Vec::new();

    for length in cds_capture.iter() {
    let mut lengthcap = length.seq.len();
    cds_length.push(lengthcap)
    }

    for length_mrna in mrna_capture.iter() {
    let mut lengthmrna = length_mrna.seq.len();
    mrna_length.push(lengthmrna)
    }


    let mut mrna_file =  File::create("mRNA.fasta").expect("file not present");
    for i in mrna_capture.iter() {
      writeln!(mrna_file, ">{:?}\n{:?}\n", i.id, i.seq);
    }

    let mut cds_file = File::create("cds.fasta").expect("file not present");
    for i in cds_capture.iter() {
      writeln!(cds_file, ">{:?}\n{:?}\n", i.id, i.seq);
    }

    let mut mrna_positive = File::create("mRNApositive.fasta").expect("file not present");
    for i in mrna_capture.iter(){
    if i.strand == "+" {
    writeln!(mrna_positive, ">{}\n{}\n", i.id, i.seq );
    }
    }

    let mut cds_positive = File::create("cds-positive.fasta").expect("file not present");
    for i in cds_capture.iter() {
    if i.strand == "+"{
        writeln!(cds_positive, ">{}\n{}\n", i.id, i.seq);
    }
    }

    let mut cds_negative = File::create("cds-negative.fasta").expect("file not present");
    for i in cds_capture {
    if i.strand == "-" {
        writeln!(cds_negative, ">{}\n{}\n", i.id, i.seq);
    }
    }
    }
