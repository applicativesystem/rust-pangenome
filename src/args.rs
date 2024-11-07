use clap::Parser;

#[derive(Debug, Parser)]
#[clap(version)]

pub struct AlignmentArgs {
    /// please provide the kmer to be searched for the origin
    pub alignment_arg: String,
    /// please provide the path to be searched for the strings containing the kmer
    pub fasta_arg: String,
}
