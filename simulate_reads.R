library(Biostrings)
library(ballgown)
library(polyester)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==3) {
    args[4] = "cdnaf"
}

fasta_file = args[3]
gtfpath = args[1]
fasta = readDNAStringSet(fasta_file)
gtf = gffRead(gtfpath)


simulate_experiment(gtf=gtf, seqpath=fasta, num_reps=1, fold_changes=1, outdir=args[2], readlen=100,
                    reads_per_transcript=1000, distr="normal", error_model="illumina5", bias=args[4])
