library(Biostrings)
library(ballgown)
library(polyester)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==3) {
    args[3] = "cdnaf"
}

fasta_file = "/pine/scr/s/h/shengjie/hg19.fa"
gtfpath = args[1]
fasta = readDNAStringSet(fasta_file)
gtf = gffRead(gtfpath)


simulate_experiment(gtf=gtf, seqpath=fasta, num_reps=1, fold_changes=1, outdir=args[2], readlen=100,
                    reads_per_transcript=2000, distr="normal", error_model="illumina5",
                    error_rate=0.0026, bias=args[3])
