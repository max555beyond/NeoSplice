import numpy as np
import os
import argparse
import collections
import itertools
from BCBio import GFF
import logging

class Gene(object):
    def __init__(self, chromosome, startposition, endposition, strand, ID):
        self.chromosome = chromosome
        self.startposition= startposition
        self.endposition = endposition
        self.strand = strand
        self.frame = '.'
        self.transcripts = []
        self.score = '.'
        self.ID = ID
        self.left_splices = set()
        self.right_splices = set()
        self.all_splices = set()


class Transcript(object):
    def __init__(self, chromosome, startposition, endposition, strand, ID, genename, cds_list, start_codon, stop_codon,
                 splice_list, cds_set):
        self.chromosome = chromosome
        self.startposition = startposition
        self.endposition = endposition
        self.strand = strand
        self.ID = ID
        self.genename = genename
        self.cds_list = cds_list
        self.start_codon = start_codon
        self.stop_codon = stop_codon
        self.splice_list = splice_list
        self.score = '.'
        self.frame = '.'
        self.cds_set = cds_set


strand_symbol = {+1: '+', -1: '-'}


def _find_annotated_transcripts(gff_in_file, chromosome):
    in_file = gff_in_file
    in_handle = open(in_file)
    limit_info = dict(gff_id=[chromosome])
    transcripts_dic = {}
    for chromosome_inf in GFF.parse(in_handle, limit_info=limit_info):
        for gene in chromosome_inf.features:
            if "protein_coding" in gene.qualifiers["gene_type"]:
                gene_name = gene.qualifiers["gene_name"][0]
                gene_data = Gene(chromosome, gene.location.start.position, gene.location.end.position, gene.strand,
                                 gene_name)
                transcripts_dic[gene_name] = gene_data
                qualified = False
                for transcript in gene.sub_features:
                    if "protein_coding" in transcript.qualifiers["transcript_type"]:
                        start_codon = []
                        stop_codon = []
                        cds_list = []
                        splice_list = []
                        strand = transcript.strand

                        for feature in transcript.sub_features:
                            if feature.type == "start_codon":
                                start_codon.append(feature)
                            elif feature.type == "stop_codon":
                                stop_codon.append(feature)
                            elif feature.type == "CDS":
                                cds_list.append(feature)

                        cds_list.sort(key=lambda x: x.location.start.position)
                        start_codon.sort(key=lambda x: x.location.start.position)
                        stop_codon.sort(key=lambda x: x.location.start.position)
                        cds_set = set([(cds.location.start.position, cds.location.end.position) for cds in cds_list])
                        i = 0

                        while i < len(cds_list) - 1:
                            splice_list.append(
                                (cds_list[i].location.end.position, cds_list[i + 1].location.start.position))
                            i += 1

                        if len(start_codon) != 0 and len(cds_list) > 3:
                            gene_data.transcripts.append(Transcript(chromosome,
                                                                    transcript.location.start.position,
                                                                    transcript.location.end.position, strand,
                                                                    transcript.qualifiers["transcript_id"][0],
                                                                    gene.qualifiers["gene_name"][0], cds_list,
                                                                    start_codon, stop_codon, splice_list, cds_set))

                            qualified = True
                            gene_data.all_splices |= set(splice_list)
                            print gene_data.all_splices
                            gene_data.left_splices |= set([splice[0] for splice in splice_list])
                            print gene_data.left_splices
                            gene_data.right_splices |= set([splice[1] for splice in splice_list])
                            print gene_data.right_splices

                if not qualified:
                    transcripts_dic.pop(gene_name)

    return transcripts_dic


def write_gtf(gene_data, tumor_outf, normal_outf):
    transcript_left_splices = gene_data.left_splices
    transcript_right_splices = gene_data.right_splices

    print transcript_left_splices
    print transcript_right_splices
    possible_novel_splices = [(a, b) for a, b in
                              itertools.product(list(transcript_left_splices), list(transcript_right_splices)) if
                              a < b and (a, b) not in gene_data.all_splices]
    print possible_novel_splices

    tumor_outf.write('\t'.join(
        [gene_data.chromosome, "test", "gene", str(gene_data.startposition + 1), str(gene_data.endposition),
         gene_data.score, strand_symbol[gene_data.strand],
         gene_data.frame,
         """gene_id "{}";""".format(gene_data.ID) + ";" + """ transcript_id "{}";""".format('')]) + '\n')

    normal_outf.write('\t'.join(
        [gene_data.chromosome, "test", "gene", str(gene_data.startposition + 1), str(gene_data.endposition),
         gene_data.score, strand_symbol[gene_data.strand],
         gene_data.frame,
         """gene_id "{}";""".format(gene_data.ID) + ";" + """ transcript_id "{}";""".format('')]) + '\n')
    # in one-based coordinate
    novel_splice_idx = np.random.choice(len(possible_novel_splices))
    novel_splice = possible_novel_splices[novel_splice_idx]
    logging.info("selected novel splice {}".format(novel_splice))

    for transcript in gene_data.transcripts:
        if novel_splice[0] <= transcript.cds_list[0].location.start.position or novel_splice[1] >= transcript.cds_list[
                -1].location.end.position or novel_splice in transcript.splice_list:
            continue

        left_index = None
        right_index = None
        left_type = None
        right_type = None

        print "transcipt name " + transcript.ID
        print "transcipt splices " + str(transcript.splice_list)
        print "exon list " + str(transcript.cds_list)
        print "left splices " + str(gene_data.left_splices)
        print "right splices " + str(gene_data.right_splices)


        tumor_outf.write('\t'.join(
            [transcript.chromosome, "test", "transcript", str(transcript.startposition + 1),
             str(transcript.endposition),
             transcript.score, strand_symbol[transcript.strand], str(transcript.frame),
             """gene_id "{}";""".format(gene_data.ID) + ";" + """ transcript_id "{}";""".format(
                 transcript.ID)]) + '\n')

        normal_outf.write('\t'.join(
            [transcript.chromosome, "test", "transcript", str(transcript.startposition + 1),
             str(transcript.endposition),
             transcript.score, strand_symbol[transcript.strand], str(transcript.frame),
             """gene_id "{}";""".format(gene_data.ID) + ";" + """ transcript_id "{}";""".format(
                 transcript.ID)]) + '\n')

        for i, exon in enumerate(transcript.cds_list):
            #not preceded by an intron
            if exon.location.start.position < novel_splice[0] <= exon.location.end.position - 1:
                left_index = i
                left_type = "exon"
            # not superceded by an intron
            if exon.location.start.position <= novel_splice[1] - 1 < exon.location.end.position - 1:
                right_index = i
                right_type = "exon"

        for i, intron in enumerate(transcript.splice_list):
            if intron[0] <= novel_splice[0] <= intron[1] - 1:
                left_index = i
                left_type = "intron"
            if intron[0] <= novel_splice[1] - 1 <= intron[1] - 1:
                right_index = i
                right_type = "intron"

        print "left index " + str(left_index)
        print "left type " + str(left_type)
        print "right index " + str(right_index)
        print "right type " + str(right_type)
        #print   novel_splice
        if right_index is None or left_index is None:
            continue

        skipped_index = range(left_index + 1, right_index)
        #print "skipped_index " + str(skipped_index)
        i = 0
        temp_list = []
        ref_list = []

        while i < len(transcript.cds_list):

            if i in skipped_index:
                ref_list.append('\t'.join(
                    [transcript.chromosome, "test", "exon",
                     str(transcript.cds_list[i].location.start.position + 1),
                     str(transcript.cds_list[i].location.end.position), '.',
                     strand_symbol[transcript.cds_list[i].strand], '.',
                     """gene_id "{}\";""".format(gene_data.ID) + """ transcript_id "{}";""".format(
                         transcript.ID)]) + '\n')

            elif i != left_index and i != right_index:
                temp_list.append('\t'.join(
                    [transcript.chromosome, "test", "exon",
                     str(transcript.cds_list[i].location.start.position + 1),
                     str(transcript.cds_list[i].location.end.position), '.',
                     strand_symbol[transcript.cds_list[i].strand], '.',
                     """gene_id "{}\";""".format(gene_data.ID) + """ transcript_id "{}";""".format(
                         transcript.ID)]) + '\n')
                ref_list.append('\t'.join(
                    [transcript.chromosome, "test", "exon",
                     str(transcript.cds_list[i].location.start.position + 1),
                     str(transcript.cds_list[i].location.end.position), '.',
                     strand_symbol[transcript.cds_list[i].strand], '.',
                     """gene_id "{}\";""".format(gene_data.ID) + """ transcript_id "{}";""".format(
                         transcript.ID)]) + '\n')
            # left splice right splice in same intron/exon
            elif i == left_index == right_index:
                temp_list.append('\t'.join(
                    [transcript.chromosome, "test", "exon",
                     str(transcript.cds_list[i].location.start.position + 1),
                     str(novel_splice[0]), '.',
                     strand_symbol[transcript.cds_list[i].strand], '.',
                     """gene_id "{}\";""".format(gene_data.ID) + """ transcript_id "{}";""".format(
                         transcript.ID)]) + '\n')
                ref_list.append('\t'.join(
                    [transcript.chromosome, "test", "exon",
                     str(transcript.cds_list[i].location.start.position + 1),
                     str(transcript.cds_list[i].location.end.position), '.',
                     strand_symbol[transcript.cds_list[i].strand], '.',
                     """gene_id "{}\";""".format(gene_data.ID) + """ transcript_id "{}";""".format(
                         transcript.ID)]) + '\n')

                i += 1

                temp_list.append('\t'.join(
                    [transcript.chromosome, "test", "exon",
                     str(novel_splice[1] + 1),
                     str(transcript.cds_list[i].location.end.position), '.',
                     strand_symbol[transcript.cds_list[i].strand], '.',
                     """gene_id "{}\";""".format(gene_data.ID) + """ transcript_id "{}";""".format(
                         transcript.ID)]) + '\n')
                ref_list.append('\t'.join(
                    [transcript.chromosome, "test", "exon",
                     str(transcript.cds_list[i].location.start.position + 1),
                     str(transcript.cds_list[i].location.end.position), '.',
                     strand_symbol[transcript.cds_list[i].strand], '.',
                     """gene_id "{}\";""".format(gene_data.ID) + """ transcript_id "{}";""".format(
                         transcript.ID)]) + '\n')

            elif i == left_index:
                temp_list.append('\t'.join(
                    [transcript.chromosome, "test", "exon",
                     str(transcript.cds_list[i].location.start.position + 1),
                     str(novel_splice[0]), '.',
                     strand_symbol[transcript.cds_list[i].strand], '.',
                     """gene_id "{}\";""".format(gene_data.ID) + """ transcript_id "{}";""".format(
                         transcript.ID)]) + '\n')
                ref_list.append('\t'.join(
                    [transcript.chromosome, "test", "exon",
                     str(transcript.cds_list[i].location.start.position + 1),
                     str(transcript.cds_list[i].location.end.position), '.',
                     strand_symbol[transcript.cds_list[i].strand], '.',
                     """gene_id "{}\";""".format(gene_data.ID) + """ transcript_id "{}";""".format(
                         transcript.ID)]) + '\n')

            elif i == right_index:
                if right_type == "intron": #and novel_splice[1] != transcript.splice_list[i][1]:
                    ref_list.append('\t'.join(
                        [transcript.chromosome, "test", "exon",
                         str(transcript.cds_list[i].location.start.position + 1),
                         str(transcript.cds_list[i].location.end.position), '.',
                         strand_symbol[transcript.cds_list[i].strand], '.',
                         """gene_id "{}\";""".format(gene_data.ID) + """ transcript_id "{}";""".format(
                             transcript.ID)]) + '\n')
                    i+=1

                temp_list.append('\t'.join(
                        [transcript.chromosome, "test", "exon",
                         str(novel_splice[1]+1),
                         str(transcript.cds_list[i].location.end.position), '.',
                         strand_symbol[transcript.cds_list[i].strand], '.',
                         """gene_id "{}\";""".format(gene_data.ID) + """ transcript_id "{}";""".format(
                             transcript.ID)]) + '\n')
                ref_list.append('\t'.join(
                        [transcript.chromosome, "test", "exon",
                         str(transcript.cds_list[i].location.start.position + 1),
                         str(transcript.cds_list[i].location.end.position), '.',
                         strand_symbol[transcript.cds_list[i].strand], '.',
                         """gene_id "{}\";""".format(gene_data.ID) + """ transcript_id "{}";""".format(
                             transcript.ID)]) + '\n')

            i += 1


        if transcript.strand == -1:
            for rec in temp_list[::-1]:
                tumor_outf.write(rec)
            for rec in ref_list[::-1]:
                normal_outf.write(rec)
        else:
            for rec in temp_list:
                tumor_outf.write(rec)
            for rec in ref_list:
                normal_outf.write(rec)


def main():
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)-15s [%(processName)s.%(levelname)s] %(message)s')
    parser = argparse.ArgumentParser(description='Utilities for simulating reads')
    parser.add_argument('sample_id', type=str, nargs='?', help='provide sample id here')
    parser.add_argument('gff_file', type=str, nargs='?', help='provide gff3 file path here')
    parser.add_argument('outdir', type=str, nargs='?', help='provide outdir here')
    parser.add_argument("num_alterations", type=int, nargs='?', help='provide number of altered genes per chromosome here')
    args = parser.parse_args()

    gff_in_file = args.gff_file
    alteration_number_lim = args.num_alterations
    tumor_GTF_output = args.outdir + "tumor_GTF_output.gtf"
    normal_GTF_output = args.outdir + "normal_GTF_output.gtf"

    tumor_outf = open(tumor_GTF_output, 'w')
    normal_outf = open(normal_GTF_output, 'w')

    annotatedTransctipts = {}

    for i in range(1,23):
        chromosome = "chr" + str(i)
        logging.info("processing {}".format(chromosome))
        annotatedTransctipts[chromosome] = _find_annotated_transcripts(gff_in_file, chromosome)
        keys_to_choose = [key for key in annotatedTransctipts[chromosome] if
                          len(annotatedTransctipts[chromosome][key].transcripts) >= 8]
        genes_to_alter = np.random.choice(keys_to_choose, alteration_number_lim, replace=False)
        logging.info("selected genes: {}".format(','.join(genes_to_alter)))

        for gene in genes_to_alter:
            write_gtf(annotatedTransctipts[chromosome][gene], tumor_outf, normal_outf)

    tumor_outf.close()
    normal_outf.close()


if __name__ == '__main__':
    main()

