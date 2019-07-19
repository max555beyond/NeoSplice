import sys
import os
import itertools
import argparse
import pysam
import collections
import esgimpl
from BCBio import GFF
import bisect
from pyfaidx import Fasta
import subprocess
import re
import logging


class CigarOp:
    """Enumeration of CIGAR codes."""
    MATCH = 0
    INS = 1
    DEL = 2
    REF_SKIP = 3
    SOFT_CLIP = 4


class Exons(object):
    def __init__(self, chrnum, startposition, endposition, strand, ID):
        self.chrnum = chrnum
        self.startposition= startposition
        self.endposition = endposition
        self.strand = strand
        self.ID = ID


class Transcripts(object):
    def __init__(self, chrnum, startposition, endposition, strand, ID, genename, cds_list, start_codon, stop_codon,
                 splice_list):
        self.chrnum = chrnum
        self.startposition = startposition
        self.endposition = endposition
        self.strand = strand
        self.ID = ID
        self.genename = genename
        self.cds_list = cds_list
        self.start_codon = start_codon
        self.stop_codon = stop_codon
        self.splice_list = splice_list


def translate(seq, output_stop=False):
    gencode = {
      'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
      'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
      'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
      'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
      'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
      'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
      'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
      'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
      'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
      'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
      'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
      'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
      'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
      'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
      'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
      'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

    translate_seq = ''.join([gencode.get(seq[3 * j:3 * j + 3], 'X') for j in range(len(seq) // 3)])
    return translate_seq if ('_' not in translate_seq or output_stop) else translate_seq[:translate_seq.index('_')]


def contain_stop_codon(seq):
    if '_' in translate(seq, output_stop=True):
        return True
    else:
        return False


def run_netMHCpan(sample, chromosome, length, HLA_string, HLA_II_string, path, netMHCpan_path, netMHCIIpan_path):
    if 8 <= length <= 11:
        subprocess.call(
            '{} -f {}{}_outcome_peptide_{}_{}.fasta -BA -l {} -xls  -xlsfile {}{}_peptide_{}_{}.xls -a {}'.format(
                netMHCpan_path, path, sample, chromosome, length, length, path, sample, chromosome, length, HLA_string),
            shell=True, stdout=None)
    elif 15 <= length <= 24:
        subprocess.call(
            '{} -f {}{}_outcome_peptide_{}_{}.fasta -length {} -xls  -xlsfile {}{}_peptide_{}_{}.xls -a {}'.format(
                netMHCIIpan_path, path, sample, chromosome, length, length, path, sample, chromosome, length,
                HLA_II_string), shell=True, stdout=None)


def _open_bam(bam_name):
    return pysam.AlignmentFile(bam_name, 'rb')


def reverse_complement(seq):
    return u''.join(reversed([complement(base) for base in seq]))


def complement(base):
    if base.upper() == 'A':
        return 'T'
    elif base.upper() == 'T':
        return 'A'
    elif base.upper() == 'G':
        return 'C'
    elif base.upper() == 'C':
        return 'G'
    else:
        raise AssertionError("Unexpected value of nucleotide")


def _find_annotated_splices(gff_in_file, chromosome):
    in_file = gff_in_file
    in_handle = open(in_file)
    limit_info = dict(gff_id=[chromosome])
    transcripts_dic = {}
    for chromosome_inf in GFF.parse(in_handle, limit_info=limit_info):
        for gene in chromosome_inf.features:
            if "protein_coding" in gene.qualifiers["gene_type"]:
                for transcript in gene.sub_features:
                    if "protein_coding" in transcript.qualifiers["transcript_type"]:
                        start_codon = []
                        stop_codon = []
                        cds_list = []
                        splice_list = []
                        strand = "+" if transcript.strand == +1 else "-"

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

                        i = 0
                        while i < len(cds_list) - 1:
                            splice_list.append(
                                (cds_list[i].location.end.position, cds_list[i + 1].location.start.position))
                            i += 1

                        if len(start_codon) != 0:
                            transcripts_dic[transcript.id] = Transcripts(chromosome, transcript.location.start.position,
                                                                         transcript.location.end.position, strand,
                                                                         transcript.qualifiers["transcript_id"],
                                                                         gene.qualifiers["gene_name"], cds_list,
                                                                         start_codon, stop_codon, splice_list)
    return transcripts_dic


def _get_structrual_edge_read(chromosome, edge, bam, direction):
    if direction == -1:
        edge = edge[:2][::-1] + edge[2:]

    read_set = set()
    # add 1 and subtract 1 to deal with ins egdes with same start and end pos
    for read in bam.fetch(chromosome, edge[0]-1, edge[1]+1):
        read_pos = 0
        current_interval_start = read.reference_start
        current_interval_end = read.reference_start

        if not read.cigartuples:
            continue

        for op, count in read.cigartuples:
            if op == CigarOp.MATCH:
                current_interval_start = current_interval_end
                current_interval_end = current_interval_end + count
                read_pos += count
            elif op == CigarOp.DEL or op == CigarOp.REF_SKIP:
                current_interval_start = current_interval_end
                current_interval_end = current_interval_end + count
            elif op == CigarOp.INS or op == CigarOp.SOFT_CLIP:
                current_interval_start = current_interval_end
                read_pos += count
            else:
                logging.warn('Unexpected cigar op {0}'.format(op))

            if edge[2] == 'splice' and current_interval_start == edge[0] and current_interval_end == edge[1] and \
                    op == CigarOp.REF_SKIP:
                read_set.add(read.query_name)
                break
            if edge[2] == 'del' and current_interval_start == edge[0] and current_interval_end == edge[1] and \
                    op == CigarOp.DEL:
                read_set.add(read.query_name)
                break
            if edge[0] == edge[1] and current_interval_start == edge[0] and current_interval_end == edge[1] and \
                    op == CigarOp.INS and read.query_sequence[read_pos - count:read_pos].upper() == edge[2]:
                read_set.add(read.query_name)
                break
    return read_set


def _get_exon_edge_read(chromosome, edge, bam, direction, genome):
    if direction == -1:
        edge = edge[:2][::-1] + edge[2:]

    read_pos = collections.defaultdict(set)
    for column in bam.pileup(chromosome, edge[0], edge[1], truncate=True):
        for read in column.pileups:
            if not read.is_refskip and not read.is_del and read.alignment.query_sequence[read.query_position] == \
                    genome[chromosome][column.reference_pos].seq:
                read_pos[column.reference_pos].add(read.alignment.query_name)

    return set.union(*read_pos.values()) if read_pos else set()


def _get_snp_read(chromosome, edge, bam, direction):
    if direction == -1:
        edge = edge[:2][::-1] + edge[2:]

    reads = set()
    for column in bam.pileup(chromosome, edge[0], edge[1], truncate=True):
        for read in column.pileups:
            if not read.is_refskip and not read.is_del and read.alignment.query_sequence[read.query_position] == edge[
                    2]:
                reads.add(read.alignment.query_name)

    return reads


def get_read(chromosome, edge, bam, direction, genome):
    if edge[2] == 'exon':
        return _get_exon_edge_read(chromosome, edge, bam, direction, genome)
    # SNP edge
    elif len(edge[2]) == 1:
        return _get_snp_read(chromosome, edge, bam, direction)
    else:
        return _get_structrual_edge_read(chromosome, edge, bam, direction)


def find_path_annotated(graph_outer, start_node_outer, target_node, annotated_splices, direction):
    possible_path_list = []
    paths = []
    possible_seq_list = []
    seqs = []

    if start_node_outer not in graph_outer or target_node not in graph_outer:
        return possible_path_list, possible_seq_list

    if direction == +1 and start_node_outer >= target_node or direction == -1 and start_node_outer <= target_node:
        return possible_path_list, possible_seq_list

    def path(graph, start_node, target_node, annotated_splices, direction):
        if start_node == target_node:
            if direction == -1:
                consecutive_path = tuple(paths[::-1])
                possible_path_list.append(consecutive_path)
                possible_seq_list.append(''.join(seqs))

            else:
                consecutive_path = tuple(paths)
                possible_path_list.append(consecutive_path)
                possible_seq_list.append(''.join(seqs))

            return

        for edge in graph.edges(start_node, data=True, keys=True):
            if edge[2] == 'splice':
                if (direction == +1 and (edge[0], edge[1]) in annotated_splices and edge[1] <= target_node) or (
                        direction == -1 and (edge[1], edge[0]) in annotated_splices and edge[1] >= target_node):
                    if direction == -1:
                        paths.append(edge[:2][::-1] + (edge[2], None))
                    else:
                        paths.append(edge[:3] + (None,))

                    path(graph, edge[1], target_node, annotated_splices, direction)
                    paths.pop()
                    return

        for edge in graph.edges(start_node, data=True, keys=True):
            # avoid annotated splices not present in the graph to be skipped
            if edge[2] == 'exon':
                if direction == +1 and not any(
                        edge[0] <= splice[0] < edge[1] for splice in annotated_splices if splice[1] <= target_node) and\
                        edge[1] <= target_node:
                    paths.append(edge[:3] + (None,))
                    seqs.append(edge[3].get("seq", ''))
                elif direction == -1 and not any(
                        edge[1] < splice[1] <= edge[0] for splice in annotated_splices if splice[0] >= target_node) and\
                        edge[1] >= target_node:
                    paths.append(edge[:2][::-1] + (edge[2], None))
                    seqs.append(reverse_complement(edge[3].get('seq', '')))
                else:
                    continue

                path(graph, edge[1], target_node, annotated_splices,direction)
                paths.pop()
                seqs.pop()
                return

    path(graph_outer, start_node_outer, target_node, annotated_splices, direction)
    return possible_path_list,possible_seq_list


def find_path_dfs(chromosome, bam, graph_outer, start_node_outer, read_set, direction, strand, genome, min_coverage,
                  upstream_limit=None):
    possible_path_list = []
    paths = []
    possible_seq_list = []
    seqs = []
    possible_path_readsets = []
    visited_ins_nodes = set()

    if upstream_limit and (
            start_node_outer >= upstream_limit and direction == +1 or start_node_outer <= upstream_limit
            and direction == -1):
        return possible_path_list, possible_seq_list, possible_path_readsets

    def path(graph, start_node, read_set, direction, strand, visited_ins_nodes):
        continue_indication = False
        edges = [edge for edge in graph.edges(start_node, data=True, keys=True) if
                 edge[1] not in visited_ins_nodes and len(
                     set.intersection(get_read(chromosome, edge, bam, direction, genome), read_set)) >= min_coverage]
        ins_edges = [edge for edge in edges if edge[0] == edge[1]]

        if ins_edges:
            ins_edges_reads = set.intersection(
                set.union(*[get_read(chromosome, edge, bam, direction, genome) for edge in ins_edges]), read_set)

            # determine if all reads containing non-ins edges also contain insertions, if so must go through is edges first
            for edge in edges[:]:
                if edge[0] != edge[1] and set.intersection(get_read(chromosome, edge, bam, direction, genome),
                                                           read_set).issubset(ins_edges_reads):
                    edges.remove(edge)

        for edge in edges:
            if upstream_limit and start_node == upstream_limit and not ins_edges:
                continue_indication = False
                break
            # path with splice or del edge extend past start codon will not be included in the result
            if upstream_limit and (
                  (edge[1] < upstream_limit) and (direction == -1) or (edge[1] > upstream_limit) and (direction == +1)):
                continue_indication = True
                continue

            new_read_set = set.intersection(get_read(chromosome, edge, bam, direction, genome), read_set)
            # Splice and deletion edges
            if edge[2] == "splice" or edge[2] == "del":
                continue_indication = True
                if direction == -1:
                    paths.append(edge[:2][::-1] + (edge[2], len(new_read_set)))
                else:
                    paths.append(edge[:3] + (len(new_read_set),))

                path(graph, edge[1], new_read_set, direction, strand, visited_ins_nodes)
                paths.pop()
            # Exon edges
            elif edge[2] == 'exon':
                continue_indication = True
                seqs.append(edge[3]["seq"])
                if direction == -1:
                    paths.append(edge[:2][::-1] + (edge[2], len(new_read_set)))
                else:
                    paths.append(edge[:3] + (len(new_read_set),))

                path(graph, edge[1], new_read_set, direction, strand, visited_ins_nodes)
                seqs.pop()
                paths.pop()
            # Insertion edges
            elif edge[0] == edge[1]:
                continue_indication = True
                seqs.append(edge[3]["seq"])
                paths.append(edge[:3] + (len(new_read_set),))
                visited_ins_nodes.add(edge[1])
                path(graph, edge[1], new_read_set, direction, strand, visited_ins_nodes)
                visited_ins_nodes.remove(edge[1])
                seqs.pop()
                paths.pop()

            # SNP edges
            elif len(edge[2]) == 1:
                continue_indication = True
                seqs.append(edge[3]["seq"])
                if direction == -1:
                    paths.append(edge[:2][::-1] + (edge[2], len(new_read_set)))
                else:
                    paths.append(edge[:3] + (len(new_read_set),))

                path(graph, edge[1], new_read_set, direction, strand, visited_ins_nodes)
                seqs.pop()
                paths.pop()

        if not continue_indication and len(paths):
            possible_path_readsets.append(read_set)
            if direction == +1:
                consecutive_path = tuple(paths)
                possible_path_list.append(consecutive_path)
                if strand == '+':
                    possible_seq_list.append(''.join(seqs))
                elif strand == '-':
                    possible_seq_list.append(''.join([reverse_complement(seq) for seq in seqs[::-1]]))
            else:
                consecutive_path = tuple(paths[::-1])
                possible_path_list.append(consecutive_path)
                if strand == '+':
                    possible_seq_list.append(''.join(seqs[::-1]))
                elif strand == '-':
                    possible_seq_list.append(''.join([reverse_complement(seq) for seq in seqs]))

    path(graph_outer, start_node_outer, read_set, direction, strand, visited_ins_nodes)
    return possible_path_list, possible_seq_list, possible_path_readsets


def _get_exon_edge(se, loc):
    i = bisect.bisect_right([exon[0] for exon in se], loc)
    if i:
        s, t, _ = se[i - 1]
        if t > loc:
            return se[i - 1]
        else:
            return None
    else:
        return None


def validate_path_edges(graph, graph_paths):
    old_graph_paths = collections.deque(graph_paths)
    new_graph_paths = []

    def dfs(graph, node, graph_paths, new_graph_paths):
        if len(graph_paths) == 0:
            return True

        for edge in graph.edges(node, keys=True, data=True):
            if edge[0] == graph_paths[0][0] and ((
                    edge[1] <= graph_paths[0][1] and edge[2] == graph_paths[0][2] == 'exon') or (
                    edge[1] == graph_paths[0][1] and edge[2] == graph_paths[0][2] != 'exon')):
                prev_edge = graph_paths.popleft()
                if edge[1] < prev_edge[1]:
                    graph_paths.appendleft((edge[1], prev_edge[1], prev_edge[2]))
                if "seq" in edge[3]:
                    new_graph_paths.append((edge[0], edge[1], edge[2], edge[3]["weight"], edge[3]["seq"]))
                else:
                    new_graph_paths.append((edge[0], edge[1], edge[2], edge[3]["weight"]))

                return dfs(graph, edge[1], graph_paths, new_graph_paths)
        else:
            return False
    # K-mer SNP edge not significant
    if old_graph_paths[0][0] not in graph:
        return []

    if dfs(graph, old_graph_paths[0][0], old_graph_paths, new_graph_paths):
        return new_graph_paths
    else:
        return []


def _find_Kmer_snps(ref_pos, sequence, ref_sequence):
    return [(ref_pos + i, nuc, ref_nuc) for i, (nuc, ref_nuc) in enumerate(zip(sequence, ref_sequence)) if
            nuc != ref_nuc]


def _map_contig_to_splice_graph(start_exon, end_exon, Kmer, chromosome, genome):
    graph_paths = []
    ref_pos = Kmer.reference_start
    read_pos = 0

    for op, count in Kmer.cigartuples:
        if op == CigarOp.MATCH:
            snps = _find_Kmer_snps(ref_pos, Kmer.query_sequence[read_pos:read_pos + count].upper(),
                                  genome[chromosome][ref_pos:ref_pos + count].seq.upper())
            if snps:
                end_pos = ref_pos + count
                for snp_pos, snp, _ in snps:
                    if ref_pos != snp_pos:
                        graph_paths.append((ref_pos, snp_pos, "exon"))
                    graph_paths.append((snp_pos, snp_pos + 1, snp))
                    ref_pos = snp_pos + 1
                if ref_pos != end_pos:
                    graph_paths.append((ref_pos, end_pos, "exon"))
                    ref_pos = end_pos
            else:
                graph_paths.append((ref_pos, ref_pos + count, "exon"))
                ref_pos += count
            read_pos += count
        elif op == CigarOp.DEL or op == CigarOp.REF_SKIP:
            if op == CigarOp.DEL:
                graph_paths.append((ref_pos, ref_pos + count, "del"))
            elif op == CigarOp.REF_SKIP:
                graph_paths.append((ref_pos, ref_pos + count, "splice"))
            ref_pos += count
        elif op == CigarOp.INS or op == CigarOp.SOFT_CLIP:
            if op == CigarOp.INS:
                graph_paths.append((ref_pos, ref_pos, Kmer.query_sequence[read_pos:read_pos + count].upper()))
            read_pos += count
        else:
            logging.warn('Unexpected cigar op {0}'.format(op))

    if graph_paths[0][2] == "exon":
        graph_paths[0] = (start_exon[0], graph_paths[0][1], graph_paths[0][2])
    if graph_paths[-1][2] == "exon":
        graph_paths[-1] = (graph_paths[-1][0], end_exon[1], graph_paths[-1][2])
    return graph_paths


def retrieve_splice_edge(read):
    ref_pos = read.reference_start
    read_pos = 0
    splices = []
    for op, count in read.cigartuples:
        if op == CigarOp.MATCH:
            ref_pos += count
            read_pos += count
        elif op == CigarOp.INS:
            read_pos += count
        elif op == CigarOp.DEL:
            # Deletion edge extends from first deleted base to next retained base
            ref_pos += count
        elif op == CigarOp.REF_SKIP:
            # Splice edge extends from first excised base to next retained base
            start = ref_pos
            end = ref_pos + count
            splices.append((start, end, "splice"))
            ref_pos += count
        elif op == CigarOp.SOFT_CLIP:
            read_pos += count
        else:
            logging.warn('Unexpected cigar op {0}'.format(op))

    return splices


def get_supported_seq(upstream_seq, downstream_seq, unique_path, frame_before, novel_splice, strand, length):
    if strand == "+":
        selected_upstream_seq = (upstream_seq + ''.join([edge[4] for edge in unique_path if
                            edge[2] != "splice" and edge[2] != "del" and edge[0] <
                            novel_splice[0]]))
        # full upstream codon length
        full_frame_len = min(length - 1, (len(selected_upstream_seq) - (3 - frame_before) % 3) // 3)
        selected_upstream_seq = selected_upstream_seq[-3 * full_frame_len - (3 - frame_before) % 3:]
        selected_downstream_seq = ''.join([edge[4] for edge in unique_path if
                            edge[2] != "splice" and edge[2] != "del" and edge[0] > novel_splice[
                            0]]) + downstream_seq
        full_frame_len_downstream = min((len(selected_downstream_seq) - frame_before) // 3, length - 1)
        selected_downstream_seq = selected_downstream_seq[:frame_before + full_frame_len_downstream * 3]
    else:
        selected_upstream_seq = (upstream_seq + reverse_complement(''.join([edge[4] for edge in unique_path if
                            edge[2] != "splice" and edge[2] != "del" and edge[1] >
                            novel_splice[1]])))
        full_frame_len = min(length - 1, (len(selected_upstream_seq) - (3 - frame_before) % 3) / 3)
        selected_upstream_seq = selected_upstream_seq[-3 * full_frame_len - (3 - frame_before) % 3:]
        selected_downstream_seq = reverse_complement(''.join([edge[4] for edge in unique_path if
                                    edge[2] != "splice" and edge[2] != "del" and edge[1] < novel_splice[
                                        1]])) + downstream_seq
        full_frame_len_downstream = min((len(selected_downstream_seq) - frame_before) // 3, length - 1)
        selected_downstream_seq = selected_downstream_seq[:frame_before + full_frame_len_downstream * 3]

    return selected_upstream_seq + selected_downstream_seq, len(selected_upstream_seq), len(selected_downstream_seq)


def get_supported_path(upstream_len, downstream_len, novel_splice, strand, graph_path):
    if strand == "-":
        graph_path = graph_path[::-1]
        upstream_len, downstream_len = downstream_len, upstream_len

    index = [edge[:3] for edge in graph_path].index(novel_splice)
    selected_upstream_path = []
    selected_downstream_path = []

    def get_len(edge):
        if edge[0] == edge[1]:
            return len(edge[2])
        elif edge[2] == "splice" or edge[2] == "del":
            return 0
        else:
            return edge[1] - edge[0]

    def get_edge(edge, remaining_len, direction):
        edge_len = get_len(edge)
        if edge[0] == edge[1]:
            if direction == +1:
                return edge[0], edge[1], edge[2][:min(remaining_len, edge_len)], edge[3]
            else:
                return edge[0], edge[1], edge[2][-min(remaining_len, edge_len):], edge[3]

        elif edge[2] == "splice" or edge[2] == "del":
            return edge

        else:
            if direction == +1:
                return edge[0], min(edge[0] + remaining_len, edge[1]), edge[2], edge[3]
            else:
                return max(edge[0], edge[1] - remaining_len), edge[1], edge[2], edge[3]

    i = index - 1
    while upstream_len > 0 and i >= 0:
        selected_upstream_path.append(get_edge(graph_path[i], upstream_len, -1))
        upstream_len -= get_len(graph_path[i])
        i -= 1

    i = index + 1
    while downstream_len > 0 and i < len(graph_path):
        selected_downstream_path.append(get_edge(graph_path[i], downstream_len, +1))
        downstream_len -= get_len(graph_path[i])
        i += 1

    selected_path = selected_upstream_path[::-1] + [graph_path[index]] + selected_downstream_path
    return tuple(selected_path) if strand == "+" else tuple(selected_path[::-1])


def output_peptide(mut_peptide, mut_sequence, result_file, fasta_file, length, upstream_len, generated_peptides,
                   novel_splice, selected_path, strand, chromosome, full_path, full_seq, gene_name, peptide_count):
    if len(mut_peptide) < length:
        return
    for i in range(0, len(mut_peptide) - length + 1):
        peptide_path = get_supported_path(upstream_len - 3 * i, length * 3 - upstream_len + 3 * i, novel_splice, strand,
                                          selected_path)
        if (mut_peptide[i:i+length], novel_splice) not in generated_peptides:
            fasta_file.write(">peptide{}".format(peptide_count) + '\n' + mut_peptide[i:i+length] + '\n')
            result_file.write(mut_peptide[i:i + length] + '\t' + mut_sequence[3 * i:3 * (
                    i + length)] + '\t' + chromosome + '\t' + str(peptide_path) + '\t' + str(full_path) + '\t' +
                    full_seq + '\t' + str(novel_splice) + '\t' + strand + '\t' + gene_name + '\n')
            generated_peptides.add((mut_peptide[i:i+length], novel_splice))
            peptide_count += 1

    return peptide_count


def combine_table(sample, neoantigen_path, length, chromosome):
    MHC_file = open(neoantigen_path + '{}_peptide_{}_{}.xls'.format(sample, chromosome, length), "r")
    dat_MHC = [row for row in MHC_file]
    MHC_file.close()
    HLA_alleles = re.split('\t+', dat_MHC[0].strip())
    if length < 15:
        NMposition = 6
    else:
        NMposition = 4
    for HLA in HLA_alleles:
        file_MHC = open(neoantigen_path + '{}_peptide_MHC_{}_{}.xls'.format(sample, chromosome, length), "w")
        file_MHC.write(('{}Nm'.format(HLA) + '\t' + '{}binding property'.format(HLA)) + '\n')
        current_line = 2

        while current_line < len(dat_MHC):
            row_split = dat_MHC[current_line].strip().split('\t')
            if 0 < float(row_split[NMposition]) <= 50:
                file_MHC.write(row_split[NMposition] + '\t' + 'strong binders' + '\n')
            elif 50 < float(row_split[NMposition]) <= 150:
                file_MHC.write(row_split[NMposition] + '\t' + 'moderate binders' + '\n')
            elif 150 < float(row_split[NMposition]) <= 500:
                file_MHC.write(row_split[NMposition] + '\t' + 'weak binders' + '\n')
            else:
                file_MHC.write(row_split[NMposition] + '\t' + 'non-binders' + '\n')
            current_line += 1
        file_MHC.close()
        NMposition += 5
        subprocess.call(
            'paste {}{}_outcome_peptide_{}_{}.txt {}{}_peptide_MHC_{}_{}.xls > {}{}temp_outcome_peptide_{}_{}.txt'.format(
                neoantigen_path, sample, chromosome, length, neoantigen_path, sample, chromosome, length,
                neoantigen_path, sample,  chromosome, length), shell=True)
        subprocess.call(
            'paste {}{}_outcome_peptide_{}_{}.txt {}{}temp_outcome_peptide_{}_{}.txt > {}{}_outcome_peptide_{}_{}.txt'.format(
                neoantigen_path, sample, chromosome, length, neoantigen_path, sample, chromosome, length,
                neoantigen_path, sample, chromosome, length), shell=True)


def main():
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)-15s [%(processName)s.%(levelname)s] %(message)s')

    parser = argparse.ArgumentParser(description='Utilities for creating and working with splice graphs.')
    parser.add_argument('sample', type=str, nargs='?', help='The sample name')
    parser.add_argument('chromosome', type=str, nargs='?', help='The chromosome')
    parser.add_argument('bam_file', type=str, nargs='?', help='provide RNA-seq bam file path here')
    parser.add_argument('gff_file', type=str, nargs='?', help='The GFF file path')
    parser.add_argument('genome_fasta', type=str, nargs='?', help='The genome FASTA file path')
    parser.add_argument('kmer_bam', type=str, nargs='?', help='The bam file storing kmer info')
    parser.add_argument('splice_graph', type=str, nargs='?', help='The path to splice graph')
    parser.add_argument('tumor_junction_file', type=str, nargs='?', help='The file storing tumor junctions')
    parser.add_argument('normal_junction_file', type=str, nargs='?', help='The file storing normal junctions')
    parser.add_argument('length', type=str, nargs='?', help='The output peptide length')
    parser.add_argument('transcript_min_coverage', type=str, nargs='?', help='The minimum transcript coverage')
    parser.add_argument('netMHCpan_path', type=str, nargs='?', help='The netMHCpan-4.0 path')
    parser.add_argument('netMHCIIpan_path', type=str, nargs='?', help='The netMHCIIpan-3.2 path')
    parser.add_argument('outdir', type=str, nargs='?', help='The output directory')

    args = parser.parse_args()
    HLA_string = "HLA-A02:01"
    HLA_II_string = "DRB1_1601"

    sample = args.sample
    chromosome = args.chromosome
    bam_file = args.bam_file
    gff_in_file = args.gff_file
    tumor_junction_file = args.tumor_junction_file
    normal_junction_file = args.normal_junction_file
    splice_graph_dir = args.splice_graph_dir
    length = int(args.length)
    min_coverage = int(args.transcript_min_coverage)
    peptide_count = 0

    netMHCpan_path = args.netMHCpan_path
    netMHCIIpan_path = args.netMHCIIpan_path
    neoantigen_path = os.path.join(args.outdir+"neoantigen_result/", sample + '/')
    genome = Fasta(args.genome_fasta)

    if not os.path.isdir(neoantigen_path) and not os.path.exists(neoantigen_path):
        os.makedirs(neoantigen_path)

    output_file = open(neoantigen_path + "{}_outcome_peptide_{}_{}.txt".format(sample, chromosome, length), 'w')
    output_file.write("Variant_peptide_sequence\tDNA_sequence\tChromosome\tPeptide_graph_path\tFull_graph_path\tFull_graph_seq\tTumor_splice\tStrand\tGene\n")
    fasta_file = open(neoantigen_path + "{}_outcome_peptide_{}_{}.fasta".format(sample, chromosome, length), 'w')

    generated_peptides = set()
    normal_set = set()
    tumor_set = set()

    with open(normal_junction_file) as f:
        for line in f:
            line_split = line.strip().split("\t")
            if int(line_split[6]) >= 3:
                normal_set.add((line_split[0], int(line_split[1]) - 1, int(line_split[2])))

    with open(tumor_junction_file) as f:
        for line in f:
            line_split = line.strip().split("\t")
            if (line_split[0], int(line_split[1]) - 1, int(line_split[2])) not in normal_set and int(line_split[6]) >= 20:
                tumor_set.add((line_split[0], int(line_split[1]) - 1, int(line_split[2])))

    splice_graph = esgimpl.EsgImpl()
    splice_graph.load_from_file(args.splice_graph)

    bam = _open_bam(bam_file)
    annotated_transcripts = _find_annotated_splices(gff_in_file, chromosome)

    kmer_dat = _open_bam(args.kmer_bam)

    unique_path_edges = collections.defaultdict(set)
    unique_path_splices = collections.defaultdict(set)
    validated_path = {}
    splice_subgraph = {}

    logging.info("reading kmers")

    for kmer in kmer_dat.fetch(chromosome):
        novel_splices = retrieve_splice_edge(kmer)
        for novel_splice in novel_splices:
            # skip if not a novel splice
            if (chromosome,) + novel_splice[:-1] not in tumor_set:
                continue

            if novel_splice in splice_subgraph:
                gene = splice_subgraph[novel_splice]
            else:
                _, gene = splice_graph._to_gene(novel_splice)
                splice_subgraph[novel_splice] = gene

            exons = [edge for edge in gene.edges(keys=True) if edge[2] == "exon"]
            se = sorted(exons, key=lambda x: x[0])
            start_exon = _get_exon_edge(se, kmer.reference_start)
            end_exon = _get_exon_edge(se, kmer.reference_end - 1)

            if start_exon is None or end_exon is None:
                continue

            path_edges = _map_contig_to_splice_graph(start_exon, end_exon, kmer, chromosome, genome)

            if tuple(path_edges) not in validated_path:
                validated_path[tuple(path_edges)] = validate_path_edges(gene, path_edges)

            path_edges = validated_path[tuple(path_edges)]

            if not path_edges:
                continue

            unique_path_edges[tuple(path_edges)].add(kmer.query_name)
            unique_path_splices[tuple(path_edges)].add(novel_splice)

    logging.info("found {} unique kmer paths".format(len(unique_path_edges)))

    for kmer_num, unique_path in enumerate(unique_path_edges):
        logging.info("processing kmer {}".format(kmer_num))

        if len(unique_path_edges[unique_path]) < min_coverage:
            continue

        novel_splices = unique_path_splices[unique_path]
        read_set = unique_path_edges[unique_path]
        unique_path = tuple([edge[:3] + (len(unique_path_edges[unique_path]),) + edge[4:] for edge in unique_path])
        for novel_splice in novel_splices:
            logging.info("processing novel splice" + str(novel_splice))

            gene = splice_subgraph[novel_splice]
            kmer_seq = ''.join([edge[4] for edge in unique_path if
                        edge[2] != "splice" and edge[2] != "del"])

            possible_transcripts = []
            for transcript in annotated_transcripts:
                if unique_path[0][0] >= annotated_transcripts[transcript].startposition and unique_path[-1][1] <= \
                        annotated_transcripts[transcript].endposition:
                    possible_transcripts.append(transcript)

            for transcript in possible_transcripts:
                full_upstream_paths = []
                full_upstream_seqs = []
                supported_upstream_paths = []
                supported_upstream_seqs = []

                logging.info("using annotated transcript {}".format(transcript))

                if annotated_transcripts[transcript].strand == "+" and annotated_transcripts[transcript].start_codon[
                    -1].location.end.position < novel_splice[0]:
                    upstream_paths, upstream_seqs, upstream_read_set = find_path_dfs(chromosome, bam, gene.reverse(),
                                                                                        unique_path[0][0],
                                                                                        read_set, -1, '+', genome,
                                                                                     min_coverage,
                                                                        annotated_transcripts[transcript].start_codon[
                                                                        -1].location.end.position)

                    for upstream_path, upstream_seq in zip(upstream_paths, upstream_seqs):
                        if upstream_path[0][0] == annotated_transcripts[transcript].start_codon[
                                -1].location.end.position:
                            full_upstream_paths.append(upstream_path)
                            full_upstream_seqs.append(upstream_seq)
                            supported_upstream_paths.append(upstream_path)
                            supported_upstream_seqs.append(upstream_seq)
                        else:
                            annotated_exon_path, annotated_seq = find_path_annotated(gene, annotated_transcripts[
                                transcript].start_codon[-1].location.end.position, upstream_path[0][0],
                                                                           annotated_transcripts[
                                                                               transcript].splice_list, +1)
                            if len(annotated_exon_path):
                                full_upstream_paths.append(annotated_exon_path[0] + upstream_path)
                                full_upstream_seqs.append(annotated_seq[0] + upstream_seq)
                                supported_upstream_paths.append(upstream_path)
                                supported_upstream_seqs.append(upstream_seq)
                    if not len(upstream_paths):
                        annotated_exon_path, annotated_seq = find_path_annotated(gene, annotated_transcripts[
                            transcript].start_codon[-1].location.end.position, unique_path[0][0],
                                                                           annotated_transcripts[
                                                                               transcript].splice_list, +1)
                        if len(annotated_exon_path):
                            full_upstream_paths.append(annotated_exon_path[0])
                            full_upstream_seqs.append(annotated_seq[0])
                            supported_upstream_paths.append(None)
                            supported_upstream_seqs.append("")

                    if len(full_upstream_paths) or unique_path[0][0] <= annotated_transcripts[transcript].start_codon[
                                                                     -1].location.end.position:
                        full_downstream_paths, full_downstream_seqs, _ = find_path_dfs(chromosome, bam, gene,
                                                                                       unique_path[-1][1],
                                                                                       read_set, +1, '+', genome,
                                                                                       min_coverage)

                        if len(full_upstream_paths) and len(full_downstream_paths):
                            for combination in itertools.product(range(len(full_upstream_paths)),
                                                                 range(len(full_downstream_paths))):

                                graph_path = full_upstream_paths[combination[0]] + tuple(
                                    [edge[:4] for edge in unique_path]) + full_downstream_paths[combination[1]]
                                graph_seq = full_upstream_seqs[combination[0]] + kmer_seq + full_downstream_seqs[
                                    combination[1]]
                                upstream_seq_temp = full_upstream_seqs[combination[0]] + ''.join(
                                    [edge[4] for edge in unique_path if
                                     edge[2] != "splice" and edge[2] != "del" and edge[0] < novel_splice[0]])

                                if contain_stop_codon(upstream_seq_temp):
                                    continue

                                frame_before = -len(upstream_seq_temp) % 3

                                supported_seq, upstream_len, downstream_len = get_supported_seq(
                                    supported_upstream_seqs[combination[0]], full_downstream_seqs[combination[1]],
                                    unique_path, frame_before, novel_splice, '+', length)
                                supported_path = get_supported_path(upstream_len, downstream_len, novel_splice, '+',
                                                                    graph_path)
                                mut_peptide = translate(supported_seq)
                                count = output_peptide(mut_peptide, supported_seq, output_file, fasta_file, length,
                                                       upstream_len, generated_peptides, novel_splice, supported_path,
                                                       "+", chromosome, graph_path, graph_seq,
                                                       annotated_transcripts[transcript].genename[0], peptide_count)

                                if count:
                                    peptide_count = count

                        elif len(full_upstream_paths) and not len(full_downstream_paths):
                            for full_upstream_path, full_upstream_seq, supported_upstream_seq in zip(
                                    full_upstream_paths, full_upstream_seqs, supported_upstream_seqs):
                                graph_path = full_upstream_path + tuple([edge[:4] for edge in unique_path])
                                graph_seq = full_upstream_seq + kmer_seq
                                upstream_seq_temp = full_upstream_seq + ''.join(
                                    [edge[4] for edge in unique_path if
                                     edge[2] != "splice" and edge[2] != "del" and edge[0] < novel_splice[0]])

                                if contain_stop_codon(upstream_seq_temp):
                                    continue

                                frame_before = -len(upstream_seq_temp) % 3

                                supported_seq, upstream_len, downstream_len = get_supported_seq(supported_upstream_seq,
                                                                                                '', unique_path,
                                                                                                frame_before,
                                                                                                novel_splice, '+',
                                                                                                length)
                                supported_path = get_supported_path(upstream_len, downstream_len, novel_splice, '+',
                                                                    graph_path)
                                mut_peptide = translate(supported_seq)
                                count = output_peptide(mut_peptide, supported_seq, output_file, fasta_file, length,
                                                       upstream_len,generated_peptides, novel_splice, supported_path,
                                                       "+", chromosome, graph_path, graph_seq,
                                                       annotated_transcripts[transcript].genename[0],peptide_count)

                                if count:
                                    peptide_count = count

                        elif unique_path[0][0] <= annotated_transcripts[transcript].start_codon[
                                -1].location.end.position and len(full_downstream_paths):
                            for full_downstream_path, full_downstream_seq in zip(full_downstream_paths,
                                                                                 full_downstream_seqs):
                                graph_path = tuple([edge[:4] for edge in unique_path if
                                     edge[0] >= annotated_transcripts[transcript].start_codon[
                                         -1].location.end.position]) + full_downstream_path
                                graph_seq = ''.join([edge[4] for edge in unique_path if
                                                    edge[2] != "splice" and edge[2] != "del" and
                                                    edge[0] >= annotated_transcripts[transcript].start_codon[
                                        -1].location.end.position]) + full_downstream_seq
                                upstream_seq_temp = ''.join([edge[4] for edge in unique_path if
                                                    edge[2] != "splice" and edge[2] != "del" and
                                                     annotated_transcripts[transcript].start_codon[
                                        -1].location.end.position <= edge[0] < novel_splice[0]])
                                # make sure start codon is not skipped by splice or deletion and no stop codon before splice
                                if graph_path[0][0] != annotated_transcripts[transcript].start_codon[
                                        -1].location.end.position or contain_stop_codon(upstream_seq_temp):
                                    continue

                                frame_before = -len(upstream_seq_temp) % 3

                                supported_seq, upstream_len, downstream_len = get_supported_seq('',
                                            full_downstream_seq, unique_path, frame_before, novel_splice, '+', length)
                                supported_path = get_supported_path(upstream_len, downstream_len, novel_splice, '+',
                                                                    graph_path)
                                mut_peptide = translate(supported_seq)
                                count = output_peptide(mut_peptide, supported_seq, output_file, fasta_file, length,
                                                       upstream_len, generated_peptides, novel_splice, supported_path,
                                                       "+", chromosome, graph_path, graph_seq,
                                                       annotated_transcripts[transcript].genename[0],peptide_count)

                                if count:
                                    peptide_count = count


                        elif unique_path[0][0] <= annotated_transcripts[transcript].start_codon[
                                    -1].location.end.position and not len(full_upstream_paths) and not len(
                                    full_downstream_paths):
                            graph_path = tuple([edge[:4] for edge in unique_path if
                                 edge[0] >= annotated_transcripts[transcript].start_codon[-1].location.end.position])
                            graph_seq = ''.join([edge[4] for edge in unique_path if
                                                 edge[2] != "splice" and edge[2] != "del" and edge[0] >=
                                                 annotated_transcripts[transcript].start_codon[
                                                     -1].location.end.position])
                            upstream_seq_temp = ''.join([edge[4] for edge in unique_path if
                                                         edge[2] != "splice" and edge[2] != "del" and
                                                         annotated_transcripts[transcript].start_codon[
                                                             -1].location.end.position <= edge[0] < novel_splice[0]])

                            if graph_path[0][0] != annotated_transcripts[transcript].start_codon[
                                    -1].location.end.position or contain_stop_codon(upstream_seq_temp):
                                continue

                            frame_before = -len(upstream_seq_temp) % 3

                            supported_seq, upstream_len, downstream_len = get_supported_seq('', '', unique_path,
                                                                                            frame_before, novel_splice,
                                                                                            '+', length)
                            supported_path = get_supported_path(upstream_len, downstream_len, novel_splice, '+',
                                                                graph_path)
                            mut_peptide = translate(supported_seq)
                            count = output_peptide(mut_peptide, supported_seq, output_file, fasta_file, length,
                                                   upstream_len, generated_peptides, novel_splice, supported_path,
                                                   "+", chromosome, graph_path, graph_seq,
                                                   annotated_transcripts[transcript].genename[0], peptide_count)

                            if count:
                                peptide_count = count

                elif annotated_transcripts[transcript].strand == "-" and annotated_transcripts[transcript].start_codon[
                    0].location.start.position > novel_splice[1]:
                    upstream_paths, upstream_seqs, upstream_read_set = find_path_dfs(chromosome, bam, gene,
                                                                                     unique_path[-1][1], read_set,
                                                                                     +1, '-', genome, min_coverage,
                                                                                     annotated_transcripts[
                                                                                         transcript].start_codon[
                                                                                         0].location.start.position)

                    for upstream_path, upstream_seq in zip(upstream_paths, upstream_seqs):
                        if upstream_path[-1][1] == annotated_transcripts[transcript].start_codon[
                                0].location.start.position:
                            full_upstream_paths.append(upstream_path)
                            full_upstream_seqs.append(upstream_seq)
                            supported_upstream_paths.append(upstream_path)
                            supported_upstream_seqs.append(upstream_seq)
                        else:
                            annotated_exon_path, annotated_seq = find_path_annotated(gene.reverse(),
                                                                               annotated_transcripts[
                                                                                   transcript].start_codon[
                                                                                   0].location.start.position,
                                                                               upstream_path[-1][1],
                                                                               annotated_transcripts[
                                                                                   transcript].splice_list, -1)
                            if len(annotated_exon_path):
                                full_upstream_paths.append(upstream_path + annotated_exon_path[0])
                                full_upstream_seqs.append(annotated_seq[0] + upstream_seq)
                                supported_upstream_paths.append(upstream_path)
                                supported_upstream_seqs.append(upstream_seq)

                    if not len(upstream_paths):
                        annotated_exon_path, annotated_seq = find_path_annotated(gene.reverse(), annotated_transcripts[
                            transcript].start_codon[0].location.start.position, unique_path[-1][1],
                                                                           annotated_transcripts[
                                                                               transcript].splice_list, -1)
                        if len(annotated_exon_path):
                            full_upstream_paths.append(annotated_exon_path[0])
                            full_upstream_seqs.append(annotated_seq[0])
                            supported_upstream_paths.append(None)
                            supported_upstream_seqs.append("")
                    if len(full_upstream_paths) or unique_path[-1][1] >= annotated_transcripts[
                                                                                   transcript].start_codon[
                                                                                   0].location.start.position:
                        full_downstream_paths, full_downstream_seqs, _ = find_path_dfs(chromosome, bam, gene.reverse(),
                                                                                       unique_path[0][0], read_set,
                                                                                       -1, '-', genome, min_coverage)
                        if len(full_upstream_paths) and len(full_downstream_paths):
                            for combination in itertools.product(range(len(full_upstream_paths)),
                                                                 range(len(full_downstream_paths))):
                                graph_path = full_downstream_paths[combination[1]] + tuple(
                                    [edge[:4] for edge in unique_path]) + full_upstream_paths[combination[0]]
                                graph_path = graph_path[::-1]
                                graph_seq = full_upstream_seqs[combination[0]] + reverse_complement(kmer_seq) + \
                                            full_downstream_seqs[combination[1]]
                                upstream_seq_temp = full_upstream_seqs[combination[0]] + reverse_complement(''.join(
                                        [edge[4] for edge in unique_path if
                                         edge[2] != "splice" and edge[2] != "del" and edge[1] > novel_splice[1]]))

                                if contain_stop_codon(upstream_seq_temp):
                                    continue

                                frame_before = -len(upstream_seq_temp) % 3

                                supported_seq, upstream_len, downstream_len = get_supported_seq(
                                    supported_upstream_seqs[combination[0]], full_downstream_seqs[combination[1]],
                                    unique_path, frame_before, novel_splice, '-', length)
                                supported_path = get_supported_path(upstream_len, downstream_len, novel_splice, '-',
                                                                    graph_path)
                                mut_peptide = translate(supported_seq)
                                count = output_peptide(mut_peptide, supported_seq, output_file, fasta_file, length,
                                                       upstream_len, generated_peptides, novel_splice, supported_path,
                                                       '-', chromosome, graph_path, graph_seq,
                                                       annotated_transcripts[transcript].genename[0], peptide_count)
                                
                                if count:
                                    peptide_count = count

                        elif len(full_upstream_paths) and not len(full_downstream_paths):
                            for full_upstream_path, full_upstream_seq, supported_upstream_seq in zip(
                                    full_upstream_paths, full_upstream_seqs, supported_upstream_seqs):
                                graph_path = tuple([edge[:4] for edge in unique_path]) + full_upstream_path
                                graph_path = graph_path[::-1]
                                graph_seq = full_upstream_seq + reverse_complement(kmer_seq)
                                upstream_seq_temp = full_upstream_seq + reverse_complement(''.join(
                                    [edge[4] for edge in unique_path if
                                     edge[2] != "splice" and edge[2] != "del" and edge[1] > novel_splice[1]]))

                                if contain_stop_codon(upstream_seq_temp):
                                    continue

                                frame_before = -len(upstream_seq_temp) % 3

                                supported_seq, upstream_len, downstream_len = get_supported_seq(supported_upstream_seq,
                                                                                                '', unique_path,
                                                                                                frame_before,
                                                                                                novel_splice, '-',
                                                                                                length)
                                supported_path = get_supported_path(upstream_len, downstream_len, novel_splice, '-',
                                                                    graph_path)
                                mut_peptide = translate(supported_seq)
                                count = output_peptide(mut_peptide, supported_seq, output_file, fasta_file, length,
                                                       upstream_len, generated_peptides, novel_splice, supported_path,
                                                       '-', chromosome, graph_path, graph_seq,
                                                       annotated_transcripts[transcript].genename[0], peptide_count)

                                if count:
                                    peptide_count = count

                        elif unique_path[-1][1] >= annotated_transcripts[transcript].start_codon[
                            0].location.start.position and len(full_downstream_paths):
                            for full_downstream_path, full_downstream_seq in zip(full_downstream_paths,
                                                                                 full_downstream_seqs):
                                graph_path = full_downstream_path + tuple(
                                    [edge[:4] for edge in unique_path if edge[1] <=
                                        annotated_transcripts[transcript].start_codon[0].location.start.position])
                                graph_path = graph_path[::-1]
                                graph_seq = reverse_complement(''.join([edge[4] for edge in unique_path if
                                                                        edge[2] != "splice" and edge[
                                                                            2] != "del" and edge[1] <=
                                                                        annotated_transcripts[transcript].start_codon[
                                                                    0].location.start.position])) + full_downstream_seq
                                upstream_seq_temp = reverse_complement(''.join([edge[4] for edge in unique_path if
                                                                                edge[2] != "splice" and edge[
                                                                                    2] != "del" and novel_splice[
                                                                                    1] < edge[1] <=
                                                                                annotated_transcripts[
                                                                                    transcript].start_codon[
                                                                                    0].location.start.position]))

                                if graph_path[0][1] != annotated_transcripts[transcript].start_codon[
                                        0].location.start.position or contain_stop_codon(upstream_seq_temp):
                                    continue

                                frame_before = -len(upstream_seq_temp) % 3

                                supported_seq, upstream_len, downstream_len = get_supported_seq('', full_downstream_seq,
                                                                                                unique_path,
                                                                                                frame_before,
                                                                                                novel_splice, '-',
                                                                                                length)

                                supported_path = get_supported_path(upstream_len, downstream_len, novel_splice, '-',
                                                                    graph_path)
                                mut_peptide = translate(supported_seq)
                                count = output_peptide(mut_peptide, supported_seq, output_file, fasta_file, length,
                                                       upstream_len, generated_peptides, novel_splice, supported_path,
                                                       '-', chromosome, graph_path, graph_seq,
                                                       annotated_transcripts[transcript].genename[0], peptide_count)

                                if count:
                                    peptide_count = count

                        elif unique_path[-1][1] >= annotated_transcripts[transcript].start_codon[
                            0].location.start.position and not len(full_upstream_paths) and not len(
                                full_downstream_paths):
                            graph_path = tuple([edge[:4] for edge in unique_path if edge[1] <=
                                        annotated_transcripts[transcript].start_codon[0].location.start.position])
                            graph_path = graph_path[::-1]
                            graph_seq = reverse_complement(''.join([edge[4] for edge in unique_path if
                                                                    edge[2] != "splice" and edge[
                                                                        2] != "del" and edge[1] <=
                                                                    annotated_transcripts[transcript].start_codon[
                                                                        0].location.start.position]))

                            upstream_seq_temp = reverse_complement(''.join([edge[4] for edge in unique_path if
                                                                            edge[2] != "splice" and edge[
                                                                                2] != "del" and novel_splice[
                                                                                1] < edge[1] <=
                                                                            annotated_transcripts[
                                                                                transcript].start_codon[
                                                                                0].location.start.position]))

                            if graph_path[0][1] != annotated_transcripts[transcript].start_codon[
                                        0].location.start.position or contain_stop_codon(upstream_seq_temp):
                                continue

                            frame_before = -len(upstream_seq_temp) % 3

                            supported_seq, upstream_len, downstream_len = get_supported_seq('', '', unique_path,
                                                                                            frame_before, novel_splice,
                                                                                            '-', length)
                            supported_path = get_supported_path(upstream_len, downstream_len, novel_splice, '-',
                                                                graph_path)
                            mut_peptide = translate(supported_seq)
                            count = output_peptide(mut_peptide, supported_seq, output_file, fasta_file, length,
                                                   upstream_len, generated_peptides, novel_splice, supported_path, '-',
                                                   chromosome, graph_path, graph_seq,
                                                   annotated_transcripts[transcript].genename[0], peptide_count)

                            if count:
                                peptide_count = count

    output_file.close()
    fasta_file.close()
    bam.close()
    kmer_dat.close()
    run_netMHCpan(sample, chromosome, length, HLA_string, HLA_II_string, neoantigen_path, netMHCpan_path,
                  netMHCIIpan_path)
    combine_table(sample, neoantigen_path, length, chromosome)
    logging.info("Done!")


if __name__ == '__main__':
    main()