import pysam
import re
import argparse
import logging


class Trie():
    """Regex::Trie in Python. Creates a Trie out of a list of words. The trie can be exported to a Regex pattern.
    The corresponding Regex should match much faster than a simple Regex union."""

    def __init__(self):
        self.data = {}

    def add(self, word):
        ref = self.data
        for char in word:
            ref[char] = char in ref and ref[char] or {}
            ref = ref[char]
        ref[''] = 1

    def dump(self):
        return self.data

    def quote(self, char):
        return re.escape(char)

    def _pattern(self, pData):
        data = pData
        if "" in data and len(data.keys()) == 1:
            return None

        alt = []
        cc = []
        q = 0
        for char in sorted(data.keys()):
            if isinstance(data[char], dict):
                try:
                    recurse = self._pattern(data[char])
                    alt.append(self.quote(char) + recurse)
                except:
                    cc.append(self.quote(char))
            else:
                q = 1
        cconly = not len(alt) > 0

        if len(cc) > 0:
            if len(cc) == 1:
                alt.append(cc[0])
            else:
                alt.append('[' + ''.join(cc) + ']')

        if len(alt) == 1:
            result = alt[0]
        else:
            result = "(?:" + "|".join(alt) + ")"

        if q:
            if cconly:
                result += "?"
            else:
                result = "(?:%s)?" % result
        return result

    def pattern(self):
        return self._pattern(self.dump())


logging.basicConfig(level=logging.DEBUG, format="%(asctime)-15s [%(processName)s.%(levelname)s] %(message)s")
parser = argparse.ArgumentParser(description="Utility for retrieving tumor specific kmers in RNA-seq reads.")
parser.add_argument("Kmer_file", type=str, nargs='?', help="provide Kmer file here")
parser.add_argument("input_bam_file", type=str, nargs='?', help="provide input bam file path here")
parser.add_argument("out_bam_file", type=str, nargs='?', help="provide output bam file path here")
args = parser.parse_args()

cigar_map = {0:'M', 1:'I', 2:'D', 3:'N', 4:'S', 5:'H', 6:'P', 7:'=', 8:'X', 9:'B'}

samfile = pysam.AlignmentFile(args.input_bam_file, "rb")
kmer_reads = pysam.AlignmentFile(args.out_bam_file, "wb", template=samfile)

trie = Trie()
with open(args.Kmer_file) as f:
    for line in f:
        line_split = line.strip().split('\t')
        trie.add(line_split[0])
print trie.pattern()
pattern = re.compile(trie.pattern())
#pattern = re.compile("\*")

#i = 0
for read in samfile.fetch():
#    i+=1
#    if i > 100000:
#        break
    if not read.is_unmapped and not read.is_secondary and read.is_proper_pair and read.is_paired and\
            not read.is_duplicate and not read.is_supplementary:

        found = pattern.search(read.query_sequence)
        if not found:
            continue
        print read.query_name
        start_index = found.start()
        end_index = found.end()
        quality_sring = read.to_string().split('\t')[10][start_index:end_index]

        print "start " + str(start_index)
        print "end " + str(end_index)
        #print len(read.query_qualities)
        #print read.query_qualities[start_index:end_index]

        print "original sequence " + str(read.query_sequence)
        print "original quality " + str(read.to_string().split('\t')[10])
        print read.to_string().split('\t')[10][start_index:end_index]
        read.query_sequence = read.query_sequence[start_index:end_index]
        print "new seq " + read.query_sequence


        current_ind = 0
        in_kmer = False
        cigarString_temp = ""
        Kmer_len = end_index - start_index
        print read.get_reference_positions(full_length=True)
        print read.get_reference_positions(full_length=True)[start_index:end_index]



        if all(ref_pos is None for ref_pos in read.get_reference_positions(full_length=True)[start_index:end_index]):
            continue

        for ref_pos in read.get_reference_positions(full_length=True)[start_index:end_index]:
            if ref_pos != None:
                read.reference_start = ref_pos
                break

        for operation, count in read.cigartuples:
            print operation, count
            print current_ind
            print in_kmer
            if current_ind >= end_index:
                break

            if cigar_map[operation] == 'N' or cigar_map[operation] == 'D':
                if in_kmer:
                    cigarString_temp += str(count) + cigar_map[operation]

            elif cigar_map[operation] == 'M' or cigar_map[operation] == 'I' or cigar_map[operation] == 'S':
                if current_ind + count > start_index:
                    cigarString_temp += str(min(end_index, current_ind + count) - max(current_ind, start_index)) + \
                                        cigar_map[operation]
                    in_kmer = True
                current_ind += count
            else:
                logging.warn('Unexpected cigar op {}'.format(cigar_map[operation]))


        print "original cigar" + read.cigarstring
        print "new cigar " + cigarString_temp
        read.cigarstring = cigarString_temp
        read.query_qualities = pysam.qualitystring_to_array(quality_sring)
        print "new qual " + quality_sring
        print "new qual " + str(read.query_qualities)
        print read.to_string()
        kmer_reads.write(read)


kmer_reads.close()
samfile.close()



