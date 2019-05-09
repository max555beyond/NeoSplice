import pysam
import regex as re
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


def main():
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

    pattern = re.compile(trie.pattern())
    print trie.pattern()
    for read in samfile.fetch():
        if not read.is_unmapped and not read.is_secondary and read.is_proper_pair and read.is_paired and\
                not read.is_duplicate and not read.is_supplementary and 'N' in read.cigarstring:

            matches = list(pattern.finditer(read.query_sequence, overlapped=True))

            if not matches:
                continue
            if len(matches) >= 2 :
                print 'happend'
            print matches
            print '\t'.join(read.to_string().split('\t'))

            for found in matches:

                start_index = found.start()
                end_index = found.end()
                print "start " + str(start_index)
                print "end " + str(end_index)
                if all(ref_pos is None for ref_pos in
                       read.get_reference_positions(full_length=True)[start_index:end_index]):
                    continue

                quality_string = read.to_string().split('\t')[10][start_index:end_index]
                print "qual " + quality_string

                kmer_read = pysam.AlignedSegment()
                print dir(kmer_read)
                kmer_read.query_name = read.query_name
                kmer_read.query_sequence = read.query_sequence[start_index:end_index]
                kmer_read.flag = read.flag
                kmer_read.reference_id = read.reference_id
                print kmer_read.query_name
                print kmer_read.query_sequence
                print kmer_read.flag
                print read.reference_id
                print kmer_read.reference_id
                for ref_pos in read.get_reference_positions(full_length=True)[start_index:end_index]:
                    if ref_pos is not None:
                        kmer_read.reference_start = ref_pos
                        break

                print kmer_read.reference_start
                kmer_read.mapping_quality = read.mapping_quality
                print kmer_read.mapping_quality

                current_ind = 0
                in_kmer = False
                cigarString_temp = ""

                for operation, count in read.cigartuples:
                    if current_ind >= end_index:
                        break

                    if cigar_map[operation] == 'N' or cigar_map[operation] == 'D':
                        if in_kmer:
                            cigarString_temp += str(count) + cigar_map[operation]

                    elif cigar_map[operation] == 'M' or cigar_map[operation] == 'I' or cigar_map[operation] == 'S':
                        if current_ind + count > start_index:
                            cigarString_temp += str(
                                min(end_index, current_ind + count) - max(current_ind, start_index)) + cigar_map[
                                                    operation]
                            in_kmer = True
                        current_ind += count
                    else:
                        logging.warn('Unexpected cigar op {}'.format(cigar_map[operation]))

                kmer_read.cigarstring = cigarString_temp
                print kmer_read.cigarstring
                kmer_read.query_qualities = pysam.qualitystring_to_array(quality_string)
                print kmer_read.query_qualities

                kmer_reads.write(kmer_read)
                print kmer_read

    kmer_reads.close()
    samfile.close()
    logging.info("Done!")


if __name__ == '__main__':
    main()