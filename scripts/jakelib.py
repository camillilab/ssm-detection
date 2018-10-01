

from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
import itertools
from collections import defaultdict
import re
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
import subprocess
import os


# Static sequence file instance
class SeqFile:

    # initialize the instance; load the sequences using BioPython SeqIO
    def __init__(self, seq_file, file_type='fasta', use_threads=1):
        with open(seq_file, 'rU') as handle:
            self.records = []
            self.names = []
            for record in SeqIO.parse(handle, file_type):
                self.records.append(record)
                self.names.append(record.name)

        self.sequence_file = seq_file
        self.type = file_type
        self.threads = use_threads

    def seq_size(self):
        i = 0
        for record in self.records:
            i += len(record)
        return i


# AlignmentFile instance for the read file. Focuses on BAM. Uses samtools and other execs for speed.
class AlignmentFile:

    # initialize the instance
    def __init__(self, alignment_file, chr_name='', sort=True, index=True, use_threads=1):
        self.alignment_file = alignment_file

        # get the number of reads and the header using samtools
        #process = subprocess.run(['samtools', 'view', '-c', alignment_file], stdout=subprocess.PIPE, encoding='utf-8')
        #self.num_reads = int(process.stdout)

        process = subprocess.run(['samtools', 'view', '-H', alignment_file],
                                 stdout=subprocess.PIPE, encoding='utf-8')
        self.header = process.stdout

        self.sorted_alignment_file = chr_name + '_sorted.bam'

        self.region_name = chr_name

        self.threads = str(use_threads)
        # sort the bam file
        # make these methods
        if sort:
            #print("Sorting...")
            # will stderr to devnull suppress that annoying bam_sort_core message?
            subprocess.run(['samtools', 'sort', self.alignment_file, '-o', self.sorted_alignment_file,
                            '-@', self.threads], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if index:
            #print("Indexing...")
            subprocess.run(['samtools', 'index', self.sorted_alignment_file])

    # an ugly workaround to keeping a region name necessary for sorting later on
    def give_name(self, name):
        self.region_name = name

    # dumps reads to disk based on position and name
    # fix start and end so that if you don't specify either, it'll grab the whole file
    def dump(self, outfile, start=-1, end=-1, name=''):

        # if end is default, set to max length of sequence
        if end == -1:
            subprocess.run(['samtools', 'view', '-h',
                            self.sorted_alignment_file, name+':'+str(start), '-@', self.threads],
                           stdout=open(outfile, 'w'),
                           encoding='utf-8', universal_newlines=True, bufsize=-1)
        else:
            subprocess.run(['samtools', 'view', '-h',
                            self.sorted_alignment_file, name+':'+str(start)+'-'+str(end), '-@', self.threads],
                           stdout=open(outfile, 'w'), encoding='utf-8', universal_newlines=True, bufsize=-1)
        return

    # retrieves reads based on position and name
    def fetch(self, start=-1, end=-1, name=''):

        # if end is default, set to max length of sequence
        if end == -1:
            process = subprocess.run(['samtools', 'view',
                                      self.sorted_alignment_file, name+':'+str(start), '-@', self.threads],
                                     stdout=subprocess.PIPE,
                                     encoding='utf-8', universal_newlines=True, bufsize=-1)
        else:
            process = subprocess.run(['samtools', 'view',
                                      self.sorted_alignment_file, name+':'+str(start)+'-'+str(end), '-@', self.threads],
                                     stdout=subprocess.PIPE, encoding='utf-8', universal_newlines=True, bufsize=-1)
        return process.stdout.split('\n')

    # filters reads based on CIGAR; returns POS, CIGAR
    def filter(self, start=-1, end=-1, name='', cigar='.*'):
        # get positional reads
        reads = self.fetch(start, end, name)
        # compile the regex
        read_regex = re.compile(cigar)
        # yield reads meeting our cigar criteria
        for read in reads:
            try:
                read_cigar = read.split('\t')[5]  # CIGAR value is the fifth column

                # if it matches, yield the read
                if read_regex.match(read_cigar):
                    yield read
            except IndexError:
                pass

    # gets average depth of region
    # may need to change so that is excludes the read region itself...
    def average_depth(self, start=-1, end=-1, name=''):

        process = subprocess.run(['samtools', 'depth', '-r', name+':'+str(start)+'-'+str(end),
                                  self.sorted_alignment_file], stdout=subprocess.PIPE,
                                 encoding='utf-8', universal_newlines=True, bufsize=-1)
        depth = 0
        data = process.stdout.split('\n')
        for line in data:
            try:
                line_depth = int(line.split('\t')[2])
                depth += line_depth
            except IndexError:
                pass
        return depth / len(data)

    # gets read length of a bam file using samtools stats and text processing
    def get_read_length(self):
        process = subprocess.run('samtools stats ' + self.alignment_file + " | grep ^RL | cut -f 2- | awk '{print $1}'",
                                 shell=True, stdout=subprocess.PIPE, encoding='utf-8',
                                 universal_newlines=True, bufsize=-1)
        return process.stdout.split('\n')[0]


# Isolated sequence tools
# Returns a list containing all unique nucleotide combinations of a given range
def generate_repeat_seeds(min_size, max_size, procs=1):
    seed_range = list(range(min_size, max_size + 1))
    all_seeds = []
    with ProcessPoolExecutor(max_workers=procs) as executor:
        results = executor.map(nt_patterns, seed_range)
        for result in results:
            for item in result:
                all_seeds.append(item)
    return all_seeds


# generate all unique, non-repeated combinations of nucleotide patterns of a given length; returns a generator
def nt_patterns(length):
    combos = itertools.product('ACGT', repeat=length)
    all_combos = []
    for combo in combos:
        if is_clean(''.join(combo)):
            all_combos.append(''.join(combo))
    return all_combos


# function evaluates a generated pattern to make sure it is unique and isn't a repeat itself (ie. ATGATG)
# should help improve clarity of results and runtime speed
def is_clean(pattern):
    # if the length of the pattern is just one, just return okay
    if len(pattern) == 1:
        return True
    # if the pattern is a homopolymer, delete
    if pattern == pattern[0] * len(pattern):
        # print("Pattern {0} is a homopolymer; deleting...".format(pattern))
        return False
    # if the pattern is an odd length, shouldn't be a problem
    if len(pattern) % 2 != 0:
        return True
    # else, try to find patterns within this pattern up to half the length
    for i in range(2, int(len(pattern)/2)+1):
        # generate permutations
        pattern_seeds = itertools.permutations('AGCT', r=i)
        # look for matches
        for seed in pattern_seeds:
            match_positions = detect_pattern(pattern, ''.join(seed))
            # enumerate the generator
            matches = []
            for j in match_positions:
                matches.append(j)
            if len(matches) == (len(pattern) / len(''.join(seed))):
                return False
    return True


# detects non-overlapping patterns in a sequence
# returns a generator
def detect_pattern(seq, pattern):
    str_seq = str(seq)
    i = str_seq.find(pattern)
    while i != -1:  # find returns -1 when there are no more matches remaining
        yield i+1  # add 1 to index to align with actual sequence, position, which starts at one.
        i = str_seq.find(pattern, i+len(pattern))


# scans for repeat sequences
def scan_for_repeats(seq, pattern, num_repeat_threshold=5):

    # pull up matches for the minimum threshold
    # TODO Alter behavior to allow for an absolute length rather than a multiple of a pattern seed?
    match_positions = detect_pattern(seq, pattern=pattern*num_repeat_threshold)
    last_pos = -1
    for pos in match_positions:
        # check to make sure that this pos isn't contiguous with the previous match
        if pos == last_pos + len(pattern*num_repeat_threshold):
            pass
        else:
            # get the corrected index at the end of our detected repeat
            base = pattern*num_repeat_threshold
            i = (pos - 1) + len(base)
            # check ahead the length of the repeat to see if another one exists
            while pattern == seq[i:i + len(pattern)]:
                base += pattern
                i += len(pattern)
            # yield the pos and the extended repeat
            yield pos, base
        # assign this pos as the last pos recorded
        last_pos = pos


# writes data to disk
def write_to_disk(data, outfile):
    with open(outfile, 'w') as handle:
        for line in data:
            handle.write(line)
    return


# breaks down a cigar string into it's components (ie. 88M2I2D10M >- [88M, 2I, 2D, 10M] using regex
def parse_cigar(cigar):

    # define regex compiler
    cig_reg = re.compile('[0-9]*[A-Z]')

    # return the matches
    return cig_reg.findall(cigar)


def find_id_read_areas(cigar_parts):
    pos = 0
    id_areas = []
    for part in cigar_parts:

        # get the code and value
        code = part[-1]
        value = int(part[:-1])

        if code == 'M':
            # match, add to counter
            pos += value
        if code == 'N':
            # unknown nt, move the counter
            pos += value
        if code == 'I':
            # insertion detected, add a tuple of (pos, code, length), and move the counter
            area = (pos, code, value)
            id_areas.append(area)
            pos += value
        if code == 'D':
            # deletion detected, add a tuple of (pos, code, length), but keep the counter here
            area = (pos, code, value)
            id_areas.append(area)
    return id_areas


# detect ssm read patterns
def find_ssm(repeat_dict, aln, refseq, procs=1):
    ssm_data = dict()
    with ProcessPoolExecutor(max_workers=procs) as executor:
        futures = {executor.submit(_find_ssm_worker, k, repeat_dict[k][0], repeat_dict[k][1], aln, refseq):
                   k for k in repeat_dict}
        for future in as_completed(futures):
            try:
                data = future.result()
                ssm_data.update(data)
            except Exception as exc:
                print('generated an exception: {0}'.format(exc))
    return ssm_data


# worker function for find_ssm
def _find_ssm_worker(pos, repeat, repeat_size, aln, refseq):

    # data container
    ssm_data = defaultdict(float)

    # fetch reads containing insertion or deletion patterns matching the length of the repeat at the area of interest
    pattern_cigar = '.*' + str(len(repeat)) + '[ID]'
    # TODO eventually auto-detect the start/end based on read size?
    reads = aln.filter(start=pos-500, end=pos+500, name=aln.region_name, cigar=pattern_cigar)

    # see if these reads suggest a site-specific SSM event
    for read in reads:

        # split the read to get a list of parameters in the read
        params = read.split('\t')

        # parse the cigar string
        read_cigar = params[5]
        read_cigar_parts = parse_cigar(read_cigar)

        # find the positions within the read that suggest I/D of our pattern length
        # TODO We are missing any reads that happen to have a correct nt in the SSM region
        id_areas = find_id_read_areas(read_cigar_parts)
        repeat_id_areas = [k for k in id_areas if k[2] == len(repeat)]

        # get the read and the reference position in which the read starts
        read_seq = params[9]
        read_pos = int(params[3])

        for id_area in repeat_id_areas:

            # does the insertion or deletion map where we would expect a SSM to be?
            id_pos, id_code, id_length = int(id_area[0]), id_area[1], int(id_area[2])
            abs_id_area = id_pos + read_pos

            # the expected site of an indel are the borders of the indel
            # we know where our pattern begins and the length of the pattern, so calculate accordingly
            expected_pos_left = int(pos)
            expected_pos_right = int(pos) + len(repeat)*repeat_size

            if abs_id_area == expected_pos_left or abs_id_area == expected_pos_right:
                # we have a positional match!

                # now see if what the read suggests is actually the repeat
                # add an average agreement (min=0, max=1) to the dict
                read_pattern_seq = read_seq[id_pos:id_pos+id_length]
                j = 0
                for i in range(0, len(repeat)):
                    if read_pattern_seq[i] == repeat[i]:
                        j += 1

                ssm_data[(aln.region_name, pos, repeat, repeat_size, id_code)] += (j / len(repeat))

                # DEBUG AREA FOR READ DETECTION
                if (repeat == 'AGTC') and (repeat_size == 32):
                    print(read_cigar+'\t'+read)

    return ssm_data


# realized that I would need to save several query shifts...
def _blast_worker(blast_program, database, record):

    print('BLASTing {0}...'.format(record[0]))
    seq = record[1]

    # result_handle = NCBIWWW.qblast(blast_program, database, seq, megablast=True)
    # blast_results = NCBIXML.parse(result_handle)

    # return the first alignment
    # result = next(blast_results)

    # return result.alignments[0]
    # write seq to fasta
    new_record = SeqRecord(seq, id=record[0])
    outfile = '_'+record[0]+'.fa'
    with open(outfile, 'w') as outhandle:
        SeqIO.write(new_record, outhandle, "fasta")

    # now query, we get an XML file, no stdout
    record_xml_file = '_'+record[0]+'.xml'
    subprocess.run([blast_program, '-query', outfile, '-out', record_xml_file,
                    '-db', database, '-task', 'megablast', '-outfmt', '5'])

    # parse the xml file
    handle = open(record_xml_file)
    blast_results = NCBIXML.parse(handle)
    # return the first alignment
    result = next(blast_results)
    return result.alignments[0]


# BLAST a sequence over the interwebs - determine the most likely parent sequence from a series of contigs
# definitely use max threads, BLASTing this way somehow takes forever?
def detect_reference_sequence_blast(blast_program, database, records, num_procs=8):

    record_dict = defaultdict(int)
    title_dict = dict()
    query_shift_dict = dict()
    alignments = list()

    with ProcessPoolExecutor(max_workers=num_procs) as executor:
        futures = {executor.submit(_blast_worker, blast_program, database, k): k for k in records}
    for future in as_completed(futures):
        try:
            alignment = future.result()
            alignments.append(alignment)
        except Exception as exc:
            print("Generated an exception: {0}".format(exc))

    for i in range(len(alignments)):

        node_name = records[i][0]
        alignment = alignments[i]
        accession_num = alignment.accession
        title = alignment.hit_def
        query_start = int(alignment.hsps[0].query_start)
        subject_start = int(alignment.hsps[0].sbjct_start)

        record_dict[accession_num] += 1
        query_shift_dict[(accession_num, node_name)] = subject_start - query_start
        title_dict[accession_num] = title

    print("Accession numbers detected:")
    for key in record_dict:
        print('{0}: {1} times'.format(key, record_dict[key]))

    # what is the highest value?
    max_value = 0
    max_key = ''
    for key in record_dict:
        if record_dict[key] > max_value:
            max_key = key
            max_value = record_dict[key]

    # print the highest key
    print("Most likely reference sequence: {1} ({0})".format(title_dict[max_key], max_key))

    # only return the dict as { chr_name : shift_val }
    return_dict = {key[1]: query_shift_dict[key] for key in query_shift_dict if key[0] == max_key}
    #print(query_shift_dict)
    #print(return_dict)

    # delete all xml and fa files generated
    for record in records:
        os.remove('_'+record[0]+'.'+'fa')
        os.remove('_'+record[0]+'.xml')

    return max_key, return_dict

