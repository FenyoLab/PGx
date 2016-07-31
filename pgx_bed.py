import sys

mapping = {}

def load_mapping(bed_full_path):
    f = open(bed_full_path)
    for l in f:
        (chrom, start, stop, protein, dunno1, strand, dunno2, dunno3, dunno4, exons_count, exon_lengths, exon_starts) = l.strip().split()
        if exon_lengths.endswith(","):
            exon_lengths = exon_lengths[:-1]
        if exon_starts.endswith(","):
            exon_starts = exon_starts[:-1]
       
        lengths = map(int, exon_lengths.split(","))
        starts = map(int, exon_starts.split(","))
        mapping[protein] = (chrom, int(start), int(stop), strand, protein, int(exons_count), lengths, starts)
    f.close()


def process(peptide, mapping_entry, aa_start):
    (chrom, start, stop, strand, protein, exon_count, lengths, starts) = mapping_entry
    if strand == "+":
        start_target = (aa_start - 1) * 3
        tot_start = 0
        feature_start = -1
        #Otherwise IDE complains about curr_exon being potentially undefined...
        curr_exon = -1
        for curr_exon in range(exon_count):
            tot_start = starts[curr_exon]
            if lengths[curr_exon] > start_target:
                feature_start = tot_start + start_target
                break
            else:
                start_target -= lengths[curr_exon]
        if feature_start < 0:
            print >> sys.stderr, "couldn't place %s in %s starting at #AA %d" % (peptide, protein, aa_start)
            return ""
        end_target = len(peptide) * 3
        block_count = 0
        block_lengths = []
        block_starts = []

        overall_feature_start = start + feature_start
        current_start = overall_feature_start

        while end_target > (start + starts[curr_exon] + lengths[curr_exon] - current_start):
            block_count += 1
            block_starts.append(current_start-overall_feature_start)
            # There will be iterations where this seems like a futile cycle...
            block_lengths.append(start + starts[curr_exon] + lengths[curr_exon] - current_start)
            end_target -= (start + starts[curr_exon] + lengths[curr_exon] - current_start)
            curr_exon += 1
            if curr_exon == exon_count:
                print >> sys.stderr, "couldn't place %s in %s starting at #AA %d" % (peptide, protein, aa_start)
                return ""
            #This is the start of the next exon... i.e. starts[curr_exon] which looks redundant but wasn't
            #on the first iteration... I'm sure this whole loop structure can be simplified...
            current_start = start + starts[curr_exon]
        block_count += 1
        block_starts.append(current_start - overall_feature_start)
        block_lengths.append(end_target)
        feature_end = current_start+end_target-1
        return "\t".join([chrom, repr(overall_feature_start), repr(feature_end+1), peptide, '1000', strand, repr(overall_feature_start), repr(feature_end+1), '0,255,0', repr(block_count), ",".join(map(str, block_lengths)), ",".join(map(str, block_starts))])
    else:
        start_target = (aa_start-1) * 3
        tot_start = 0
        feature_start = -1
        #Otherwise IDE complains about curr_exon being potentially undefined...
        curr_exon = -1
        for curr_exon in range(exon_count-1, -1, -1):
            tot_start = starts[curr_exon] + lengths[curr_exon]
            if lengths[curr_exon] > start_target:
                feature_start = tot_start - start_target
                break
            else:
                start_target -= lengths[curr_exon]
        if feature_start < 0:
            print >> sys.stderr, "couldn't place %s in %s starting at #AA %d" % (peptide, protein, aa_start)
            return ""
        end_target = len(peptide) * 3
        block_count = 0
        block_lengths = []
        block_starts = []

        overall_feature_start = start + feature_start
        current_start = overall_feature_start

        while end_target > (current_start - start - starts[curr_exon]):
            block_count += 1
            block_starts = [starts[curr_exon]] + block_starts
            # There will be iterations where this seems like a futile cycle...
            block_lengths = [(current_start - start - starts[curr_exon])] + block_lengths
            end_target -= (current_start - start - starts[curr_exon])
            curr_exon -= 1
            if curr_exon == -1:
                print >> sys.stderr, "couldn't place %s in %s starting at #AA %d" % (peptide, protein, aa_start)
                return ""
                #This is the start of the next exon... i.e. starts[curr_exon] which looks redundant but wasn't
            #on the first iteration... I'm sure this whole loop structure can be simplified...
            current_start = start + starts[curr_exon] + lengths[curr_exon]
        block_count += 1
        feature_end = current_start-end_target
        block_starts = [0] + map(lambda x: x-(feature_end-start), block_starts)
        block_lengths = [end_target] + block_lengths
        return "\t".join([chrom, repr(feature_end), repr(overall_feature_start), peptide, '1000', strand, repr(feature_end), repr(overall_feature_start), '0,255,0', repr(block_count), ",".join(map(str, block_lengths)), ",".join(map(str, block_starts))])


def batch_process(f):
    print >> sys.stdout, 'track name=peptides description="Peptides identified by Mass Spectrometry" useScore=0 itemRgb="On"'
    for l in f:
        vals = l.strip().split("\t")
        peptide = vals[0]
        protein = vals[1]
        aa_start = int(vals[2])
        lines = set()
        if not protein in mapping:
            print >> sys.stderr, "Hit in a protein (%s) which is not mapped by the reference bed file of the '%s' proteome." % (protein, proteome[:-1])
            continue
        else:
            retval = process(peptide, mapping[protein], aa_start)
            if (retval == "") or (retval in lines):
                continue
            else:
                lines.add(retval)
                print >> sys.stdout, retval

import time
start_time = time.time()

proteome = sys.argv[1]
if len(sys.argv) == 3:
    infile = open(sys.argv[2])
else:
    infile = sys.stdin

if not proteome.endswith("/"):
    proteome += "/"
load_mapping(proteome + "proteome.bed")
batch_process(infile)
# there is no harm in closing stdin... http://effbot.org/pyfaq/why-doesn-t-closing-sys-stdout-stdin-stderr-really-close-it.htm
infile.close()
end_time = time.time()
print >> sys.stderr, "peptides mapped in %.3f seconds" % (end_time-start_time)
