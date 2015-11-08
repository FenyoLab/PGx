import time

start = time.time()

import urllib
import gzip

import subprocess

refseq_dir = "refseq_human/"


def download_refseq_human_proteome():
    fasta_url = "ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.protein.faa.gz"
    features_url = "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.protein.gpff.gz"
    seq_ali_url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refSeqAli.txt.gz"
    if not os.path.exists(refseq_dir):
        os.makedirs(refseq_dir)

    urllib.urlretrieve(fasta_url, refseq_dir+"proteome.fasta.gz")
    fasta_out = open(refseq_dir + "proteome.fasta", 'w')
    fasta_in = gzip.open(refseq_dir+"proteome.fasta.gz", 'rb')
    newline = False
    for l in fasta_in:
        if l.startswith(">"):
            if newline:
                fasta_out.write("\n")
            else:
                newline = True
            fasta_out.write(l.strip())
            fasta_out.write("\n")
        else:
            fasta_out.write(l.strip())
    fasta_in.close()
    fasta_out.close()

    urllib.urlretrieve(features_url, refseq_dir + "proteome.gpff.gz")
    pgx_util.gunzip(proteome_dir + "proteome.gpff.gz")
    urllib.urlretrieve(seq_ali_url, refseq_dir + "refSeqAli.txt.gz")
    pgx_util.gunzip(proteome_dir + "refSeqAli.txt.gz")
    p = subprocess.Popen(["perl", "refseq_gpff_to_bed.pl"], cwd=refseq_dir)
    p.wait()


download_refseq_human_proteome()
stop = time.time()

import datetime
today = datetime.date.today()

f = open(refseq_dir + "download_notes.txt", 'w')
print >> f, today.strftime("human proteome downloaded on the %d, %h %Y")
f.close()

# Reformated to ensure that there is one protein per row!
print "Downloaded human proteome, gpff, SeqAli files and generated a .bed file in %.2f seconds" % (stop-start)
