import time
import sys

start = time.time()

import cPickle


def index(proteome_path, q=4):
    #
    # Proteins in the target_file must all occupy one line exactly!
    #
    proteins = {}
    peptides = {}

    acc = ""
    code = 0
    f = open(proteome_path + "proteome.fasta")
    for l in f:
        if l.startswith(">"):
            acc = l.strip()[1:]
        else:
            #The entire protein must be on a single line!!!
            code += 1
            seq = l.strip()
            proteins[code] = (acc, seq)
            for i in range(len(seq)-q+1):
                pep = seq[i:i+q].replace('L', 'I')
                if not pep in peptides:
                    peptides[pep] = set()
                peptides[pep].add(code)
    db_name = proteome_path + "proteome.pickle"
    f = open(db_name, 'wb')
    cPickle.dump(q, f)
    cPickle.dump(proteins, f)
    cPickle.dump(peptides, f)
    f.close()


proteome_dir = sys.argv[1]
if not proteome_dir.endswith("/") :
    proteome_dir += "/"
index(proteome_dir)

stop = time.time()

print >> sys.stderr, "Indexing of %s finished in %.3f seconds" % (proteome_dir,stop-start)
