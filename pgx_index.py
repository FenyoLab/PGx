import time
import sys

start = time.time()

import cPickle

def index(proteome_path, q=4):
    proteins = {}
    peptides = {}
    code_container = [0] # This only exists because closures are not possible in Python2.7 (no "nonlocal" keyword)

    def process_entry(acc,seq):
        #
        # This closure exists only to deal with the fact that there is no mandated endfile delimiter in FASTA files
        # (and I refused to write the same code twice in the line processing for-loop). If anybody has a more elegant
        # way to do this, I'm all ears! (deposit in pull-request/suggestion-box)
        #
        if seq == "" or acc == "": # This is both useful as a sanity check and eliminates if statements for the external for loop.
            return
        #
        # Python 2.7 closures are broken: no "nonlocal" keyword
        # see: https://www.python.org/dev/peps/pep-3104/
        code = code_container[0] # 
        code += 1                # 
        code_container[0] = code # 
        proteins[code] = (acc, seq)
        for i in range(len(seq)-q+1):
            pep = seq[i:i+q].replace('L', 'I')
            if not pep in peptides:
                peptides[pep] = set()
            peptides[pep].add(code)
    
    acc = ""
    seq = ""
    f = open(proteome_path + "proteome.fasta")
    for l in f:
        if l.startswith(">"):
            process_entry(acc,seq)
            seq = ""
            acc = l.strip()[1:]
        else:
            nu_seq = l.strip()
            seq += nu_seq
    process_entry(acc,seq)

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
