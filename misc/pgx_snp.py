import re
import sys

q = 0
proteins = None
peptides = None

def lookup(aPeptide):
    global q
    global peptides
    global proteins

    target = aPeptide.replace("L", "I")
    output = []

    counter = 0
    scores = {}
    for i in range(len(target)-q+1):
        counter += 1
        if target[i:i+q] in peptides:
            candidates = peptides[target[i:i+q]]
            for protein in candidates:
                if not protein in scores:
                    scores[protein] = 0
                scores[protein] += 1
    meta_score = {}
    for i in range(1,len(target)-q+2):
        meta_score[i] = 0
    for code in scores:
        meta_score[scores[code]] += 1
    for i in range(1,len(target)-q+2):
        print >> sys.stderr, "%d: %d" % (i, meta_score[i])    
    for code in scores:
        if scores[code] >= counter - q:
            transeq = proteins[code][1].replace("L", "I")
            for offset in range(0,len(transeq)-len(target)+1):
                mismatch = 0
                for internal in range(len(target)):
                    if transeq[offset + internal] != target[internal]:
                        mismatch += 1
                if mismatch < 2:
                    output.append((proteins[code][0], offset+1, proteins[code][1][offset:(offset+len(target))]))
    return output


import time
start = time.time()
proteome = sys.argv[1]
if len(sys.argv) == 3:
    infile = fopen(sys.argv[2])
else:
    infile = sys.stdin

import cPickle

if not proteome.endswith("/"):
    proteome += "/"
f = open(proteome + 'proteome.pickle', 'rb')
q = cPickle.load(f)
proteins = cPickle.load(f)
peptides = cPickle.load(f)
f.close()

for l in infile:
    pep = l.strip().split()[0]
    matches = lookup(pep)
    for match in matches:
        print >> sys.stdout, "%s\t%s\t%d" % (pep, match[0], match[1])
# there is no harm in closing stdin... http://effbot.org/pyfaq/why-doesn-t-closing-sys-stdout-stdin-stderr-really-close-it.htm
infile.close()
stop = time.time()
print >> sys.stderr, "query processed in %.3f seconds" % (stop-start)
