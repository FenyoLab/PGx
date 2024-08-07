import re
import sys


def lookup(aPeptide):
    global q
    global peptides
    global proteins

    target = aPeptide.replace("L", "I")
    output = []

    if target[0:q] in peptides:
        candidates = peptides[target[0:q]]
        for i in range(1, len(target)-q+1):
            if target[i:i+q] in peptides:
                candidates = candidates.intersection(peptides[target[i:i+q]])
            else:
                candidates = set()
                break
    else:
        candidates = set()

    for code in candidates:
        #Obviously, the proteins should be pre-I/L transformed...
        transeq = proteins[code][1].replace("L", "I")
        for m in re.finditer('(?=%s)' % target, transeq):
            output.append((proteins[code][0], m.start()+1))
    return output


if __name__ == "__main__":

    q = 0
    proteins = None
    peptides = None

    import time
    start = time.time()
    proteome = sys.argv[1]
    if len(sys.argv) == 3:
        infile = open(sys.argv[2])
    else:
        infile = sys.stdin

    import pickle

    if not proteome.endswith("/"):
        proteome += "/"
    f = open(proteome + 'proteome.pickle', 'rb')
    q = pickle.load(f)
    proteins = pickle.load(f)
    peptides = pickle.load(f)
    f.close()

    for l in infile:
        pep = l.strip().split()[0]
        matches = lookup(pep)
        for match in matches:
            print("%s\t%s\t%d" % (pep, match[0], match[1]), file=sys.stdout)

    # there is no harm in closing stdin... http://effbot.org/pyfaq/why-doesn-t-closing-sys-stdout-stdin-stderr-really-close-it.htm
    infile.close()
    stop = time.time()
    print("query processed in %.3f seconds" % (stop-start), file=sys.stderr)
