from __future__ import division

def read_windows(fin, ws):
    while True:
        s = fin.read(ws)
        if not s: break

        window = s.strip()
        window = window.replace("\n", "")

        yield window

def calc_gc_content(window):
    c_count = window.count('C')
    g_count = window.count('G')
    L = len(window)

    return((c_count+g_count) / L)

def find_gc_content(fin):
    gc_content = []
    fin.readline()
    for window in read_windows(fin, 100):
        gc = calc_gc_content(window)
        gc_content.append(gc)

    return(gc_content)

def write_to_file(gc_content, fout_name):
    fout= open(fout_name, "w")

    for gc in gc_content:
        fout.write(str(gc) + "\n")

def main():
    fin_name = "data/Saccharomyces_cerevisiae.R64-1-1.75.dna_rm.chromosome.III.fa"
    fout_name = "data/cg.dat"

    with open(fin_name) as fin:
        gc_content = find_gc_content(fin)
        write_to_file(gc_content, fout_name)

    return

if __name__ == "__main__":
    main()
