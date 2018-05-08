from subprocess import call
import numpy as np

def extract(sequence_string):
    f = open("dot.ps", 'r')
    
    record = False
    
    start = "%start of base pair probability data\n"
    end = "0.9500000 lbox"
    
    bpp = []

    for line in f:
        if end in line:
            record = False
        if record:
            BPP = [int(line.split()[0]), int(line.split()[1]), float(line.split()[2])]
            bpp.append(BPP)
        if start == line:
            record = True

    np.savetxt('vienna_bpp_data/{}.txt'.format(sequence_string), np.array(bpp))
    f.close()

def main():
    seq = open("experimental/R101_Sequence.txt", 'r')
    
    for ii in range(130):#12089):
        SEQUENCE = seq.readline().rstrip()
        sequence = open("sequence_input.txt", 'w')
        sequence.write(SEQUENCE)
        sequence.close()
        
        call(["RNAfold", "sequence_input.txt", "-p", "--noGU"])
        
        extract(SEQUENCE)
    
    seq.close()

main()
