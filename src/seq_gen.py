"""
"""

import RNA
import numpy as np
import os
import itertools
import argparse
import multiprocess as mp
from folding_wrapper import *
import matplotlib.pyplot as plt



def fold_with_rafft(sequence) :
    cmd = "python RAFFT/rafft.py -s {}"
    p = os.popen(cmd.format(sequence))
    rst = p.read().split()
    if len(rst) > 0 :
        return rst[1]
    else :
        print("Error during the folding")
        return None

def pprafft_fold (sequences) :
    pool = mp.Pool(mp.cpu_count())
    result = numpy.array(pool.map(fold_with_rafft, sequences))
    pool.close()
    return result


def genNeighboors(ref_seq, dist, nb_seq) :
    data = []
    list_ = []
    alpha = ["A","U","C","G"]
    for tpl in itertools.combinations(range(len(ref_seq)),dist):
        s = np.array(list(ref_seq))
        for i in range(dist) :
            set_ =list(set(alpha)-set([s[tpl[i]]]))
            s[i] = set_[np.random.randint(0,3)]
        data.append((ref_seq,''.join(s)))
        list_ += [''.join(s)]
        if len(data)> nb_seq :
            return data, list_
    return data, list_

def neutrality(listOfStructures) :
    return len(set(listOfStructures))/(len(listOfStructures)*1.)

def parse_arguments():
    """Parsing command line
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--length', '-l', help="sequence length")
    parser.add_argument('--sequence_number', '-sn', help="sequence number to generate for each mutation distance ")
    parser.add_argument('--max_dis', '-md', help="maximum mutation allow")

    return parser.parse_args()



def main() :
    args = parse_arguments()
    ref_seq = RNA.random_string(int(args.length),"AUCG")
    dt, seqs = genNeighboors(ref_seq,int(args.max_dis), int(args.sequence_number))

    print("Reference sequence = {}".format(ref_seq))

    nt_data = []
    nt_data_rafft = []
    for d in range(int(args.max_dis)) :
        folded = ppRNAfold(seqs)[0]
        nt_data.append(neutrality(folded))

        folded = pprafft_fold(seqs)
        nt_data_rafft.append(neutrality(folded))
        _,seqs = genNeighboors(ref_seq,d, int(args.sequence_number))

        print("mutant at distance {}".format(d))

    plt.plot(range(int(args.max_dis)), nt_data,'o-', "RNAfold")
    plt.plot(range(int(args.max_dis)), nt_data_rafft,'o-',label='RAFFT')
    plt.xlabel("Steps (Mutations)")
    plt.ylabel("Neutrality")
    plt.savefig("../images/neutrality.pdf")
    plt.show()




if __name__ == "__main__" :
    main()
