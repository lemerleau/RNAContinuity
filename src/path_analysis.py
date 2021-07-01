"""
  @author: Nono Saha Cyrille Merleau
  @email: nonosaha@mis.mpg.de
"""
#importing necessary libraries
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import RNA
import argparse as agp
import os
import json




def parse_arguments():

    parser = agp.ArgumentParser(description=__doc__, formatter_class=agp.RawTextHelpFormatter)
    parser.add_argument('--log_folder', '-lf', help="Folder name where the populations were stored. The folder should be in ../data/ken_rafft/<folder_name>")
    return parser.parse_args()


def main() :
    print("Hello welcome to my")
    args = parse_arguments()
    root_path = '../data/ken_rafft/'+str(args.log_folder)
    linage = pd.read_csv(root_path+"/evol_path.csv").values
    target = linage[0][4]
    print(linage)
    gen = 13057
    data = []
    list_dt = [gen]
    distinct_strcs = []
    sequences_path = []
    lst_seq = [linage[0][3]]
    for i in range(len(linage)-1):
        line = linage[i]
        if RNA.hamming_distance(line[4], linage[i+1][4]) == 0 :
            lst_seq.append(linage[i+1][3])
            list_dt.append(gen-i-1)
        else :
            print(line[4], RNA.hamming_distance(line[4], target))
            distinct_strcs.append(line[4])
            data.append(list_dt)
            list_dt = [gen-i-1]
            sequences_path.append([(line[3], line[4]),(linage[i+1][3], linage[i+1][4])])
            lst_seq = [linage[i+1][3]]

        #print(i)
    assert len(sequences_path) == len(data)
    Y_axis = []
    Y = []
    for i in range(len(data)):
        Y_axis +=[60-1.5*i]*len(data[::-1][i])
        Y +=[60-1.5*i]
    print(data, len(data))
    if np.array_equal(np.array(range(1,gen+1))[::-1],np.concatenate(data)) :
        print("Good data")
    else :
        print("Bad Data")

    print("Number of distinct structures: ", len(distinct_strcs), len(set(distinct_strcs)))

    ds = list(set(linage[:,4]))
    time_ds = {}
    labels = {}

    for i in range(len(ds)) :
        time_ds[ds[i]]=[]
        labels[ds[i]] = r'$S_{'+str(i)+"}$"

    path_data = {}
    for i in range(len(distinct_strcs)) :
        path_data[labels[distinct_strcs[i]]]=(sequences_path[i])

    with open(root_path+"/labels.json", 'w') as js_file:
         json.dump(path_data,js_file)

    for i in range(gen) :
        df = pd.read_csv(root_path+'/genselected'+str(i)+'.csv')
        strcs = list(set(df['structure'].values))
        for s in ds :
            if s in strcs :
                time_ds[s].append(i)

        print(i)

    for i in range(len(distinct_strcs)) :
        k = distinct_strcs[i]
        plt.plot(time_ds[k],[Y[::-1][i]]*len(time_ds[k]),'|',color='red')
        plt.text(time_ds[k][-1]+2,Y[::-1][i], s=labels[k], fontsize=7)

    plt.plot(np.concatenate(data)[::-1], Y_axis)
    plt.savefig("../images/rnafold_continuity.pdf")
    plt.savefig("../images/rnafold_continuity.png")
    plt.show()





if __name__ =="__main__" :
    main()
