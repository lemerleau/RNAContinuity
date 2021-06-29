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



#Get in a specific generation a backtracking number that corresponds to an unique id ref_id.
def getBTNum(log_folder, gen, ref_id) :
    data = pd.read_csv(log_folder+"genmutated"+str(gen)+".csv")
    strings = np.array(data['id'].values)
    print(ref_id, gen)
    return data[data['id']==ref_id].values[0][0]

#Get the parent (the selected RNA sequence from which the mutant dirived from) at generation gen and at line bt_num.
def getParent(log_folder, gen, bt_num) :
    data = pd.read_csv(log_folder+"genselected"+str(gen)+".csv")
    parent = data.values[bt_num]
    return parent

#Get an RNA sequence from a population at generation 'gen' where the unique id= ref_id.
def getSeq_fromID(log_folder, gen,ref_id) :
    data = pd.read_csv(log_folder+"genmutated"+str(gen)+".csv")
    rst = data[data['id']==ref_id]['sequence'].values
    return rst[0]

# Get from the log folder, the evolutionnay linage from a ref RNA sequence, target structure.
def backtrackevolution(log_folder,  target, max_gen,ref_id= ''):
    gen = max_gen
    result = [[gen,-1,getSeq_fromID(log_folder, gen,ref_id), target,0.,ref_id]]
    bt_num = getBTNum(log_folder, gen, ref_id)
    print ("Start backtracking evolution from gen {}.......".format(gen))
    ref_index,ref_seq,ref_str, ref_fitness,id_ = getParent(log_folder, gen, bt_num)
    result.append([gen, bt_num, ref_seq,ref_str, ref_fitness,id_])
    print (gen, bt_num, ref_seq,ref_str, ref_fitness,id_)

    for i in range(gen-1, -1, -1) :

        ref_index,ref_seq,ref_str, ref_fitness,id_  = getParent(log_folder, i, bt_num)
        bt_num = getBTNum(log_folder, i, id_)
        result.append([i, bt_num, ref_seq,ref_str, ref_fitness,id_])
        print (i, bt_num, ref_seq,ref_str, ref_fitness, id_)
    print ("backtracking finished")
    return result

def getNumberGen(folder) :
    try :
        files = os.listdir(folder)
        min_gen = len(files)
    except :
        print ("Folder doesn't exist")
    return int(min_gen/2)-1

def getWinner(log_folder, gen):
    data = pd.read_csv(log_folder+"/genselected"+str(gen)+".csv")

    ref = data[data['fitness']==0.0]
    if len(ref)>0:
        return ref.values
    return None, None

def parse_arguments():

    parser = agp.ArgumentParser(description=__doc__, formatter_class=agp.RawTextHelpFormatter)
    #parser.add_argument('--target', '-t', help="Target secondary structure")
    #parser.add_argument('--max_gen', '-mt', help="Max number of generations", default=100)
    #parser.add_argument('--ref_id', '-id', help="Reference id corresponding to the sequence ref_sequence that folds into the target at gen. mt", default=None)
    #parser.add_argument('--ref_sequence', '-rs', help="Reference sequence: the sequence the folds into the target at gen. mt", default=None)
    parser.add_argument('--log_folder', '-lf', help="Folder name where the populations were stored. The folder should be in ../data/ken_rafft/<folder_name>")
    return parser.parse_args()

def main() :
    args = parse_arguments()
    print("Ploting...")

    print("1. Average distance to the target")
    root_path = '../data/ken_rafft/'+str(args.log_folder)
    mf = []
    mxf = []
    mentr = []
    gen = getNumberGen(root_path)
    print('Number of generations is: ', gen)

    ref = getWinner(root_path, gen)[0]
    print("The ref sequence is: ", ref)
    for i in range(gen+1) :
        print("Generation i = ", i)
        df = pd.read_csv(root_path+'/genmutated'+str(i)+'.csv')
        mf += [np.mean(np.array(df['fitness'].values, dtype=float))]
        mxf += [np.min(np.array(df['fitness'].values, dtype=float))]
        #mentr += [np.mean(np.array(df['entrop'].values, dtype=float))]


    ref_id =  ref[-1]
    ref_seq = ref[1]
    target =ref[2]

    #
    #print(evolutionary_path)

    try:
        linage = pd.read_csv(root_path+"/evol_path.csv").values
    except FileNotFoundError as fe:
        print("ERROOOOE")
        evolutionary_path = backtrackevolution(root_path+'/', target,gen,ref_id)
        linage =pd.DataFrame(np.array(evolutionary_path))
        linage.to_csv(root_path+"/evol_path.csv")
        linage = linage.values
    print(len(linage), linage.shape)
    data = {}
    ds = set(linage[:,4])
    print(ds, len(ds))


    print("Loading the time spent by each structure in the evolutionary path...")
    time_ds = {}

    for s in ds :
        time_ds[s]=[]

    for i in range(gen+1) :
        df = pd.read_csv(root_path+'/genmutated'+str(i)+'.csv')
        strcs = list(df['structure'].values)
        for s in ds :
            if s in strcs :
                time_ds[s].append(i)

    init_time = []
    for k in time_ds.keys() :
        init_time +=[(k,time_ds[k])]

    init_time.sort(key=lambda elt: elt[1])


    init_ts = {}
    print("linage length =========", len(linage))
    for i in range(len(linage))[::-1]:
        if linage[i][4] not in init_ts.keys() :
            print(RNA.hamming_distance(linage[i][4],target))
            init_ts[linage[i][4]] = gen-i

    green_data = []
    ep_time = sorted(init_ts.values())
    Y_axis = []
    for i in range(len(ep_time)-1) :
        green_data +=range(ep_time[i], ep_time[i+1])
        Y_axis +=[48-i*1.6]*len(range(ep_time[i], ep_time[i+1]))

    ##### Ploting...#######
    print("Number of green point: ", len(green_data))
    fig, ax = plt.subplots()
    """
    figure = plt.figure(constrained_layout=True, figsize=(6,4))
    gs = figure.add_gridspec(nrows=1, ncols=1, left=0.05, right=0.48, wspace=0.05)
    ax = figure.add_subplot(gs[0,0])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)


    #plt.title("RAFFT")
    ax.plot(mf,color='orange', label="Mean")
    ax.plot(mxf,color='blue', label="Max")
    plt.legend()
    #plt.plot(var)
    plt.ylabel("Averange distance to the target")
    plt.xlabel(r"Time($t$)")

    #ax = figure.add_subplot(gs[0,1])
    #plt.title("RAFFT")
    #ax.plot(mentr,color='blue')
    #plt.ylabel("Averange entropy of sequences")
    #plt.xlabel(r"Time($t$)")
    """
    i = 0
    print(init_time)
    for cplt in init_time :
        k = cplt[0]
        plt.plot(time_ds[k],[48-i*1.6]*len(time_ds[k]),'|',color='red')
        i = i+1

    plt.plot(green_data, Y_axis, color="green")
    plt.plot(mf, color='black')
    plt.savefig('../images/continuity_rafftken.pdf')

    #plt.savefig('../images/entropy.pdf')
    plt.show()





if __name__=="__main__" :
    main()
