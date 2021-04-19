"""
  @author: Nono Saha Cyrille Merleau
  @email: nonosaha@mis.mpg.de
"""
#importing necessary libraries
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import RNA



#utility functions
def getBTNum(log_folder, gen, ref_id) :
    data = pd.read_csv(log_folder+"genmutated"+str(gen)+".csv")
    strings = np.array(data['id'].values)
    return data[data['id']==ref_id].values[0][0]

def getParent(log_folder, gen, bt_num) :
    data = pd.read_csv(log_folder+"genselected"+str(gen)+".csv")
    parent = data.values[bt_num]
    return parent

def getSeq_fromID(log_folder, gen,ref_id) :
    data = pd.read_csv(log_folder+"genmutated"+str(gen)+".csv")
    rst = data[data['id']==ref_id]['sequence'].values
    return rst[0]


def backtrackevolution(log_folder,  target, max_gen,ref_id= ''):
    gen = max_gen
    result = [[gen,-1,getSeq_fromID(log_folder, gen,ref_id), target,0.,ref_id]]
    bt_num = getBTNum(log_folder, gen, ref_id)
    print ("Start backtracking evolution from gen {}.......".format(gen))
    ref_index,ref_seq,ref_str, ref_fitness,id_ = getParent(log_folder, gen, bt_num)
    result.append([gen, bt_num, ref_seq,ref_str, ref_fitness,id_])
    print (gen, bt_num, ref_seq,ref_str, ref_fitness,id_)

    for i in range(gen-1, -1, -1) :

        bt_num = getBTNum(log_folder, i, id_)
        ref_index,ref_seq,ref_str, ref_fitness,id_  = getParent(log_folder, i, bt_num)

        result.append([i, bt_num, ref_seq,ref_str, ref_fitness,id_])
        print (i, bt_num, ref_seq,ref_str, ref_fitness)
    print ("backtracking finished")
    return result


def main() :

    print("Ploting...")

    print("1. Average distance to the target")
    root_path = '../data/exp0/'
    mf = []
    for i in range(2000) :
        df = pd.read_csv(root_path+'/genselected'+str(i)+'.csv')
        mf += [np.mean(np.array(df['fitness'].values, dtype=float))]

    """
    ref_id = "910a5780-c9ad-43da-87cf-e22e132b4645"
    ref_seq = "GAGCUACUGGCGCGACAACCGGCGUACGUCCAUGUCACGGGCGAUAUAUGGAAUACUCACUUCCAAUAGCUCUCCU"
    target = '((((((...((((........)))).(((((.......))))).....(((((.......))))).))))))....'
    gen = 2662
    evolutionary_path = backtrackevolution(root_path, target,gen,ref_id)
    linage = np.array(evolutionary_path)
    data = {}
    ds = set(linage[:,3])


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
    print("=========", len(evolutionary_path))
    for i in range(len(evolutionary_path))[::-1]:
        if evolutionary_path[i][3] not in init_ts.keys() :
            print(RNA.hamming_distance(evolutionary_path[i][3],target))
            init_ts[evolutionary_path[i][3]] = gen-i

    green_data = []
    ep_time = sorted(init_ts.values())
    Y_axis = []
    for i in range(len(ep_time)-1) :
        green_data +=range(ep_time[i], ep_time[i+1])
        Y_axis +=[46-i*1.8]*len(range(ep_time[i], ep_time[i+1]))

    ##### Ploting...#######
    fig, ax = plt.subplots()
    """
    plt.plot(mf,color='black')
    #plt.plot(var)
    plt.ylabel("Averange distance to the target")
    plt.xlabel(r"Time($t$)")
    """
    i = 0
    for cplt in init_time :
        k = cplt[0]
        plt.plot(time_ds[k],[46-i*1.8]*len(time_ds[k]),'|',color='red')
        i = i+1

    plt.plot(green_data, Y_axis, color="green")
    plt.savefig('../images/continuity.pdf')
    """
    plt.show()





if __name__=="__main__" :
    main()
