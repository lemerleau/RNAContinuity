'''
    @Author: Nono Saha Cyrille Merleau
    @Email:  nonosaha@mis.mpg.de
'''



#Here are necessary libraries to import 
import numpy 
import RNA as rna 
import random 
import matplotlib.pyplot as plt
import RNAEvolution
import Landscape
import pandas 
import os

def getNumberGen(folder) : 
    try : 
        files = os.listdir(folder)
        min_gen = len(files) 
    except : 
        print "Folder doesn't exist"
    return min_gen/2


def getBTNum(log_folder, gen, ref_seq) : 
    data = pandas.read_csv(log_folder+"gen"+str(gen)+".csv")
    bt_num = None 
    for ind in data.values : 
        if str(ind[1]) == ref_seq : 
            bt_num = ind[0] 
        
    return bt_num 

def getParent(log_folder, gen, bt_num) :
    data = pandas.read_csv(log_folder+"prev_gen"+str(gen)+".csv")
    return data.values[bt_num][1], data.values[bt_num][2],data.values[bt_num][4]

def backtrackevolution(log_folder,  target, found_sequence= 'GGUCCGCCAGCUCUUGC'): 
    gen = getNumberGen(log_folder) 
    print gen
    result = [[gen,None,found_sequence,target,1.0]]
    bt_num = getBTNum(log_folder, gen-1, found_sequence)
    print bt_num
    print "Start backtracking the evolution......."
    ref_seq,ref_str, ref_fitness = getParent(log_folder, gen-1, bt_num)
    result.append([gen-1, bt_num, ref_seq,ref_str, ref_fitness])
    print (gen-1, bt_num, ref_seq,ref_str, ref_fitness)
    for i in range(gen-2, -1, -1) : 
    
        bt_num = getBTNum(log_folder, i, ref_seq)
        print bt_num, i, ref_seq
        ref_seq,ref_str, ref_fitness = getParent(log_folder, i, bt_num)
        result.append([i, bt_num, ref_seq,ref_str, ref_fitness])
        print (i, bt_num, ref_seq,ref_str, ref_fitness)
    print "backtracking finished"
    return result

def getMeanFitness(log_folder): 
    gen = getNumberGen(log_folder) 
    means = []
    for i in range(gen): 
        data = pandas.read_csv(log_folder+"/gen"+str(i)+".csv")
        means.append(numpy.mean(data.values[:, 4]))
    return means

def getWinner(log_folder): 
    gen = getNumberGen(log_folder)
    data = pandas.read_csv(log_folder+"gen"+str(gen-1)+".csv") 

    for row in data.values : 
        if row[4] == 1.0 : 
            return row[2],row[1]
    return None, None
def getEnsDiversity(log_folder): 
    gen = getNumberGen(log_folder) 
    means = []
    EnsDivs = []
    for i in range(gen): 
        data = pandas.read_csv(log_folder+"/gen"+str(i)+".csv")
        ens_div = []
        for seq in data.values[:, 1] : 
            #print seq
            os.system("echo "+ str(seq)+"| RNAfold -p |tail -1|cut -d ' ' -f 11 > rest.txt" )
            dt = pandas.read_csv('rest.txt', header=None)
            ens_div.append(dt.values[0][0])
        means.append(numpy.mean(ens_div))
        EnsDivs.append(ens_div)
    return EnsDivs, means

def ComputeEnsDiversity(listOfSeq): 
    ens_div = []
    for seq in listOfSeq : 
        cmd = "echo "+str(seq)+"| RNAfold -p |tail -1|cut -d ' ' -f 11 > rest.txt" 
        os.system(cmd)
        dt = pandas.read_csv('rest.txt', header=None)
        ens_div.append(dt.values[0][0])
    return ens_div

#Main function 
def main() : 
    folders = range(0,31,1)
    success = len(folders)
    for fold in folders : 
        target, found_sequence = getWinner("../Logs/BenchMark/BT/15/"+str(fold)+"/1/")
        print "Found sequence == ", found_sequence, target
        if found_sequence == None : 
            success = success - 1 
            continue

        
        #bt_fitness = numpy.array(backtrackevolution("../Logs/BenchMark/BT/44/"+str(fold)+"/0/", target=target, found_sequence=found_sequence))[:,4]
        bt_sequence = numpy.array(backtrackevolution("../Logs/BenchMark/BT/15/"+str(fold)+"/1/", target=target, found_sequence=found_sequence))[:,2]
        bt_ED = ComputeEnsDiversity(bt_sequence)

        #means = getMeanFitness("../Logs/BenchMark/BT/44/"+str(fold)+"/0/")
        EnsDivs, means = getEnsDiversity("../Logs/BenchMark/BT/15/"+str(fold)+"/1/")
        #plt.subplot(1,1, 1)
        #plt.plot(means, '-b',label="Mean Ens Div")
        #plt.plot(numpy.flip(bt_ED), '-r', label="Back tracking Ens Div")
        #plt.ylabel("Fitness", fontweight="bold")
        #plt.xlabel("Generation", fontweight="bold")
        #plt.legend(loc="upper left", shadow=True, fontsize='5')
        #plt.savefig("../Logs/BenchMark/BT/33/"+str(fold)+"/F4_bt_evolution.eps")
        #plt.savefig("../Logs/BenchMark/BT/33/"+str(fold)+"/F4_bt_evolution.png")
        #plt.clf()
        
        #Compute and plot the rank of ensemble diversity 
        from scipy.stats import rankdata
        bt_ED = numpy.flip(bt_ED)
        bt_ED = bt_ED[1:]
        bt_ED_ranks = []
        mean_ranks = []
        print len(bt_ED), len(EnsDivs)
        for i in range(len(bt_ED)) : 
            ranks = rankdata(EnsDivs[i])
            mean_ranks.append(numpy.mean(ranks))
            bt_ED_ranks.append(ranks[list(EnsDivs[i]).index(bt_ED[i])])

        plt.plot(mean_ranks)
        plt.plot(bt_ED_ranks)

    #plt.savefig("../Logs/BenchMark/BT/15/ED_bt_evolution.eps")
    #plt.savefig("../Logs/BenchMark/BT/15/ED_bt_evolution.png")
    plt.show()

    print "Backtracked on ", success, " success"

if __name__== "__main__" : 

    main()
