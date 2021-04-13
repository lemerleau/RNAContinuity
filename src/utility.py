
"""


"""
import pandas
import random
import numpy
import os
import subprocess
import multiprocess
import time
import scipy
import argparse
import RNA
import uuid
from folding_wrapper import *


def pphamming(listOfStructures, landscape) :
    pool = multiprocess.Pool(multiprocess.cpu_count())
    dists = pool.map(landscape.fitness,listOfStructures)
    pool.close()
    return dists

def ppens_defect(listOfSequences, landscape) :
    pool = multiprocess.Pool(multiprocess.cpu_count())
    dists = pool.map(landscape.ens_defect,listOfSequences)
    pool.close()
    return dists

def get_bp_position(structure) :
    position = RNA.ptable(structure)
    position = list(position[1:])

    base_paire_pos = []
    for i in range(len(position)) :
        if position[i] != 0 :
            if (position[i]-1,i) in base_paire_pos :
                continue;
            else :
                base_paire_pos.append((i,position[i]-1))

    return base_paire_pos

#Logging population
def bt_save_population(prev_pop, population,gen, root_path) :
    data = []
    prev_data = []
    for i in range(len(population)) :
        data.append([population[i].RNA_seq, population[i].RNA_structure,population[i].mfe, population[i].fitness])
        prev_data.append([prev_pop[i].RNA_seq, prev_pop[i].RNA_structure,prev_pop[i].mfe, prev_pop[i].fitness])


    dataFrame = pandas.DataFrame(data)
    prev_dataFrame = pandas.DataFrame(prev_data)
    prev_dataFrame.to_csv(root_path+"/prev_gen"+str(gen)+".csv")
    dataFrame.to_csv(root_path+"/gen"+str(gen)+".csv")


def save_population(population,gen, root_path, file_prefix) :

    population.to_csv(root_path+"/gen"+file_prefix+str(gen)+".csv")



#Evaluate the energy
def ppeval(listOfSeqs, target) :
    task = uuid.uuid4()
    with open("tmp/rnaeval_in"+str(task), "w") as file_ :
        for seq in listOfSeqs :
            file_.write(seq.strip()+"\n"+target.strip()+"\n")
        file_.close()

    os.system("RNAeval -j --infile=tmp/rnaeval_in"+str(task)+" |tr -d A-Z,'(',')'|cut -d ' ' -f 2- > tmp/result_"+str(task))
    with open("tmp/result_"+str(task), "r") as file_ :
        eval_ = file_.read().split()
    os.remove("tmp/result_"+str(task))
    os.remove("tmp/rnaeval_in"+str(task))
    return list(numpy.array(eval_, dtype=float))



def ppfold(listOfSeqs,tool) :

    if tool =="v" :
        return ppRNAfold(listOfSeqs)
    if tool =="l" :
        return pplinearFold(listOfSeqs)
    if tool=="c" :
        return ppcontextFold(listOfSeqs)
