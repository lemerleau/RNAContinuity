

import RNA
import numpy as np
import utility as ut
from pandas import DataFrame
import Landscape
import matplotlib.pyplot as plt
import os
import multiprocess as mp
import uuid

def fold_with_nussinov(sequence) :
    cmd = "python RAFFT/nussinov.py -s {}"
    p = os.popen(cmd.format(sequence))
    rst = p.read().split()
    if len(rst) > 0 :
        #print(rst)
        return rst[2]
    else :
        print("Error during the nussinov folding")
        return None


def ppnussinov_fold (sequences) :
    result = []
    for seq in sequences :
        result += [fold_with_nussinov(seq)]

    #pool = mp.Pool(mp.cpu_count())
    #result = np.array(pool.map(fold_with_nussinov, sequences))
    #pool.close()
    return result


def get_bp_position(structure) :
    pk_pairs, bp, nbp= [] , [], []

    pairs = {
        'bp' : [],
        'pk' : [],
        'nbp' : []
    }
    for i,elt in enumerate(structure) :
        if elt =='[' :
            pk_pairs.append(i)
        elif elt ==']' :
            pairs['pk'] += [(pk_pairs.pop(),i)]

        elif elt == '(' :
            bp += [i]
        elif elt == ')' :
            pairs['bp'] += [(bp.pop(),i)]
        else :
            pairs['nbp'] += [i]
    return pairs



def init(pop_size, seq_length, nucleotide_set,landscape) :

    pop = [RNA.random_string(seq_length,''.join(nucleotide_set))]*pop_size

    pop_str = [fold_with_nussinov(pop[0])]*pop_size
    fitness = [RNA.hamming_distance(pop_str[0],landscape.target)]*pop_size

    return {pop[0]:(pop_str[0],RNA.hamming_distance(landscape.target,pop_str[0]))},DataFrame(np.array([pop,pop_str,fitness]).T, columns=['sequence', 'structure','fitness'])


def reproduce(pop,new_pop_strc, fitnesses) :

    weights  = 1./(1+np.array(fitnesses, dtype = float))
    selected = np.random.choice(range(len(pop)),size=len(pop),p=weights/sum(weights))
    new_seqs = np.array(pop)
    return DataFrame(np.array([new_seqs[selected],np.array(new_pop_strc)[selected],np.array(fitnesses)[selected]]).T,
                                                        columns=['sequence', 'structure','fitness'])

def nthHarmonic(N,s) :

    harmonic = 1.00
    for i in range(2, N + 1) :
        harmonic += 1 / i**s

    return harmonic

def zipf_rvgen(low, high, N, size, a=6.) :
    choices = np.array(range(low, high+1))
    probs = choices**(-a) / nthHarmonic(N,a)
    return np.random.choice(choices,p=np.array(probs)/sum(probs),size=size)

def gen_point_mutation_dist(size,pos,c):

    bp_pos = len(pos['bp']+pos['pk'])
    nbp_pos = len(pos['nbp'])
    dist = {}
    if c !=None :
        if 0<c<7.5:
            if bp_pos > 0 :
                dist['bp'] = zipf_rvgen(1,bp_pos, bp_pos, size, c)
            else :
                dist['bp'] = []
            dist['nbp'] = zipf_rvgen(1,nbp_pos, nbp_pos, size, c)
    else :
        if bp_pos > 0 :
            dist['bp'] = np.ones(size, dtype=int)
        else :
            dist['bp'] = []
        dist['nbp'] = np.ones(size, dtype=int)

    return dist



def mutate(pop, pos,  distribution) :
    base_paire = ["GC","CG","AU","UA", "GU", "UG"]
    nucleotides = ["A", "G","U","C"]
    bp_pos = pos["bp"] + pos['pk']
    nbp_indexes = pos['nbp']
    mutants = []

    bp_points = distribution['bp']
    nb_points = distribution['nbp']

    for i in range(len(pop)) :
        mutant = np.array(list(pop[i]))

	    #Mutate a non base pair position
        nbp_indexes_to_mutate = np.random.choice(list(nbp_indexes), nb_points[i], replace=False)
        mutant[nbp_indexes_to_mutate]= np.random.choice(nucleotides, nb_points[i])


        #Mutate a base pair position
        if len(bp_pos) > 0 :

            bp_choice = np.random.choice(base_paire,bp_points[i])
            bp_indexes_to_mutate = np.random.choice(range(len(bp_pos)), bp_points[i],replace=False)
            bp_to_mutate = np.array(bp_pos)[bp_indexes_to_mutate]

            for j in range(bp_points[i]) :
                mutant[bp_to_mutate[j][0]] =bp_choice[j][0]
                mutant[bp_to_mutate[j][1]] =bp_choice[j][1]

        mutants.append("".join(mutant))

    return mutants


def crossover(parent1, parent2) :
    return

def evolution(params) :

    current_pop = params['init_pop']
    db = params['db']
    mean_fitness = [np.mean(np.array(current_pop["fitness"], dtype = float))]
    max_fitness = [np.min(np.array(current_pop["fitness"], dtype = float))]
    mu_dist = gen_point_mutation_dist(params['pop_size'], params['pos'], params['c'])

    for i in range(params['time']) :
        if i%1 == 0 :
            print('Generation {}, and mean fitness {}'.format(i,mean_fitness[i]))

        #mutation
        mutated = mutate(current_pop['sequence'].values, params['pos'], mu_dist)

        #Evaluation
        #new_pop_strc = ut.ppfold(mutated)[0]
        new_pop_strc = ppnussinov_fold(mutated)
        new_fitnesses = ut.pphamming(new_pop_strc, params['landscape'])
        current_pop = reproduce(mutated,new_pop_strc,new_fitnesses)
        mean_fitness.append(np.mean(np.array(current_pop["fitness"], dtype = float)))
        max_fitness.append(np.min(np.array(current_pop["fitness"], dtype = float)))
        #uncomment if you would to break once the target is found
        if 0. in np.array(current_pop["fitness"], dtype = float) :
            print(current_pop[current_pop['fitness']=='0'].values, i)
            break

    return current_pop,mean_fitness, max_fitness



def main() :

    nucleotides = ["A", "U" , "G", "C"]
    pop_size = 100
    target= '..((((((.........))))))..((((.....))))'
    rate = None
    length = len(target)
    ldscape = Landscape.Landscape(target)
    time =50

    print("Starting RNA evolution with target {} of length {}".format(target,len(target)))
    db,init_pop = init(pop_size,length,['A'],ldscape)
    params = {
    'c' : rate,
    "pop_size" : pop_size,
    'init_pop': init_pop,
    'landscape': ldscape,
    'nucleotides': nucleotides,
    'time': time,
    'db': db,
    "pos" : get_bp_position(target)
    }
    print(init_pop)
    data_min = []
    data_mean = []

    for c  in np.arange(1, 2,0.5) :
        params['c'] = c
        print ("Start evolution with c = {}......".format(params['c']))
        pool = mp.Pool(mp.cpu_count())
        results = pool.map(evolution, [params]*3)
        print(results)

    """
    print(data_mean)
    print(data_min)
    figure = plt.figure(constrained_layout=True, figsize=(8,5))
    gs = figure.add_gridspec(nrows=1, ncols=2, left=0.05, right=0.48, wspace=0.05)
    ax = figure.add_subplot(gs[0,0])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)


    plt.title("Mean Fitness")

    ax.plot(data_mean[0][0], color='orange', label=r"$c = 1.5$")
    ax.plot(data_mean[0][1],color='blue', label=r"$c = 6.5$")
    for line in data_mean :
        ax.plot(line[0], color='orange', linewidth=0.3)
        ax.plot(line[1],color='blue', linewidth=0.3)
    plt.legend()
    #plt.plot(var)
    plt.ylabel("Averange distance to the target")
    plt.xlabel(r"Time($t$)")


    ax = figure.add_subplot(gs[0,1])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)


    plt.title("Max Fitness")
    ax.plot(data_min[0][0], color='orange', label=r"$c = 1.5$")
    ax.plot(data_min[0][1],color='blue', label=r"$c = 6.5$")
    for line in data_min :
        ax.plot(line[0], color='orange', linewidth=0.3)
        ax.plot(line[1],color='blue', linewidth=0.3)
    plt.legend()
    #plt.plot(var)
    plt.ylabel("Maximum distance to the target")
    plt.xlabel(r"Time($t$)")
    plt.savefig('../images/Levy_nussinov.pdf')
    #print("the mean fitness plot is saved in images/mean_fitness.pdf")
    plt.show()
    """
    #mutate(init_pop,rate, nucleotides)
    #print(ken_rafft(RNA.random_string(length, 'AUGC')))
    #print(fold_with_nussinov(RNA.random_string(30, 'AUGC')))




if __name__=="__main__" :
    main()
