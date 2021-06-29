

import RNA
import numpy as np
import utility as ut
from pandas import DataFrame
import Landscape
import matplotlib.pyplot as plt


def init(pop_size, seq_length, nucleotide_set,landscape) :

    pop = [RNA.random_string(seq_length,''.join(nucleotide_set))]*pop_size
    pop_str = [RNA.fold(pop[0])[0]]*pop_size
    fitness = ut.pphamming(pop_str,landscape)

    return {pop[0]:(pop_str[0],RNA.hamming_distance(landscape.target,pop_str[0]))},DataFrame(np.array([pop,pop_str,fitness]).T, columns=['sequence', 'structure','fitness'])


def reproduce(pop,new_pop_strc, fitnesses) :
    weights  = 1./(0.01+np.array(fitnesses, dtype = float)/len(pop[0]))
    selected = np.random.choice(range(len(pop)),size=len(pop),p=weights/sum(weights))
    new_seqs = np.array(pop)
    return DataFrame(np.array([new_seqs[selected],np.array(new_pop_strc)[selected],np.array(fitnesses)[selected]]).T,
                                                        columns=['sequence', 'structure','fitness'])


def mutate(pop, rate,nucleotides) :

    mutants = []
    for i in range(len(pop)) :
        mutant = np.array(list(pop[i]))
        nb_mut = np.random.poisson(rate*len(mutant))
        nbp_indexes_to_mutate = np.random.choice(range(len(mutant)), nb_mut, replace=False)
        mutant[nbp_indexes_to_mutate]= np.random.choice(nucleotides, nb_mut)

        mutants.append("".join(mutant))

    return mutants

def crossover(parent1, parent2) :
    return

def evolution(params) :

    current_pop = params['init_pop']
    db = params['db']
    mean_fitness = [np.mean(np.array(current_pop["fitness"], dtype = float))]

    for i in range(params['time']) :

        if i%10 == 0 :
            print('Generation {}, and mean fitness {}'.format(i,mean_fitness[i]))

        #mutation
        mutated = mutate(current_pop['sequence'].values, params['rate'], params['nucleotides'])

        #Evaluate
        new_pop_strc = []
        new_fitnesses = []
        for s in mutated :
            try:
                ss = db[s]
            except Exception as e:
                ss = RNA.fold(s)[0]
                db[s] = (ss,RNA.hamming_distance(ss,params['landscape'].target))
                ss = db[s]
            new_pop_strc += [ss[0]]
            new_fitnesses += [float(ss[1])]
        #print(len(new_pop_strc))
        #new_pop_strc = [RNA.fold(seq)[0] for seq in mutated]
        #new_pop_strc = ut.ppfold(mutated)[0]

        #new_fitnesses = ut.pphamming(new_pop_strc, params['landscape'])
        current_pop = reproduce(mutated,new_pop_strc,new_fitnesses)

        mean_fitness.append(np.mean(np.array(current_pop["fitness"], dtype = float)))

        #uncomment if you would to break once the target is found
        #if 0. in np.array(current_pop["fitness"], dtype = float) :
            #break

    return current_pop,mean_fitness



def main() :

    nucleotides = ["A", "U" , "G", "C"]
    pop_size = 1000
    target= '((((((...((((........)))).(((((.......))))).....(((((.......))))).))))))....'
    #target = ".((((....((((((((....)))).((((........)))).....((((((....)))))).))))))))...."
    rate = 0.001
    length = len(target)
    ldscape = Landscape.Landscape(target)
    time = 300

    print("Starting RNA evolution with target of length {}".format(len(target)))
    db,init_pop = init(pop_size,length,nucleotides,ldscape)
    params = {
    'rate' : rate,
    'init_pop': init_pop,
    'landscape': ldscape,
    'nucleotides': nucleotides,
    'time': time,
    'db': db
    }
    print(init_pop)
    print ("Begin mutation......")
    best_pop,mf = evolution(params)
    print(mf)
    plt.plot(mf)
    plt.xlabel(r'Time($t$)')
    plt.ylabel('Mean Fitness')
    plt.savefig('../images/mean_fitness.pdf')
    print("the mean fitness plot is saved in images/mean_fitness.pdf")
    plt.show()
    #mutate(init_pop,rate, nucleotides)




if __name__=="__main__" :
    main()
