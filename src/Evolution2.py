
import RNA
import numpy as np
import utility as ut
from pandas import DataFrame
import Landscape
import matplotlib.pyplot as plt
import os
import multiprocess as mp
import uuid

def fold_with_rafft(sequence) :
    cmd = "python RAFFT/rafft.py -s {}"
    p = os.popen(cmd.format(sequence))
    rst = p.read().split()
    if len(rst) > 0 :
        return rst[1]
    else :
        print("Error during the folding")
        return None

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

def ken_rafft(sequence) :
    id_ = uuid.uuid4()
    cmd = "python RAFFT/rafft.py -s {} -ms 20 --verbose > {}.out"
    p_rafft = os.popen(cmd.format(sequence,id_))
    p_rafft.close()
    ken_cmd = "python RAFFT/utility/kinetic.py {}.out -ns 10000 "
    p = os.popen(ken_cmd.format(id_))
    rst = p.read().split('\n')
    p.close()

    rst.remove("")
    os.remove('{}.out'.format(id_))
    return rst[-1].split()[0]



def pprafft_fold (sequences) :
    pool = mp.Pool(mp.cpu_count())
    result = numpy.array(pool.map(ken_rafft, sequences))
    pool.close()
    return result

def init(pop_size, seq_length, nucleotide_set,landscape) :

    pop = [RNA.random_string(seq_length,''.join(nucleotide_set))]*pop_size
    #RNA.read_parameter_file('rna_turner1999.par')
    #pop_str = [RNA.fold(pop[0])[0]]*pop_size

    pop_str = [ken_rafft(pop[0])]*pop_size
    fitness = [RNA.hamming_distance(pop_str[0],landscape.target)]*pop_size
    entropies = [entropy(pop[0])]*pop_size

    return {pop[0]:(pop_str[0],RNA.hamming_distance(landscape.target,pop_str[0]),entropies[0])},DataFrame(np.array([pop,pop_str,fitness,entropies]).T, columns=['sequence', 'structure','fitness','entrop'])

def entropy(s) :
    fc = RNA.fold_compound(s)
    fc.pf()
    P_ij= np.array(fc.bpp())
    P_ij = P_ij[P_ij>0]
    return -np.sum(P_ij*np.log(P_ij))


def reproduce(pop,new_pop_strc, fitnesses, alpha=0.5) :
    if alpha == 0 :
        weights  = 1./(0.01+np.array(fitnesses[0], dtype = float)/len(pop[0]))
    else :
        fits1  = 1./(0.01+np.array(fitnesses[0], dtype = float)/len(pop[0]))
        fits1 = fits1/sum(fits1)
        fits2  = np.array(fitnesses[1], dtype = float)/len(pop[0])

        fits2 = fits2/sum(fits2)

        weights = alpha*fits1 + (1-alpha)*fits2
    selected = np.random.choice(range(len(pop)),size=len(pop),p=weights/sum(weights))
    new_seqs = np.array(pop)
    return DataFrame(np.array([new_seqs[selected],np.array(new_pop_strc)[selected],np.array(fitnesses[0])[selected],np.array(fitnesses[1])[selected]]).T,
                                                        columns=['sequence', 'structure','fitness','entrop'])



def mutate(pop, rate,nucleotides) :

    mutants = []
    for i in range(len(pop)) :
        mutant = np.array(list(pop[i]))
        nb_mut = np.random.poisson(rate*len(mutant))
        nbp_indexes_to_mutate = np.random.choice(range(len(mutant)), nb_mut, replace=False)
        mutant[nbp_indexes_to_mutate]= np.random.choice(nucleotides, nb_mut)

        mutants.append("".join(mutant))

    return mutants

def mutate2(pop, rate,nucleotides) :

    mutants = []
    for i in range(len(pop)) :
        mutant = np.array(list(pop[i]))
        nb_mut = np.random.rand(len(mutant))
        nbp_indexes_to_mutate = np.where(nb_mut<0.001)
        if len(nbp_indexes_to_mutate[0])>0 :
            for idx in nbp_indexes_to_mutate[0] :
                mutant[idx]= np.random.choice(list(set(nucleotides)-set([mutant[idx]])), 1)[0]
                break

        mutants.append("".join(mutant))

    return mutants

def crossover(parent1, parent2) :
    return

def evolution(params) :

    current_pop = params['init_pop']
    db = params['db']
    mean_fitness = [np.mean(np.array(current_pop["fitness"], dtype = float))]
    mean_entrop = [np.mean(np.array(current_pop["entrop"], dtype = float))]

    for i in range(params['time']) :
        #ut.save_population(current_pop,i,'../data/exp0','selected')
        if i%1 == 0 :
            print('Generation {}, and mean fitness {}, mean entropy {}'.format(i,mean_fitness[i],mean_entrop[i]))

        #mutation
        mutated = mutate(current_pop['sequence'].values, params['rate'], params['nucleotides'])

        #Evaluate
        new_pop_strc = []
        new_fitnesses = []
        new_entropies = []
        for s in mutated :
            try:
                ss = db[s]
            except Exception as e:
                #ss = fold_with_nussinov(s)
                #RNA.read_parameter_file('rna_turner1999.par')
                ss = RNA.fold(s)[0]
                db[s] = (ss,RNA.hamming_distance(ss,params['landscape'].target), 0)

                ss = db[s]
            new_pop_strc += [ss[0]]
            new_fitnesses += [float(ss[1])]
            new_entropies += [float(ss[2])]
        #print(len(new_pop_strc))
        #new_pop_strc = [RNA.fold(seq)[0] for seq in mutated]
        #new_pop_strc = ut.ppfold(mutated)[0]
        #new_pop_strc = pprafft_fold(mutated)
        #new_fitnesses = ut.pphamming(new_pop_strc, params['landscape'])
        #ut.save_population(DataFrame(np.array([mutated,new_pop_strc,new_fitnesses, new_entropies]).T,columns=['sequence', 'structure','fitness','entrop']),i,'../data/exp0','mutated')
        current_pop = reproduce(mutated,new_pop_strc,[new_fitnesses,new_entropies],alpha=0.)
        mean_fitness.append(np.mean(np.array(current_pop["fitness"], dtype = float)))
        mean_entrop += [np.mean(np.array(current_pop["entrop"], dtype = float))]
        #uncomment if you would to break once the target is found
        if 0. in np.array(current_pop["fitness"], dtype = float) :
            break

    return current_pop,mean_fitness



def main() :

    nucleotides = ["A", "U" , "G", "C"]
    pop_size = 1000
    target= '(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))....'
    rate = 0.001
    length = len(target)
    ldscape = Landscape.Landscape(target)
    time =2000

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
    #plt.savefig('../images/mean_fitness.pdf')
    #print("the mean fitness plot is saved in images/mean_fitness.pdf")
    plt.show()
    #mutate(init_pop,rate, nucleotides)
    #print(ken_rafft(RNA.random_string(length, 'AUGC')))
    #print(fold_with_nussinov(RNA.random_string(30, 'AUGC')))




if __name__=="__main__" :
    main()
