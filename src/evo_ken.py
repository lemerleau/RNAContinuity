

import RNA
from numpy import array, mean, random, log,where
import utility as ut
from pandas import DataFrame, read_csv
import Landscape
from os import popen,remove, listdir, path, mkdir
from multiprocess import Pool,cpu_count
from uuid import uuid4
from argparse import ArgumentParser,RawTextHelpFormatter

def fold_with_rafft(sequence) :
    cmd = "python RAFFT/rafft.py -s {}"
    p = popen(cmd.format(sequence))
    rst = p.read().split()
    if len(rst) > 0 :
        return rst[1]
    else :
        print("Error during the folding")
        return None

def ken_rafft(sequence) :
    id_ = uuid4()
    cmd = "python RAFFT/rafft.py -s {} -ms 20 --verbose > tmp/{}.out"
    p_rafft = popen(cmd.format(sequence,id_))
    p_rafft.close()
    ken_cmd = "python RAFFT/utility/kinetic.py tmp/{}.out -ns 10000 "
    p = popen(ken_cmd.format(id_))
    rst = p.read().split('\n')
    p.close()

    rst.remove("")
    remove('tmp/{}.out'.format(id_))
    return rst[-1].split()[0]


def pprafft_fold (sequences) :
    pool = Pool(cpu_count())
    result = array(pool.map(ken_rafft, sequences))
    pool.close()
    return result

def init(pop_size, seq_length, nucleotide_set,landscape) :

    pop = [RNA.random_string(seq_length,''.join(nucleotide_set))]*pop_size
    #RNA.read_parameter_file('rna_turner1999.par')
    #pop_str = [RNA.fold(pop[0])[0]]*pop_size

    pop_str = [ken_rafft(pop[0])]*pop_size
    fitness = [RNA.hamming_distance(pop_str[0],landscape.target)]*pop_size
    #entropies = [entropy(pop[0])]*pop_size
    ids = [uuid4() for i in range(pop_size)]

    return {pop[0]:(pop_str[0],RNA.hamming_distance(landscape.target,pop_str[0]))},DataFrame(array([pop,pop_str,fitness,ids]).T, columns=['sequence', 'structure','fitness','id'])

def entropy(s) :
    fc = RNA.fold_compound(s)
    fc.pf()
    P_ij= array(fc.bpp())
    P_ij = P_ij[P_ij>0]
    return -sum(P_ij*log(P_ij))


def reproduce(population, size) :
    fitnesses = population['fitness'].values
    seqs = population['sequence'].values
    weights  = 1./(0.01+array(fitnesses, dtype = float)/len(seqs[0]))
    selected = list(random.choice(range(len(seqs)),size=size,p=weights/sum(weights)))+list(where(fitnesses == min(fitnesses))[0][:int(size*0.05)])

    selected = array(selected)
    #print(where(fitnesses == min(fitnesses))[0][:5])
    #print(sorted(fitnesses)[::1], selected)
    ids = population['id'].values
    pop_strc = population['structure'].values
    return DataFrame(array([seqs[selected],pop_strc[selected],fitnesses[selected],ids[selected]]).T,columns=['sequence', 'structure','fitness','id'])


def poisson_mutation(pop, rate,nucleotides) :

    mutants = []
    for i in range(len(pop)) :
        mutant = array(list(pop[i]))
        nb_mut = random.poisson(rate*len(mutant))
        nbp_indexes_to_mutate = random.choice(range(len(mutant)), nb_mut, replace=False)
        mutant[nbp_indexes_to_mutate]= random.choice(nucleotides, nb_mut)

        mutants.append("".join(mutant))

    return mutants

def uniform_mutation(pop, rate,nucleotides) :

    mutants = []
    seqs = pop['sequence'].values
    ids = list(pop['id'].values)
    for i in range(len(seqs)) :
        mutant = array(list(seqs[i]))
        nb_mut = random.rand(len(mutant))
        nbp_indexes_to_mutate = where(nb_mut<0.001)
        if len(nbp_indexes_to_mutate[0])>0 :
            for idx in nbp_indexes_to_mutate[0] :
                mutant[idx]= random.choice(list(set(nucleotides)-set([mutant[idx]])), 1)[0]
                ids[i] = uuid4()
                break

        mutants.append("".join(mutant))

    return mutants, ids

def crossover(parent1, parent2) :
    return

def evolution(params) :

    current_pop = params['init_pop']
    db = params['db']
    mean_fitness = [mean(array(current_pop["fitness"], dtype = float))]
    #mean_entrop = [mean(array(current_pop["entrop"], dtype = float))]

    for i in range(params['init_gen'],params['time']) :

        if i%1 == 0 :
            print('Generation {}, and mean fitness {}'.format(i,mean_fitness[i-params['init_gen']]))

        #mutation
        mutated, ids = uniform_mutation(current_pop, params['rate'], params['nucleotides'])

        #Evaluate
        if params['ft'] == 'vrn' :
            #new_pop_strc = ut.ppfold(mutated)[0] #uncomment it to fold the population in parallel.

            new_pop_strc = []
            new_fitnesses = []
            for s in mutated :
                try:
                    ss = db[s]
                except Exception as e:
                    #ss = ken_rafft(s)
                    #RNA.read_parameter_file('rna_turner1999.par')
                    ss = RNA.fold(s)[0]
                    db[s] = (ss,RNA.hamming_distance(ss,params['landscape'].target))

                    ss = db[s]
                new_pop_strc += [ss[0]]
                new_fitnesses += [float(ss[1])]
        else :
            new_pop_strc = pprafft_fold(mutated)
            new_fitnesses = ut.pphamming(new_pop_strc, params['landscape'])

        #print(len(new_pop_strc))
        #new_pop_strc = [RNA.fold(seq)[0] for seq in mutated]


        mutated_pop =DataFrame(array([mutated,new_pop_strc,new_fitnesses,ids]).T,columns=['sequence', 'structure','fitness','id'])
        ut.save_population(mutated_pop,i,'../data/ken_rafft/'+str(params['lf']),'mutated')
        current_pop = reproduce(mutated_pop, params['pop_size'])
        ut.save_population(current_pop,i,'../data/ken_rafft/'+str(params['lf']),'selected')
        mean_fitness.append(mean(array(current_pop["fitness"], dtype = float)))
        #mean_entrop += [mean(array(current_pop["entrop"], dtype = float))]
        #uncomment if you would to break once the target is found
        if 0. in array(current_pop["fitness"], dtype = float) :
            break

    return current_pop,mean_fitness


def parse_arguments():

    parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    parser.add_argument('--target', '-t', help="Target secondary structure")
    parser.add_argument('--pop_size', '-n', help="Population size", default=100)
    parser.add_argument('--mutation_rate', '-mu', help="Peer nucleotide mutation rate", default=0.001)
    parser.add_argument('--time', '-g', help="Number of generations", default=100)
    parser.add_argument('--log_folder', '-lf', help="Folder name where to save the populations")
    parser.add_argument('--folding_tool', '-ft', default = 'vrn', help="Folding tool to use: RNAfold (vrn) or RAFFT+kinetic ansazt (rafft)")

    return parser.parse_args()


def main() :

    args = parse_arguments()

    nucleotides = ["A", "U" , "G", "C"]
    pop_size = int(args.pop_size)
    #target= '((((((...((((........)))).(((((.......))))).....(((((.......))))).))))))....'
    target = str(args.target)
    rate = float(args.mutation_rate)
    length = len(target)
    ldscape = Landscape.Landscape(target)
    time =int(args.time)
    root_dir = "../data/ken_rafft/"
    if not path.exists(root_dir) :
        mkdir(root_dir)
        mkdir(root_dir+str(args.log_folder))
    elif not path.exists(root_dir+str(args.log_folder)):
        mkdir(root_dir+str(args.log_folder))

    print("Starting RNA evolution with target of length {}".format(len(target)))
    files = listdir(root_dir+str(args.log_folder+"/"))
    print(len(files),files)
    nb_sel = len(list(filter(lambda elt: 'genselected' in elt, files)))

    db,init_pop = init(pop_size,length,nucleotides,ldscape)
    params = {
    'rate' : rate,
    'init_pop': init_pop,
    'landscape': ldscape,
    'nucleotides': nucleotides,
    'time': time,
    'db': db,
    'lf' : args.log_folder,
    'pop_size': pop_size,
    'init_gen': 0,
    'ft' : args.folding_tool
    }
    params['init_pop'] = init_pop

    if nb_sel > 0 :
        init_pop2 = read_csv("../data/ken_rafft/"+str(args.log_folder+"/genselected"+str(nb_sel-1)+".csv"))
        params['init_gen'] = nb_sel-1
        params['init_pop'] = init_pop2
        print (init_pop2)

    print(init_pop)
    print ("Begin Evolution......")
    best_pop,mf = evolution(params)



if __name__=="__main__" :
    main()
