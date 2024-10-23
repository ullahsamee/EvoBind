import numpy as np
import copy
import pdb

def search_aa(seqlen, peptide_sequence, searched_seqs):
    '''
    '''
    restypes = np.array(['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V' ])

    #Go through a shuffled version of the positions and aas
    for pi in np.random.choice(np.arange(seqlen),seqlen, replace=False):
        for aa in np.random.choice(restypes,len(restypes), replace=False):
            new_seq = copy.deepcopy(peptide_sequence)
            new_seq = new_seq[:pi]+aa+new_seq[pi+1:]
            if new_seq in searched_seqs:
                #Recursion - next will use new seq as peptide_sequence
                search_aa(seqlen, new_seq, searched_seqs)
            else:
                return new_seq

def mutate_sequence(peptide_sequence, sequence_scores):
    '''Mutate the amino acid sequence randomly
    '''

    restypes = np.array(['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V' ])

    seqlen = len(peptide_sequence)
    searched_seqs = sequence_scores['sequence']
    #Mutate seq
    seeds = [peptide_sequence]
    #Go through a shuffled version of the positions and aas
    for seed in seeds:
        for pi in np.random.choice(np.arange(seqlen),seqlen, replace=False):
            for aa in np.random.choice(restypes,len(restypes), replace=False):
                new_seq = copy.deepcopy(seed)
                new_seq = new_seq[:pi]+aa+new_seq[pi+1:]
                if new_seq in searched_seqs:
                    continue
                else:
                    # Check if the sequence contains exactly two 'C's
                    if new_seq.count('C') == 2:
                        # Ensure the two 'C's are not adjacent to each other
                        c_indices = [i for i, x in enumerate(new_seq) if x == 'C']
                        if abs(c_indices[0] - c_indices[1]) > 1:
                            return new_seq
                        else:
                            seeds.append(new_seq)

        seeds.append(new_seq)
