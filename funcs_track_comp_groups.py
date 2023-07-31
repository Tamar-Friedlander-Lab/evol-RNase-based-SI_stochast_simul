import numpy as np

def parent_hap_n_iters_before(haplotypes_dict, current_hap, parents_haps):

    hap = current_hap
    if hap in haplotypes_dict:
        #current hap was born before more than 10 generations
        parents_hap = hap
        return parents_hap
    else:
        # current hap was born before less than 10 generations
        if parents_haps[hap] in haplotypes_dict:
            parents_hap = parents_haps[hap]
            return parents_hap
        else:
            parent_not_in_previous_iter = True
            parents_hap = parents_haps[hap]
            while parent_not_in_previous_iter:
                if parents_hap in haplotypes_dict:
                    parent_not_in_previous_iter = False
                    return parents_hap
                else:
                    try:
                        parents_haps[parents_hap]
                    except:
                        return np.nan
                    parents_hap = parents_haps[parents_hap]
                    parent_not_in_previous_iter = True