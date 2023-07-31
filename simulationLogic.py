from typing import Dict, Any, Tuple, Union
import numpy as np
import pandas as pd
from copy import deepcopy
from itertools import combinations
import simulationVars as Model
import pickle
import UNIPROT_calcs

# ----------------------------------------------------------------
#      INITIATE THE POPULATION
# ----------------------------------------------------------------

def GetRadomAAs():
    return np.random.choice(4, Model.seq_len, p=UNIPROT_calcs.calc_mutation_probs_bio_classes())

def createAllele(bio_class_ar):
    #create one allele object
    allele = []
    indexInRNasesPool = None #index of RNase in RNases pool
    indexInSLFsPool = None
    previousAllele = None

    for i in range(Model.seq_len):
        aaType = Model.AA_letter(bio_class_ar[i])

        aa = Model.AA(aaType)
        allele.append(aa)
    alleleObj = Model.Allele(allele, indexInRNasesPool, indexInSLFsPool, previousAllele)
    return alleleObj


def createAlleles():
    '''
    this function create an alleles- a list of AAs
    output: alleles - an object with the fields:
            al_pool - one allele as a list of 18 AA (18 int values
                      in the range of 0-3 opptional bio-classes.)
    '''
    alleles = []
    AAs_ar = GetRadomAAs()
    alleleObj = createAllele(AAs_ar)
    alleles.append(alleleObj)
    return alleles


def initiate_population_given_prob(haplotypes_dict, prob):

    '''
    'haplotypes_obj' - a list of objects from the type of 'haplotypes'
    'prob' - a list of diploid probabilities. Its length is n_h*(n_h-1)/2

    output:
    'diploid_count_dict' - a dictionary, keys - the two haplotypes indices (in acsending order), values - the amount in population
    'haplotype_count_dict' - a dictionary, keys - the haplotid index, values - the amount in population
    'haplotypes_dict' - the same as in the input the amount of haplotypes in the population was updated
    '''

    num_combs = int(Model.num_types_of_haplotypes * (Model.num_types_of_haplotypes - 1) / 2)
    rand_diploids_indices = np.random.choice(num_combs, Model.popu_size, p = prob)

    #initiate dicts:
    # a dictionary of unique genotypes (diploids) and their amounts
    # each is diploid is a combination of two values: haplotype1 and haplotype2 indices
    diploid_count_dict: Dict[Tuple[Any, Any], Union[int, Any]] = {}
    # a dictionary of unique haplotype indices and their amounts
    haplotype_count_dict: Dict[Tuple[Any], Union[int, Any]] = {}

    unique_diploids_indices, _, _, counts = np.unique(rand_diploids_indices, axis=0, return_index=True, return_inverse=True,
                                                                                                         return_counts=True)
    counts = counts.tolist()
    num_types_of_diploids = len(unique_diploids_indices)
    combs = list(combinations(np.arange(Model.num_types_of_haplotypes), 2))

    for comb in combs:
        diploid_count_dict[comb] = 0

    for i in range(Model.num_types_of_haplotypes):
        haplotype_count_dict[i] = 0

    # update 'haplotype_count_dict' a dict of keys - haplotype indices,
    # values - the number of copies of each haplotype.
    for i in range(num_types_of_diploids):
        comb_ind = combs[unique_diploids_indices[i]]
        diploid_count_dict[comb_ind] = counts[i]
        haplotype_count_dict[comb_ind[0]] = haplotype_count_dict[comb_ind[0]] + counts[i]
        haplotype_count_dict[comb_ind[1]] = haplotype_count_dict[comb_ind[1]] + counts[i]

    for hap_ind in haplotype_count_dict.keys():
        del haplotypes_dict[hap_ind].Male_part
        del haplotypes_dict[hap_ind].Female_part

    return diploid_count_dict, haplotype_count_dict, haplotypes_dict

def initiate_haplotypes():

    # initiate a set of haplotypes, by the following constrains:
    # 1. create a set of RNases by
    # either - SI haplotypes (RNase + non compatible SLFs) and save only the RNases (remove the SLFs)
    # or create only the RNases
    # 2. create SI interactions with a set of SLFs for every RNase, use the constraint by Fully-Cross compatibility
    # either - add SLFs that are interacting with one and only one RNase to the SI haplotypes
    # or set a bank of SLFs, add only non-compatible SLFs to SI haplotypes, do it untill all the haplotypes were fully cross-compatible
    # the set of haplotypes should be built by a combination of one from 1 and one from 2.

    seq_length = Model.seq_len
    all_SLFs_list = []
    all_RNase_list = []
    SLF_object = []
    RNase_object = []
    all_RNase_ar = []
    tmp_RNase = []

    if Model.init_constrained_RNases:
        for n in range(Model.num_types_of_haplotypes):
            SC = True
            # get into this while loop as long as a new randomly picked Hi is SC;
            # its Ri is fertilized by any of the rest of its SLFs
            # if so create a new haplotype.
            while SC:
                tmp_RNase = createAlleles()  # Ri RNase
                # allele is an array of 18 AA indices
                current_RNase_list = []
                for i in range(seq_length): current_RNase_list.append(tmp_RNase[0].Allele[i].Letter.value)
                current_RNase_ar = np.array(current_RNase_list)
                interaction_vals = np.zeros(Model.slf_amount_in_haploid)
                for jj in range(Model.slf_amount_in_haploid):
                    tmp_SLF = createAlleles()  # Fj SLF
                    # allele is an array of 18 AA indices
                    current_SLF_list = []
                    for i in range(seq_length): current_SLF_list.append(tmp_SLF[0].Allele[i].Letter.value)
                    current_SLF_ar = np.array(current_SLF_list)
                    E_pool_Ri_Fj = Calc_Interaction_Energy(current_RNase_ar, current_SLF_ar, Model.eij, Model.seq_len)
                    interaction_vals[jj] = E_pool_Ri_Fj
                is_fertilized = interaction_vals < Model.interaction_thresh  # True- fertilized by any of the SLF already created, False- cannot be fertilized
                SC = is_fertilized.any()

            tmp_RNase[0].IndexInRNasesPool = n
            RNase_object.extend(deepcopy(tmp_RNase))
            all_RNase_list.append(current_RNase_list)
            all_RNase_ar = np.array(all_RNase_list)

    else:
        for n in range(Model.num_types_of_haplotypes):
            tmp_RNase = createAlleles()# Ri RNase
            # allele is an array of 18 AA indices
            current_RNase_list = []
            for i in range(seq_length): current_RNase_list.append(tmp_RNase[0].Allele[i].Letter.value)

            tmp_RNase[0].IndexInRNasesPool = n
            RNase_object.extend(deepcopy(tmp_RNase))
            all_RNase_list.append(current_RNase_list)
            all_RNase_ar = np.array(all_RNase_list)

    if Model.one_on_one: #
        for n in range(Model.num_types_of_haplotypes):
            not_fertilizes_Ri = True
            fertilizes_others = True

            # get into this while loop as long as a SLF Fi which was created in 'createAlleles' does not fertilize Ri
            # or fertilizes any of the rest R0 - R(i-1) RNases
            # and randomly pick a new Fi
            slf_ind = 0
            while not_fertilizes_Ri or fertilizes_others:
                tmp_SLF = createAlleles()  # Fj SLF
                # allele is an array of 18 AA indices
                current_SLF_list = []
                for i in range(seq_length): current_SLF_list.append(tmp_SLF[0].Allele[i].Letter.value)
                current_SLF_ar = np.array(current_SLF_list)

                E_pool_Ri_Fj = Calc_Interaction_Energy(all_RNase_ar[n], current_SLF_ar, Model.eij, Model.seq_len)
                interaction_strength = E_pool_Ri_Fj
                not_fertilizes_Ri = interaction_strength >= Model.interaction_thresh #True- not fertilizes, False- ferilizes

                # make sure that this SLF does not fertilize any of the previous RNases
                all_RNase_ar_exclude_n = np.delete(all_RNase_ar, (n), axis=0)
                E_pool_Ri_Fj = Calc_Interaction_Energy(all_RNase_ar_exclude_n, current_SLF_ar, Model.eij, Model.seq_len)
                interaction_strength = E_pool_Ri_Fj

                fertilizes_others = interaction_strength < Model.interaction_thresh #True-    False-
                fertilizes_others = fertilizes_others.any()
                slf_ind = slf_ind + 1

            tmp_SLF[0].IndexInSLFsPool = n
            SLF_object.extend(deepcopy(tmp_SLF))
            all_SLFs_list.append(current_SLF_list)

        haplotypes_dict = build_haplotypes(RNase_object, SLF_object)
        for jj in range(Model.added_slfs):
            fertilize_any_R = True
            while fertilize_any_R:
                tmp_SLF = createAlleles()
                # allele is an array of 18 AA indices
                current_SLF_list = []
                for i in range(seq_length): current_SLF_list.append(tmp_SLF[0].Allele[i].Letter.value)
                current_SLF_ar = np.array(current_SLF_list)
                E_pool_Ri_Fj = Calc_Interaction_Energy(all_RNase_ar, current_SLF_ar, Model.eij, Model.seq_len)
                interaction_strength = E_pool_Ri_Fj
                fertilize_any_R = interaction_strength < Model.interaction_thresh  # True-    False-
                fertilize_any_R = fertilize_any_R.any()
            tmp_SLF[0].IndexInSLFsPool = n+jj+1
            for hap_ind, hap in haplotypes_dict.items():
                hap.Male_part.append(deepcopy(tmp_SLF[0]))
            current_SLF_list = []
            for i in range(seq_length): current_SLF_list.append(tmp_SLF[0].Allele[i].Letter.value)
            all_SLFs_list.append(current_SLF_list)

        all_SLFs_ar = np.array(all_SLFs_list)
        for hap_ind, hap in haplotypes_dict.items():
            flat_list = [item for sublist in
                            [hap.SLFsAncestralIndices, np.arange(Model.num_types_of_haplotypes,
                            Model.num_types_of_haplotypes + Model.added_slfs, 1)] for item in sublist]
            hap.SLFsAncestralIndices = flat_list
            hap.SLFsIndices = flat_list
            hap.SLFsIndices_unique = np.array(flat_list)
            hap.SLFs_len = len(flat_list)
            hap.SLFs_unique_len = len(flat_list)

    else:
        slfs_new_inds_dict = {}
        num_slfs = 1000
        hap_dict = {}
        for n in range(num_slfs):
            tmp_SLF = createAlleles()
            tmp_SLF[0].IndexInSLFsPool = n
            SLF_object.extend(deepcopy(tmp_SLF))

        # build haplotypes
        current_SLFs_list = []

        for ii_rnase in range(Model.num_types_of_haplotypes):
            hap_dict[ii_rnase] = []

        RNases_fertilized_bool_arr = np.zeros((Model.num_types_of_haplotypes), dtype = bool)
        RNase_fertilized_bool_arr = np.zeros((Model.num_types_of_haplotypes-1), dtype=bool)
        fully_compatible = False
        all_SLFs_list = []
        while not fully_compatible:
            for jj_slf, tmp_SLF in enumerate(SLF_object): # iterate over all the SLFs
                current_SLF_list = []
                for i in range(seq_length): current_SLF_list.append(tmp_SLF.Allele[i].Letter.value)
                current_SLF_ar = np.array(current_SLF_list)
                all_SLFs_list.append(current_SLF_ar)
                if tmp_SLF.IndexInSLFsPool not in slfs_new_inds_dict:
                    slfs_new_inds_dict[tmp_SLF.IndexInSLFsPool] = []
                slfs_new_inds_dict[tmp_SLF.IndexInSLFsPool].append(len(current_SLFs_list)-1)
                for ii_rnase, tmp_RNase in enumerate(RNase_object):  # iterate over all the RNases
                    current_RNase_list = []
                    for i in range(seq_length): current_RNase_list.append(tmp_RNase.Allele[i].Letter.value)
                    current_RNase_ar = np.array(current_RNase_list)
                    E_pool_Ri_Fj = Calc_Interaction_Energy(current_RNase_ar, current_SLF_ar, Model.eij, Model.seq_len)
                    fertilizes_Ri = E_pool_Ri_Fj[0] < Model.interaction_thresh  # True- not fertilizes, False- ferilizes
                    if not fertilizes_Ri: #this SLF does not fertilizes this RNase
                        hap_dict[ii_rnase].append(deepcopy(tmp_SLF))#add

                for ii_rnase, tmp_RNase in enumerate(RNase_object):
                    if not RNases_fertilized_bool_arr[ii_rnase]:
                        current_RNase_list = []
                        for i in range(seq_length): current_RNase_list.append(tmp_RNase.Allele[i].Letter.value)
                        current_RNase_ar = np.array(current_RNase_list)
                        ii_hap_tmp = 0
                        for hap_ind, tmp_SLFs in hap_dict.items():
                            if hap_ind != ii_rnase:
                                current_SLFs_ar = np.zeros((len(tmp_SLFs), Model.seq_len))

                                for jj_slf, slf in enumerate(tmp_SLFs):
                                    current_SLF_list = []
                                    for i in range(seq_length): current_SLF_list.append(slf.Allele[i].Letter.value)
                                    current_SLF_ar = np.array(current_SLF_list)
                                    current_SLFs_ar[jj_slf,:] = current_SLF_ar
                                E_pool_Ri_Fj = Calc_Interaction_Energy(current_RNase_ar, current_SLFs_ar, Model.eij,
                                                           Model.seq_len)
                                RNase_fertilized_bool_arr[ii_hap_tmp] = any((E_pool_Ri_Fj < Model.interaction_thresh)[0])
                                ii_hap_tmp += 1
                        if all(RNase_fertilized_bool_arr): #check if the RNase is fertilized by all the haplotypes
                            RNases_fertilized_bool_arr[ii_rnase] = True
                if all(RNases_fertilized_bool_arr):
                    fully_compatible = True
                    break

        # haplotypes = []
        haplotypes_dict = {}
        for hap_ii in range(Model.num_types_of_haplotypes):
            haplotype = Model.Haplotype(RNase_object[hap_ii], hap_dict[hap_ii],
                        hap_ii, hap_ii, [], [], 0, 0, 'Initiation', 'False',
                                        [], [], [], [], [], [])

            # haplotype.Haplotype_index_in_pool = hap_ii
            haplotype.RNasesIndices = []
            haplotype.RNasesIndices.append(haplotype.Female_part.IndexInRNasesPool)
            haplotype.SLFsIndices = []
            for slf_ii, slf in enumerate(hap_dict[hap_ii]):
                slf_ind = slf.IndexInSLFsPool
                haplotype.SLFsIndices.append(slf_ind)
            haplotype.PreviousHaplotype = []
            haplotype.IS_SC = 'False'
            haplotype.SLFs_len = len(haplotype.SLFsIndices)
            haplotype.RNases_len = len(haplotype.RNasesIndices)
            haplotype.SLFsIndices_unique = np.unique(haplotype.SLFsIndices)
            haplotype.SLFsAncestralIndices = haplotype.SLFsIndices
            haplotype.RNasesAncestralIndices = haplotype.RNasesIndices
            haplotypes_dict[hap_ii] = haplotype

        all_SLFs_ar = np.array(all_SLFs_list)

    return haplotypes_dict, all_RNase_ar, all_SLFs_ar

def build_haplotypes(RNase_object, SLF_object):

    #build an Haplotype object that contain a female part Ri and a male part that contain all SLFs but Fi (the one that fertilizes Ri).
    # e.g. R0 F1 F2 F3
    #      R1 F0 F2 F3 etc.

    haplotypes_dict = {}
    num_of_alleles = Model.slf_amount_in_haploid

    for ind in range(Model.num_types_of_haplotypes):
        allele_indices_list = list(range(Model.num_types_of_haplotypes))
        allele_indices_list.remove(ind)  # all SLFs but the i'th, the one that fertilizes the i'th RNase

        alleles = []
        for k in range(num_of_alleles):
            tmp_all = deepcopy(SLF_object[allele_indices_list[k]])
            alleles.append(tmp_all)

        haplotype = Model.Haplotype(RNase_object[ind], alleles,
                    ind, ind, [], [], 0, 0, 'Initiation', 'False',
                                [], [], [], [], [], [])

        haplotype.Haplotype_index_in_pool = ind
        haplotype.RNasesIndices = []
        haplotype.RNasesIndices.append(haplotype.Female_part.IndexInRNasesPool)
        haplotype.SLFsIndices = []
        for k in range(Model.slf_amount_in_haploid):
            haplotype.SLFsIndices.append(haplotype.Male_part[k].IndexInSLFsPool)
        haplotype.PreviousHaplotype = []
        haplotype.BirthIter = 0
        haplotype.BirthReason = 'Initiation'
        haplotype.RNases_len = len(haplotype.RNasesIndices)
        haplotype.SLFs_unique_len = len(np.unique(haplotype.SLFsIndices))
        haplotype.SLFs_len = len(haplotype.SLFsIndices)
        haplotype.SLFsIndices_unique = np.unique(haplotype.SLFsIndices)
        haplotype.RNasesAncestralIndices = haplotype.RNasesIndices
        haplotype.SLFsAncestralIndices = haplotype.SLFsIndices
        haplotypes_dict[ind] = haplotype

    return haplotypes_dict


def is_compatible(E_Ri_Fj_dict, unique_SLFs_indices, unique_RNases_indices):
    # at least one of the SLFs can detoxify the RNase
    is_compatible = False
    for ii_slfs, slf_ind in enumerate(unique_SLFs_indices):
        if E_Ri_Fj_dict[int(unique_RNases_indices),slf_ind] < Model.interaction_thresh:
            is_compatible = True
            break

    return is_compatible

# -----------------------------------------------------------
#         CALCULATE INTERACTIONS
# -----------------------------------------------------------

def Calc_Interaction_Energy(female_al_pool, male_al_pool, eij, seq_len):

    '''
    this function calculates the interaction intensities between all the pairs of alleles in the pool.
    the Interaction between the two alleles is calculated by summing the interaction correspondences along the pair of alleles (i and j)
    The total interaction energy - equals the sum of interactions between the corresponding amino acids of the two protein surfaces

    input  - the Alleles pool of RNnases and SLFs
    output - E_pool_Ri_Fj an array of Interaction Energies between two alleles

    An Allele vector (RNase (R-allele)/ SLF (F-allele)) looks like this:
    Fj = {Fj1, Fj2, ... , FjL}
    Ri = {Ri1, Ri2, ... , RiL}
    L - is the length of the allele- a sequence of amino acids

    note - this function is called by 'initiate_haplotypes' only.
    Further needs of interaction calculations are calculated and saved as re-updating interactions dict
    '''

    #make sure that array shape is (num of alleles, seq_len)
    if len(male_al_pool) == 0:
        if female_al_pool.shape[0] == female_al_pool.size:
            female_al_pool = female_al_pool.reshape(1, seq_len)

    elif len(female_al_pool) == 0:
        if male_al_pool.shape[0] == male_al_pool.size:
            male_al_pool = male_al_pool.reshape(1, seq_len)

    else:
        if female_al_pool.shape[0] == female_al_pool.size:
            female_al_pool = female_al_pool.reshape(1, seq_len)

        if male_al_pool.shape[0] == male_al_pool.size:
            male_al_pool = male_al_pool.reshape(1, seq_len)

    num_female_types = int(np.size(female_al_pool)/seq_len)  # the amount of different alleles in the pool
    num_male_types = int(np.size(male_al_pool)/seq_len)  # the amount of different alleles in the pool
    E_pool_Ri_Fj = np.zeros((num_female_types, num_male_types))

    for ii in range(num_female_types):
        Ri = female_al_pool[ii]
        for jj in range(num_male_types):  # (ii, num_male_types):
            Fj = male_al_pool[jj]
            if any(np.isnan(male_al_pool[jj])):
                sum_eRi_Fj = np.empty((1))
                sum_eRi_Fj[:] = np.NaN
                E_pool_Ri_Fj[ii][jj] = sum_eRi_Fj[0]

            else:
                Fj = Fj.astype(int)
                sum_eRi_Fj = np.sum(eij[Ri, Fj], axis=0)

                E_pool_Ri_Fj[ii][jj] = sum_eRi_Fj

    return E_pool_Ri_Fj


def calc_female_fertilization_weights(diploid_count_dict, haplotypes_count_dict, delta_g_r_dict):
    # calculate fertilization weights for all the diploids in the population
    # based on the frequencies of the haplotypes that are able to fertilize them.
    # the purpose is to weight differently females that are able to be fertilized by whole population
    # and those that are able to be fretilized by a subset of the population.
    # In earlier vertion, females that are fertilized by a subset of the population and female that are fertilized by
    # the whole populaion have the same influence on the offspring set,
    # thus causing the population of non fertilizing haplotypes to be wiped out
    # and thus the population diversity becomes smaller

    diploid_fertilization_weights_dict = {}
    diploids_fertilized_by_haplotype_dict = {}

    for dip in diploid_count_dict.keys():

        for hap in haplotypes_count_dict.keys():

            g_r = tuple([item for sublist in [dip, [hap]] for item in sublist])
            is_C = delta_g_r_dict[g_r]
            if is_C:  # is_fertilized:
                if dip not in diploids_fertilized_by_haplotype_dict:  # the dict does not have this key yet
                    diploids_fertilized_by_haplotype_dict[dip] = []
                if hap not in diploids_fertilized_by_haplotype_dict[dip]:
                    diploids_fertilized_by_haplotype_dict[dip].append(hap)

    for dip in diploid_count_dict.keys():
        diploid_fertilization_weights_dict[dip] = 0
    for dip, hap_lst in diploids_fertilized_by_haplotype_dict.items():
        tmp_sum = 0
        for hap_ind in hap_lst:
            tmp_sum += haplotypes_count_dict[hap_ind]
        p_fertilize_haps =  tmp_sum/2/Model.popu_size
        Wi = 1 - (1-p_fertilize_haps)**Model.num_of_fertilizations_trials
        diploid_fertilization_weights_dict[dip] = Wi

    return diploid_fertilization_weights_dict

def calc_Rg_r_3D(offspring_dict, diploid_count_dict, delta_g_g_dict, cond1_g_g_dict, cond2_g_g_dict):

    # calculate R
    # given a genotype g and a foreign pollen with haploid genotype h / or a self pollen with haploid genotype g1 or g2
    # calculate the probability to produce an offspring with diploid genotype r

    R_g_g_r_dict = {}
    R_g_h_r_dict = {}

    # self fertilization
    for diploid_ind, diploid_i in enumerate(diploid_count_dict):
        diploid_i_lst = list(diploid_i)
        if delta_g_g_dict[diploid_i] == 0.5: # g is half self compatible
            is_C1 = cond1_g_g_dict[diploid_i] == 1
            is_C2 = cond2_g_g_dict[diploid_i] == 1
            if is_C1:
                offspring = [diploid_i[0], diploid_i[0]]
                g_h_r = tuple([item for sublist in [diploid_i_lst, [diploid_i[0]], offspring] for item in sublist])
                R_g_g_r_dict[g_h_r] = 0.5
                offspring = diploid_i_lst
                g_h_r = tuple([item for sublist in [diploid_i_lst, [diploid_i[0]], offspring] for item in sublist])
                R_g_g_r_dict[g_h_r] = 0.5
            if is_C2:
                offspring = [diploid_i[1], diploid_i[1]]
                g_h_r = tuple([item for sublist in [diploid_i_lst, [diploid_i[1]], offspring] for item in sublist])
                R_g_g_r_dict[g_h_r] = 0.5
                offspring = diploid_i_lst
                g_h_r = tuple([item for sublist in [diploid_i_lst, [diploid_i[1]], offspring] for item in sublist])
                R_g_g_r_dict[g_h_r] = 0.5
        if delta_g_g_dict[diploid_i] == 1: # g is fully self compatible
            if diploid_i[0] == diploid_i[1]: # g is homozygous,
                offspring = diploid_i_lst #the offspring and the parent are identicle
                g_h_r = tuple([item for sublist in [diploid_i_lst, [diploid_i[0]], offspring] for item in sublist])
                R_g_g_r_dict[g_h_r] = 1
            else: # g is heterozygous
                offspring = diploid_i_lst # g is heterozygous the offspring and the parent are identicle
                g_h_r = tuple([item for sublist in [diploid_i_lst, [diploid_i[0]], offspring] for item in sublist])
                R_g_g_r_dict[g_h_r] = 0.5
                g_h_r = tuple([item for sublist in [diploid_i_lst, [diploid_i[1]], offspring] for item in sublist])
                R_g_g_r_dict[g_h_r] = 0.5
                offspring = [diploid_i[0], diploid_i[0]] # g is heterozygous, the offspring is homozygous
                g_h_r = tuple([item for sublist in [diploid_i_lst, [diploid_i[0]], offspring] for item in sublist])
                R_g_g_r_dict[g_h_r] = 0.25
                offspring = [diploid_i[1], diploid_i[1]]
                g_h_r = tuple([item for sublist in [diploid_i_lst, [diploid_i[1]], offspring] for item in sublist])
                R_g_g_r_dict[g_h_r] = 0.25

    # foreign fertilizations
    for offspring_key, parents_values in offspring_dict.items():
        for parent in parents_values:
            haplotype_ind = parent[1]
            diploid_i = parent[0]

            g_h_r = tuple([item for sublist in [diploid_i, [haplotype_ind], list(offspring_key)] for item in sublist])
            if diploid_i[0] == diploid_i[1]:  # homozygosity
                if diploid_i == offspring_key:  # offspring equals diploid # i.e. g = S0S0, h = S0, r = S0S0
                    if diploid_i[0] == haplotype_ind:
                        R_g_h_r_dict[g_h_r] = 1

                else:  # the offspring r and the genotype g are not the same. # i.e. g = S0S0, h = S1, r = S0S1
                    if tuple([diploid_i[0], haplotype_ind]) == offspring_key or tuple(
                            [haplotype_ind, diploid_i[0]]) == offspring_key:
                        R_g_h_r_dict[g_h_r] = 1

            else:  # heterozygosity
                # the offspring r and the genotype g are not the same. # i.e. g = S0S1, h = S2, r = S0S2
                if tuple([diploid_i[0], haplotype_ind]) == offspring_key \
                        or tuple([haplotype_ind, diploid_i[0]]) == offspring_key \
                        or tuple([diploid_i[1], haplotype_ind]) == offspring_key \
                        or tuple([haplotype_ind, diploid_i[1]]) == offspring_key:
                    R_g_h_r_dict[g_h_r] = 0.5

    return R_g_g_r_dict, R_g_h_r_dict

def calc_pai_genotype_self_and_out(E_Ri_Fj_dict, haplotypes_dict, diploid_count_dict, haplotypes_count_dict,
                       delta_g_g_dict, delta_g_r_dict, cond1_g_g_dict, cond2_g_g_dict):

    '''
    input:
    E_Ri_Fj_dict - a dict with all alleles (RNases and SLFs in pool) interaction energies
    haplotypes_dict - a dict of all unique haplotypes' object
    diploid_count_dict - a dict of all unique diploid amounts
    haplotypes_count_dict - a dict of all unique haplotype amounts
    delta_g_g_dict -
    delta_g_r_dict -
    cond1_g_g_dict -
    cond2_g_g_dict -


    output:

    delta_g_g_dict, delta_g_r_dict, pai_g_g_self, pai_g_h_out, cond1_g_g_dict, cond2_g_g_dict


    pai_g_g_self - is the probability that an individual with a genotype g  is self-fertilized
                    ordered by the order of the diploids in "diploids_obj"
                   size array(number of unique genotypes X 1)

    pai_g_h_out - the probability that an individual with genotype g is outcrossed with haploid pollen h
                    ordered by the order of the diploids in "diploids_obj" and the haplotypes in "haplotypes_obj"
                  size array (number of unique genotypes X 2*(number of unique genotypes))

    assistance calculations:
    pr_multiply_delta_g_r - at this version if a haploid r in haplotypes_obj is one of the haplotypes in the
    fertilized genotype and it is its only time it appears. this haplotype r is summed in the summation over pr_multiply_delta_g_r
    '''

    len_genotype_collect = len(diploid_count_dict)
    len_haplotype_collect = len(haplotypes_count_dict)

    pai_g_g_self = np.zeros([len_genotype_collect, 1])
    pai_g_h_out = np.zeros([len_genotype_collect, len_haplotype_collect])

    for dip_count_ind, dip_keys in enumerate(diploid_count_dict):
        pr_multiply_delta_g_r = np.zeros([len_haplotype_collect])
        R_g1 = haplotypes_dict[dip_keys[0]].RNasesIndices
        R_g2 = haplotypes_dict[dip_keys[1]].RNasesIndices
        unique_RNases = np.unique([R_g1[0], R_g2[0]])
        if not dip_keys in delta_g_g_dict:

            # ----------------------------------------------
            #  check if the genotype is self compatible
            # ----------------------------------------------
            unique_SLFs1 = haplotypes_dict[dip_keys[0]].SLFsIndices_unique
            unique_SLFs2 = haplotypes_dict[dip_keys[1]].SLFsIndices_unique

            fertilization_cond1 = is_compatible(E_Ri_Fj_dict, unique_SLFs1, unique_RNases)
            fertilization_cond2 = is_compatible(E_Ri_Fj_dict, unique_SLFs2, unique_RNases)

            if fertilization_cond1 and fertilization_cond2: #Full self compatible
                can_be_self_fertilized = True
                delta_g_g = 1
                delta_g_g_dict[dip_keys] = delta_g_g
            elif bool(fertilization_cond1 and not fertilization_cond2): #half self compatible
                can_be_self_fertilized = True #half self compatible
                delta_g_g = 0.5
                delta_g_g_dict[dip_keys] = delta_g_g
                cond1_g_g_dict[dip_keys] = 1
                cond2_g_g_dict[dip_keys] = 0
            elif bool(not fertilization_cond1 and fertilization_cond2):
                can_be_self_fertilized = True  # half self compatible
                delta_g_g = 0.5
                delta_g_g_dict[dip_keys] = delta_g_g
                cond1_g_g_dict[dip_keys] = 0
                cond2_g_g_dict[dip_keys] = 1
            else: #self incompatible
                can_be_self_fertilized = False
                delta_g_g = 0
                delta_g_g_dict[dip_keys] = delta_g_g
        else:
            if delta_g_g_dict[dip_keys] == 0:
                can_be_self_fertilized = False
        # ----------------------------------------------
        #  check foreign fertilization
        # ----------------------------------------------

        can_be_fertilized = []
        sum_haplotypes_count_dict = sum(haplotypes_count_dict.values())
        for hap_ind, hap_keys in enumerate(haplotypes_count_dict): #.items():  # loop over all the haplotypes
            g_r = tuple([item for sublist in [dip_keys, [hap_keys]] for item in sublist])
            if g_r not in delta_g_r_dict:
                #if g_r not in delta_g_r_dict:
                unique_SLFs = haplotypes_dict[hap_keys].SLFsIndices_unique
                fertilization_cond = is_compatible(E_Ri_Fj_dict, unique_SLFs, unique_RNases)

                if fertilization_cond == True: # this haplotype fertilizes the genotype
                    delta_g_r_dict[g_r] = 1
                    can_be_fertilized.append(1)
                else:
                    delta_g_r_dict[g_r] = 0
                    can_be_fertilized.append(0)
            else:
                can_be_fertilized.append(delta_g_r_dict[g_r])
            hap_vals = haplotypes_count_dict[hap_keys]
            pr = hap_vals / sum_haplotypes_count_dict
            pr_multiply_delta_g_r[hap_ind] = pr * delta_g_r_dict[g_r]

        # at least one of the non-self haplotypes is enough for enabeling the fertilization
        can_be_fertilized = any(can_be_fertilized)

        if not can_be_fertilized and not can_be_self_fertilized : #skip the rest of PAI calculations for this genotype
            continue
        pai_denominator = Model.alpha * delta_g_g_dict[dip_keys] + (1 - Model.alpha) * np.sum(pr_multiply_delta_g_r)
        # calculate the probability that an individual with a genotype g (the female part of a diploid) is self-fertilized:
        pai_g_g_self[dip_count_ind] = Model.alpha * delta_g_g_dict[dip_keys]  / pai_denominator

        # calculate the probability that an individual with genotype g is outcrossed with haploid pollen h:
        pai_g_h_out_numerator = ((1 - Model.alpha) * pr_multiply_delta_g_r)
        pai_g_h_out[dip_count_ind] = pd.Series(pai_g_h_out_numerator) / pai_denominator

    if Model.debug_mode == 1:
        sum_pai_g = np.zeros((len_genotype_collect))
        is_sum_pai_g_0_1 = np.empty((len_genotype_collect))
        for diploid_ind in range(len_genotype_collect):
            sum_pai_g[diploid_ind] = pai_g_g_self[diploid_ind] + sum(pai_g_h_out[diploid_ind])
            is_sum_pai_g_0_1[diploid_ind] = abs(sum_pai_g[diploid_ind] - 1) < 1e-10 or sum_pai_g[diploid_ind] == 0
        if not all(is_sum_pai_g_0_1 == 1):
            print("pai self and out are incorrect")

    return delta_g_g_dict, delta_g_r_dict, pai_g_g_self, pai_g_h_out, cond1_g_g_dict, cond2_g_g_dict

def calc_x_r(offspring_dict, diploid_count_dict, haplotypes_count_dict,
             pai_g_g_self, pai_g_h_out, R_g_g_r_dict, R_g_h_r_dict, diploid_fertilization_weights_dict):
    # this function calculates x_r - the frequencies array of all current genotypes in the next generation

    x_g = diploid_count_dict
    diploidType_tuple = tuple(diploid_count_dict)
    haplotypeType_tuple = tuple(haplotypes_count_dict)
    sum_R_g_h = {}
    sum_R_g_g = {}
    num_offspring = len(offspring_dict)
    x_r = np.zeros((num_offspring,1))
    x_r_self_sum = np.zeros((num_offspring, 1))
    x_r_foreign_sum = np.zeros((num_offspring, 1))

    #the SELF term
    for diploid_ind, diploid_i in enumerate(diploid_count_dict.keys()):
        if pai_g_g_self[diploid_ind] != 0:
            diploid_i_lst = list(diploid_i)
            pai_g_g_self_one_diploid = pai_g_g_self[diploid_ind]
            if isinstance(pai_g_g_self_one_diploid, (np.ndarray)):
                pai_g_g_self_one_diploid = pai_g_g_self_one_diploid[0]

            x_g_one_diploid = x_g[diploid_i] / Model.popu_size  # frequency of a certain genotype in population
            if diploid_i[0] == diploid_i[1]:  # homozygous
                offspring = diploid_i_lst
                g_h_r = tuple([item for sublist in [diploid_i_lst, [diploid_i[0]], offspring] for item in sublist])
                if g_h_r in R_g_g_r_dict:
                    offspring_ind = tuple(offspring_dict.keys()).index(tuple(offspring))
                    calc_tmp = (1 - Model.delta[0]) * x_g_one_diploid * pai_g_g_self_one_diploid * R_g_g_r_dict[g_h_r]
                    x_r_self_sum[offspring_ind] += calc_tmp

            else:
                offspring = [diploid_i[0], diploid_i[0]]
                g_h_r = tuple([item for sublist in [diploid_i_lst, [diploid_i[0]], offspring] for item in sublist])
                if g_h_r in R_g_g_r_dict:
                    offspring_ind = tuple(offspring_dict.keys()).index(tuple(offspring))
                    calc_tmp = (1 - Model.delta[0]) * x_g_one_diploid * pai_g_g_self_one_diploid * R_g_g_r_dict[g_h_r]
                    x_r_self_sum[offspring_ind] += calc_tmp

                offspring = diploid_i_lst
                g_h_r = tuple([item for sublist in [diploid_i_lst, [diploid_i[0]], offspring] for item in sublist])
                if g_h_r in R_g_g_r_dict:
                    offspring_ind = tuple(offspring_dict.keys()).index(tuple(offspring))
                    calc_tmp = (1 - Model.delta[0]) * x_g_one_diploid * pai_g_g_self_one_diploid * R_g_g_r_dict[g_h_r]
                    x_r_self_sum[offspring_ind] += calc_tmp

                offspring = [diploid_i[1], diploid_i[1]]
                g_h_r = tuple([item for sublist in [diploid_i_lst, [diploid_i[1]], offspring] for item in sublist])
                if g_h_r in R_g_g_r_dict:
                    offspring_ind = tuple(offspring_dict.keys()).index(tuple(offspring))
                    calc_tmp = (1 - Model.delta[0]) * x_g_one_diploid * pai_g_g_self_one_diploid * R_g_g_r_dict[g_h_r]
                    x_r_self_sum[offspring_ind] += calc_tmp

                offspring = diploid_i_lst
                g_h_r = tuple([item for sublist in [diploid_i_lst, [diploid_i[1]], offspring] for item in sublist])
                if g_h_r in R_g_g_r_dict:
                    offspring_ind = tuple(offspring_dict.keys()).index(tuple(offspring))
                    calc_tmp = (1 - Model.delta[0]) * x_g_one_diploid * pai_g_g_self_one_diploid * R_g_g_r_dict[g_h_r]
                    x_r_self_sum[offspring_ind] += calc_tmp

    for offspring_ind, (offspring_key, parents_values) in enumerate(offspring_dict.items()):
        for parent in parents_values:
            diploid_i = parent[0]
            haplotype_i = parent[1]
            x_g_one_diploid = x_g[diploid_i] / Model.popu_size  # frequency of a certain genotype in population
            g_h_r = tuple([item for sublist in [diploid_i, [haplotype_i], offspring_key] for item in sublist])
            g_h = tuple([item for sublist in [diploid_i, [haplotype_i]] for item in sublist])
            if g_h not in sum_R_g_h:
                sum_R_g_h[g_h] = 0
            sum_R_g_h[g_h] += R_g_h_r_dict[g_h_r]
            if diploid_i not in sum_R_g_g:
                sum_R_g_g[diploid_i] = 0

            diploid_ind = diploidType_tuple.index(diploid_i)
            haplotype_ind = haplotypeType_tuple.index(haplotype_i)

            W_diploid_i = diploid_fertilization_weights_dict[diploid_i]
            calc_tmp \
                = x_g_one_diploid * pai_g_h_out[diploid_ind, haplotype_ind] \
                  * R_g_h_r_dict[g_h_r]*W_diploid_i

            x_r_foreign_sum[offspring_ind] += calc_tmp
        x_r[offspring_ind] = x_r_self_sum[offspring_ind] + x_r_foreign_sum[offspring_ind]
    return x_r

def calc_freq_diploid_next_gen(E_Ri_Fj_dict, diploid_count_dict, haplotypes_count_dict, haplotypes_dict, \
               delta_g_g_dict, delta_g_r_dict, cond1_g_g_dict, cond2_g_g_dict):
    # this function calculates x_r - the frequencies array of all current genotypes in the next generation

    delta_g_g_dict, delta_g_r_dict, pai_g_g_self, pai_g_h_out, cond1_g_g_dict, cond2_g_g_dict\
        = calc_pai_genotype_self_and_out\
        (E_Ri_Fj_dict, haplotypes_dict, diploid_count_dict, haplotypes_count_dict, \
         delta_g_g_dict, delta_g_r_dict, cond1_g_g_dict, cond2_g_g_dict)
    x_g = diploid_count_dict

    # ------------------------------------------------------------------------------------------------------
    #      a list of weights for all females based on the population proportion that fertilzes each female
    # ------------------------------------------------------------------------------------------------------
    diploid_fertilization_weights_dict = calc_female_fertilization_weights\
        (diploid_count_dict, haplotypes_count_dict, delta_g_r_dict)
    #--------------------------------------------------------
    #    create a list of potential offspring
    #--------------------------------------------------------
    # create an array of indices of all potential offsprings
    # each one is named by two indices 1. Haplotype1 2. Haplotype2
    offspring_dict = create_offspring_dict (diploid_count_dict, haplotypes_count_dict)

    # ---------------------------------------------------------------------
    #          calculate W - the mean fitness in the population
    # ---------------------------------------------------------------------
    num_diploids = len(diploid_count_dict)
    num_haplotypes = len(haplotypes_count_dict)
    W_self_term = np.zeros((num_diploids, 1))
    W_foreign_term = np.zeros((num_diploids, num_haplotypes))

    diploidType_tuple = tuple(diploid_count_dict)
    haplotypeType_tuple = tuple(haplotypes_count_dict)
    for genotype_ind, diploid_i in enumerate(diploidType_tuple):

        x_g_one_diploid = x_g[diploid_i] / Model.popu_size # frequency of a certain genotype in population
        W_diploid_i = diploid_fertilization_weights_dict[diploid_i]
        W_self_term[genotype_ind] =  (1 - Model.delta[0]) *x_g_one_diploid * pai_g_g_self[genotype_ind]

        for haplotype_ind, haplotype_i in enumerate(haplotypeType_tuple):
            # the second term in W ( the double summation over the types of g and type of haplotypes)
            W_foreign_term[genotype_ind, haplotype_ind] = x_g_one_diploid * \
              pai_g_h_out[genotype_ind, haplotype_ind]*W_diploid_i

    W = np.sum(W_self_term) + np.sum(W_foreign_term)

    #----------------------------------------------------------------------
    #             calculate Rgg_r and R_gh_r arrays
    #----------------------------------------------------------------------
    R_g_g_r_dict, R_g_h_r_dict = calc_Rg_r_3D(offspring_dict,
                 diploid_count_dict, delta_g_g_dict, cond1_g_g_dict, cond2_g_g_dict)
    #----------------------------------------------------------------------------------
    #           calculate x_r - the frequency of genotype r in the next generation
    #----------------------------------------------------------------------------------

    x_r_arr = calc_x_r(offspring_dict, diploid_count_dict, haplotypes_count_dict,
           pai_g_g_self, pai_g_h_out, R_g_g_r_dict, R_g_h_r_dict, diploid_fertilization_weights_dict)
    x_r_arr = x_r_arr /W
    x_r_arr = x_r_arr / sum(x_r_arr)
    return x_r_arr, offspring_dict, delta_g_g_dict, delta_g_r_dict, cond1_g_g_dict, cond2_g_g_dict

def build_the_next_generation(offspring_freq_arr, offspring_dict):

    '''
    the function receives the current the dict of all potential offsprings and their probabilities
    and build a new population of N new diploids

    Input:
    offspring_freq_arr   - offspring frequencies (X_r, sum(X_r) = 1)
    offspring_dict       - dictionary, keys - offspring diploids, values - diploids and haplotypes as parents

    output:
    diploid_count_dict    - dictionary of number of diploids in (EXCLUDING EMPTY OFFSPRING DIPLOIDS!!!)
    haplotypes_count_dict - dictionary of number of haplotypes in (INCLUDING EMPTY OFFSPRING HAPLOTYPES!!!)
    '''

    # a dictionary of unique genotypes (diploids) each one contains two values:
    # the haplotype1 and haplotype2 indices
    diploid_count_dict: Dict[Tuple[Any, Any], Union[int, Any]] = {}
    # a dictionary of unique haplotypes that contains the haplotypes indices
    haplotypes_count_dict: Dict[Tuple[Any], Union[int, Any]] = {}

    offspring_freq_arr = np.reshape(offspring_freq_arr, offspring_freq_arr.size)
    len_offspring_dict = len(offspring_dict)
    rand_offsprings_indices = np.random.choice(len_offspring_dict, Model.popu_size, p=offspring_freq_arr)
    unique_offsprings_indices, count = np.unique(rand_offsprings_indices, axis=0,
                                                         return_counts=True)

    # amount of every one of the offspring in the next generation
    offspring_count = np.zeros((len_offspring_dict, ), dtype=int)

    for a,b in zip(unique_offsprings_indices, count):
        offspring_count[a] = b

    offspring_haps_arr = np.unique([item for sublist in list(offspring_dict.keys()) for item in sublist])
    for hap_ind in offspring_haps_arr:
        haplotypes_count_dict[hap_ind] = 0
    for offspring_ind, offspring in enumerate(offspring_dict): #loop over all the offspring
        #fill the dict diploid_count_dict with new offspring that have more than 0 copies.
        if offspring_count[offspring_ind] == 0: # skip offspring that do not exist
             continue
        diploid_count_dict[offspring] = offspring_count[offspring_ind]

        haplotypes_count_dict[offspring[0]] += offspring_count[offspring_ind]
        haplotypes_count_dict[offspring[1]] += offspring_count[offspring_ind]
    haplotypes_count_dict = {x: y for x, y in haplotypes_count_dict.items() if y != 0}

    return diploid_count_dict, haplotypes_count_dict

def create_offspring_dict (diploid_count_dict, haplotypes_count_dict):
    offspring_dict = {}
    d_keys_list = list(diploid_count_dict.keys())
    h_keys_list = list(haplotypes_count_dict.keys())
    for d in d_keys_list:
        for h in h_keys_list:
            if d[0] <= h:
                inds = tuple([d[0],h])
            else:
                inds = tuple([h,d[0]])
            #inds = tuple(np.sort((d[0], h), 0))
            if inds not in offspring_dict:
                offspring_dict[inds] = []
            offspring_dict[inds].append([d, h])

            if d[1] <= h:
                inds = tuple([d[1],h])
            else:
                inds = tuple([h,d[1]])
            #inds = tuple(np.sort((d[1], h), 0))
            if inds not in offspring_dict:
                offspring_dict[inds] = []
            offspring_dict[inds].append([d, h])

    return offspring_dict

def init_AA_interactions_to_dict(all_RNase_ar, all_SLF_ar):

    #calculate interactions as a dict, on the following cycles this dict will be updated and not re-calculated
    E_Ri_Fj_dict = {}
    AA_rnases_dict = {}
    AA_slfs_dict = {}
    for rnase_ind, AA_rnase in enumerate(all_RNase_ar):
        AA_rnases_dict[rnase_ind] = AA_rnase
    for slf_ind, AA_slf in enumerate(all_SLF_ar):
        if not any(np.isnan(AA_slf)):
            AA_slf = AA_slf.astype(int)
            AA_slfs_dict[slf_ind] = AA_slf
    for rnase_ind, AA_rnase in enumerate(all_RNase_ar):
        for slf_ind, AA_slf in enumerate(all_SLF_ar):
            if not any(np.isnan(AA_slf)):
                AA_slf = AA_slf.astype(int)
                E_Ri_Fj_dict[rnase_ind, slf_ind] = np.sum(Model.eij[AA_rnase, AA_slf])

    return E_Ri_Fj_dict, AA_rnases_dict, AA_slfs_dict

def Update_AA_Interaction_dict(E_Ri_Fj_dict, AA_rnases_dict, AA_slfs_dict, haplotypes_count_dict, haplotypes_dict):

    # given E_Ri_Fj_dict from previous cycles, add interactions for new RNase-SLF couples.
    E_Ri_Fj_update_dict = {}
    RNases_count_dict = {}
    for keys, values in haplotypes_count_dict.items():
        if values > 0:
            RNases_hap = haplotypes_dict[keys].RNasesIndices
            for ii_rnase in RNases_hap:
                if ii_rnase not in RNases_count_dict:
                    RNases_count_dict[ii_rnase] = 0
                RNases_count_dict[ii_rnase] += values

    SLFs_count_dict = {}
    for keys, values in haplotypes_count_dict.items():
        if values > 0:
            SLFs_hap = haplotypes_dict[keys].SLFsIndices

            for ii_slf in SLFs_hap:
                if ii_slf not in SLFs_count_dict:
                    SLFs_count_dict[ii_slf] = 0
                SLFs_count_dict[ii_slf] += values

    for rnase_ind in RNases_count_dict:
        rnase_val = AA_rnases_dict[rnase_ind]
        for slf_ind in SLFs_count_dict:
            slf_val = AA_slfs_dict[slf_ind]
            if tuple([rnase_ind, slf_ind]) not in E_Ri_Fj_dict:
                E_Ri_Fj_update_dict[rnase_ind, slf_ind] = np.sum(Model.eij[rnase_val, slf_val], axis=0)
            else:
                E_Ri_Fj_update_dict[rnase_ind, slf_ind] = E_Ri_Fj_dict[rnase_ind, slf_ind]

    return E_Ri_Fj_update_dict

def save_files(path, AA_slfs_dict, AA_rnases_dict, E_Ri_Fj_dict):

    with open(path + '/AA_slfs_dict.pkl', 'wb') as output:
        pickle.dump(AA_slfs_dict, output, pickle.HIGHEST_PROTOCOL)

    with open(path + '/AA_rnases_dict.pkl', 'wb') as output:
        pickle.dump(AA_rnases_dict, output, pickle.HIGHEST_PROTOCOL)

    # with open(path + '/haplotypes_dict.pkl', 'wb') as output:
    #     pickle.dump(haplotypes_dict, output, pickle.HIGHEST_PROTOCOL)

    # with open(path + '/iters_haplotypes_dict.pkl', 'wb') as output:
    #     pickle.dump(iters_haplotypes_dict, output, pickle.HIGHEST_PROTOCOL)
    #
    # with open(path + '/iters_haplotypes_count_dict.pkl', 'wb') as output:
    #     pickle.dump(iters_haplotypes_count_dict, output, pickle.HIGHEST_PROTOCOL)
    with open(path + '/E_Ri_Fj_dict.pkl', 'wb') as output:
        pickle.dump(E_Ri_Fj_dict, output, pickle.HIGHEST_PROTOCOL)

