import simulationVars as Model
import simulationLogic as Logic
import mutation_module as Mutation
import allele_duplication_module as Duplicate
import allele_deletion_module as Delete

import UNIPROT_calcs
import numpy as np
import pickle
from decimal import Decimal
import time
import os
from numpy import random
from typing import Dict, Any, Tuple, Union

def run_simulation(iterationIndex):

    #initiate bioclass probabilities
    bioclass_uniprot_freqs = UNIPROT_calcs.calc_mutation_probs_bio_classes()

    #define unique dir name
    year = time.strftime("%Y")
    month = time.strftime("%m")
    month_int = int(month)
    day = time.strftime("%d")
    day_int = int(day)
    hour = time.strftime("%H")
    hour_int = int(hour)
    minute = time.strftime("%M")
    minute_int = int(minute)
    second = time.strftime("%S")
    second_int = int(second)
    timestr = year + "-" + month + "-" + day + "-" + hour + "-" +  minute + "-" + second + "-" + str(iterationIndex)
    seed_strftime = int(month_int * 1e8 + day_int * 1e6 + hour_int * 1e5 + minute_int * 1e2 + second_int)
    if Model.use_seed == 1:
        random.seed(seed_strftime)

    directory = timestr
    print(directory)
    parent_dir = "output_data_10haplotypes"
    path = os.path.join(parent_dir, directory)
    os.mkdir(path)
    if Model.use_seed == 1:
        np.random.seed(None)
        st0 = np.random.get_state()
        with open(path + '/init_state_RNG.pkl', 'wb') as output:
            pickle.dump(st0, output, pickle.HIGHEST_PROTOCOL)

    log_file = open(path + '/log_file_'
                                + timestr
                                + '.txt', "a")

    #write to log file
    log_file.write("Directory: " +directory + "\n"
                   "population size: " + str(Model.popu_size) + "\n" 
                    "num of iteration: " + str(Model.num_of_iters) + "\n"
                    "allele length: " + str(Model.seq_len) + "\n"
                    "initialized num of haplotypes: " + str(Model.num_types_of_haplotypes) + "\n"
                    "initialized num of SLFs in haplotype: " + str(Model.slf_amount_in_haploid) + "\n"
                    "num of AA bio-classes: " + str(Model.alphabet_size) + "\n"
                    "alpha: " + str(Model.alpha) + "\n"
                    "delta: " + str(Model.delta[0]) + "\n"
                    "interaction threshold: " + str(Model.interaction_thresh) + "\n"
                    "p mutation: " +"{:.2E}".format(Decimal(str(Model.p_mutation))) + "\n"
                    "p duplication: " + "{:.2E}".format(Decimal(str(Model.p_duplicate))) + "\n"
                    "p deletion: " + "{:.2E}".format(Decimal(str(Model.p_delete))) + "\n"
                    "5'th class for STOP CODON? - NO \n"  
                    "initialization of bio-classes frequencies: " + Model.initiate_bioclasses_freq + "\n" 
                    "initialized as SC?" + Model.SC_type + "\n"
                    "k(number of female fertilization trials - new model)" + str(Model.num_of_fertilizations_trials) + "\n"   
                    "initialization the population - " + Model.initiation_RNases_costrained_to_log_file + "\n"                                                                                                          
                    "" + Model.initiation_text_to_log_file + "\n" 
                    "RNases duplication is " + Model.enable_dup + "\n"
                    "RNases deletion is " + Model.enable_del + "\n"
                    "seed file in init_state_RNG.pkl, for using seed type: with open('/init_state_RNG.pkl', 'rb') as input:) \n"
                    "and - np.random.set_state(st0) \n"
                    "see: https://stackoverflow.com/questions/32172054/how-can-i-retrieve-the-current-seed-of-numpys-random-number-generator")

    log_file.close()

    # a dict of hap - parent hap mapping
    parents_haps = {}
    # -----------------------------------------------------
    #      Inintiate the population
    #------------------------------------------------------

    haplotypes_init_dict, all_RNase_ar, all_SLF_ar = Logic.initiate_haplotypes()

    E_Ri_Fj_dict, AA_rnases_dict, AA_slfs_dict \
        = Logic.init_AA_interactions_to_dict(all_RNase_ar, all_SLF_ar)

    with open(path + '/initial_unique_haplotypes_dict.pkl', 'wb') as output:
        pickle.dump(haplotypes_init_dict, output, pickle.HIGHEST_PROTOCOL)

    num_of_init_diploids = int(Model.num_types_of_haplotypes * (Model.num_types_of_haplotypes - 1) / 2)

    # if union probability for all the diploids
    p = 1/num_of_init_diploids * np.ones(num_of_init_diploids)
    diploid_count_dict, haplotypes_count_dict, haplotypes_dict \
        = Logic.initiate_population_given_prob(haplotypes_init_dict, p)

    # save the initial population using pickle
    with open(path + '/initial_haplotypes_dict.pkl', 'wb') as output:
        pickle.dump(haplotypes_dict, output, pickle.HIGHEST_PROTOCOL)

    with open(path + '/initial_diploid_count_dict.pkl', 'wb') as output:
        pickle.dump(diploid_count_dict, output, pickle.HIGHEST_PROTOCOL)

    with open(path + '/initial_haplotypes_count_dict.pkl', 'wb') as output:
        pickle.dump(haplotypes_count_dict, output, pickle.HIGHEST_PROTOCOL)

    with open(path + '/AA_rnases_dict.pkl', 'wb') as output:
        pickle.dump(AA_rnases_dict, output, pickle.HIGHEST_PROTOCOL)
    with open(path + '/AA_slfs_dict.pkl', 'wb') as output:
        pickle.dump(AA_slfs_dict, output, pickle.HIGHEST_PROTOCOL)

    with open(path + '/E_Ri_Fj_dict.pkl', 'wb') as output:
        pickle.dump(E_Ri_Fj_dict, output, pickle.HIGHEST_PROTOCOL)

    #------------------------------------------------------------
    #          initializations
    # -----------------------------------------------------------
    print("start generation # 0")
    iter_ind = 0
    freq_SC_haps = 0

    #initiate selection module dicts
    delta_g_g_dict = {}
    delta_g_r_dict = {}
    cond1_g_g_dict = {}
    cond2_g_g_dict = {}

    more_than_one_original_ancestor= True
    d_iters = 500
    first_iter_of_one_ancestor = 0
    while [more_than_one_original_ancestor \
           or iter_ind < Model.num_of_iters + first_iter_of_one_ancestor]\
            and freq_SC_haps < Model.SC_freq_thresh:

        #-------------------------------------------------------
        #                    duplicate alleles
        #-------------------------------------------------------

        haplotypes_dict, diploid_count_dict, haplotypes_count_dict\
            = Duplicate.duplicate_population(iter_ind, haplotypes_dict, diploid_count_dict, haplotypes_count_dict)

        # -------------------------------------------------------
        #                    delete alleles
        # -------------------------------------------------------

        haplotypes_dict, diploid_count_dict, haplotypes_count_dict \
            = Delete.delete_population(iter_ind, haplotypes_dict, diploid_count_dict,
                                             haplotypes_count_dict)

        # -------------------------------------------------------
        #                    mutate population
        # -------------------------------------------------------
        haplotypes_dict, diploid_count_dict, haplotypes_count_dict, \
            AA_rnases_dict, AA_slfs_dict = \
            Mutation.mutate_population(iter_ind, bioclass_uniprot_freqs, haplotypes_dict,
           diploid_count_dict, haplotypes_count_dict, AA_rnases_dict, AA_slfs_dict)

        # -----------------------------------------------------------------
        #                  update interactions
        # -----------------------------------------------------------------
        E_Ri_Fj_dict = Logic.Update_AA_Interaction_dict\
            (E_Ri_Fj_dict, AA_rnases_dict, AA_slfs_dict, haplotypes_count_dict, haplotypes_dict)

        # -------------------------------------------------------
        #   calculate diploids frequencies in the next generation
        # -------------------------------------------------------
        offspring_freq_arr, offspring_dict, delta_g_g_dict, delta_g_r_dict, cond1_g_g_dict, cond2_g_g_dict \
            = Logic.calc_freq_diploid_next_gen\
            (E_Ri_Fj_dict, diploid_count_dict, \
            haplotypes_count_dict, haplotypes_dict, \
             delta_g_g_dict, delta_g_r_dict, cond1_g_g_dict, cond2_g_g_dict)

        # -------------------------------------------------------
        #                  build next generation
        # -------------------------------------------------------
        diploid_count_dict, haplotypes_count_dict = \
            Logic.build_the_next_generation(offspring_freq_arr, offspring_dict)

        # --------------------------------------------------------------------------------
        #    update haplotypes_dict
        # --------------------------------------------------------------------------------

        haps_lst = haplotypes_count_dict.keys()
        tmpdict = haplotypes_dict
        haplotypes_dict = dict((k, tmpdict[k]) for k in haps_lst if k in tmpdict)

        #--------------------------------------------------------------------------------
        #    update delta_g_g_dict, delta_g_r_dict, cond1_g_g_dict, cond2_g_g_dict
        #--------------------------------------------------------------------------------

        diploid_lst = diploid_count_dict.keys()
        tmpdict = delta_g_g_dict
        delta_g_g_dict = dict((k, tmpdict[k]) for k in diploid_lst if k in tmpdict)
        tmpdict = cond1_g_g_dict
        cond1_g_g_dict = dict((k, tmpdict[k]) for k in diploid_lst if k in tmpdict)
        tmpdict = cond2_g_g_dict
        cond2_g_g_dict = dict((k, tmpdict[k]) for k in diploid_lst if k in tmpdict)

        g_r = {}
        for dip_keys in diploid_count_dict.keys():
            for hap_keys in haplotypes_count_dict:
                g_r_key = tuple([item for sublist in [dip_keys, [hap_keys]] for item in sublist])
                g_r[g_r_key] = []

        lst = g_r.keys()
        tmpdict = delta_g_r_dict
        delta_g_r_dict = dict((k, tmpdict[k]) for k in lst if k in tmpdict)

        # ------------------------------------------------------
        #      count SC haplotypes
        # ------------------------------------------------------
        num_SC_haps = 0
        for hap_ind, hap in haplotypes_dict.items():
            rnase_ind = hap.RNasesIndices[0]
            rnase_AA = AA_rnases_dict[rnase_ind]
            comp_cond = False
            for slf_ind in hap.SLFsIndices:
                slf_AA = AA_slfs_dict[slf_ind]
                if tuple([rnase_ind, slf_ind]) not in E_Ri_Fj_dict:
                    E = np.sum(Model.eij[rnase_AA, slf_AA], axis=0)
                else:
                    E = E_Ri_Fj_dict[rnase_ind, slf_ind]
                if E < Model.interaction_thresh:
                    # interacting
                    comp_cond = True
                    break
            if comp_cond:
                hap.IS_SC = 'True'
                num_SC_haps += haplotypes_count_dict[hap.Haplotype_index_in_pool]

        freq_SC_haps = num_SC_haps/Model.popu_size/2

        #-------------------------------------------------------------------------
        #     count the number of original ancestor RNases left in the population
        #-------------------------------------------------------------------------
        ancestor_rnases_arr = np.zeros(len(haplotypes_dict), dtype = int)
        for hap_ii, [hap_ind, hap] in enumerate(haplotypes_dict.items()):
            ancestor_rnases_arr[hap_ii] = hap.RNasesAncestralIndices[0]
        if len(np.unique(ancestor_rnases_arr)) == 1:
            more_than_one_original_ancestor = False

        #------------------------------------------------------------------------
        #       save data
        #------------------------------------------------------------------------
        if not more_than_one_original_ancestor:
            first_iter_of_one_ancestor = iter_ind
            d_iters = 10
            for current_hap in haplotypes_dict:
                parents_haps[current_hap] = haplotypes_dict[current_hap].PreviousHaplotype
            if iter_ind % d_iters == 0:
                with open(path + '/parents_haps.pkl', 'wb') as output:
                    pickle.dump(parents_haps, output, pickle.HIGHEST_PROTOCOL)

        if iter_ind % d_iters == 0:
            print("start generation # %d" % iter_ind)
            Logic.save_files(path, AA_slfs_dict, AA_rnases_dict, E_Ri_Fj_dict)

            with open(path + '/diploid_count_dict_' + str(iter_ind) + '.pkl', 'wb') as output:
                pickle.dump(diploid_count_dict, output, pickle.HIGHEST_PROTOCOL)

            with open(path + '/haplotypes_count_dict_' + str(iter_ind) + '.pkl', 'wb') as output:
                pickle.dump(haplotypes_count_dict, output, pickle.HIGHEST_PROTOCOL)

            with open(path + '/haplotypes_dict_' + str(iter_ind) + '.pkl', 'wb') as output:
                pickle.dump(haplotypes_dict, output, pickle.HIGHEST_PROTOCOL)

        iter_ind += 1

    Logic.save_files(path, AA_slfs_dict, AA_rnases_dict, E_Ri_Fj_dict)

    with open(path + '/diploid_count_dict_' + str(iter_ind) + '.pkl', 'wb') as output:
        pickle.dump(diploid_count_dict, output, pickle.HIGHEST_PROTOCOL)

    with open(path + '/haplotypes_count_dict_' + str(iter_ind) + '.pkl', 'wb') as output:
        pickle.dump(haplotypes_count_dict, output, pickle.HIGHEST_PROTOCOL)

    with open(path + '/haplotypes_dict_' + str(iter_ind) + '.pkl', 'wb') as output:
        pickle.dump(haplotypes_dict, output, pickle.HIGHEST_PROTOCOL)
