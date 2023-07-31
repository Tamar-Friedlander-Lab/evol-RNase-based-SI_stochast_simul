import simulationVars as Model
import simulationLogic as Logic
import mutation_module as Mutation
import allele_duplication_module as Duplicate
import allele_deletion_module as Delete
import pickle
import os
import UNIPROT_calcs
import numpy as np

def run_simulation(iterationIndex):
    bioclass_uniprot_freqs = UNIPROT_calcs.calc_mutation_probs_bio_classes()
    os.chdir("..")
    os.chdir("upper_data_dir")
    directory = "data_dir"

    subdirectories = ['2022-11-12-15-57-45-11'] #specific run
    subdirectory = subdirectories[0] #[iterationIndex]
    path = directory + "/" + subdirectory
    dir_rep = 'rerun_' + str(iterationIndex)
    path2 = path + dir_rep #os.path.join(path, dir_rep)
    os.mkdir(path2)

    print(path2)
    parents_haps = {}
    more_than_one_original_ancestor = False
    last_iter = 92310 #103170 #104400 #101540
    rnase_ind = 21728

    with open(path + '/haplotypes_dict_' + str(last_iter) + '.pkl', 'rb') as input:
        haplotypes_dict = pickle.load(input)
    with open(path + '/AA_slfs_dict.pkl', 'rb') as input:
        AA_slfs_dict = pickle.load(input)
    with open(path + '/AA_rnases_dict.pkl', 'rb') as input:
        AA_rnases_dict = pickle.load(input)
    with open(path + '/E_Ri_Fj_dict.pkl', 'rb') as input:
            E_Ri_Fj_dict = pickle.load(input)
    with open(path + '/diploid_count_dict_' + str(last_iter) + '.pkl', 'rb') as input:
        diploid_count_dict = pickle.load(input)
    with open(path + '/haplotypes_count_dict_' + str(last_iter) + '.pkl', 'rb') as input:
        haplotypes_count_dict = pickle.load(input)

    last_iter = last_iter + 1
    iter_ind = last_iter
    freq_SC_haps = 0

    delta_g_g_dict = {}
    delta_g_r_dict = {}
    cond1_g_g_dict = {}
    cond2_g_g_dict = {}

    rnase_still_alive = True
    while more_than_one_original_ancestor and rnase_still_alive \
            and freq_SC_haps < Model.SC_freq_thresh:

        # -------------------------------------------------------
        #                    duplicate alleles
        # -------------------------------------------------------

        haplotypes_dict, diploid_count_dict, haplotypes_count_dict \
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

        E_Ri_Fj_dict = Logic.Update_AA_Interaction_dict(E_Ri_Fj_dict, AA_rnases_dict, AA_slfs_dict,
                                                        haplotypes_count_dict, haplotypes_dict)

        # -------------------------------------------------------
        #   calculate diploids frequencies in the next generation
        # -------------------------------------------------------

        offspring_freq_arr, offspring_dict, delta_g_g_dict, delta_g_r_dict, cond1_g_g_dict, cond2_g_g_dict \
            = Logic.calc_freq_diploid_next_gen \
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

        # --------------------------------------------------------------------------------
        #    update delta_g_g_dict, delta_g_r_dict, cond1_g_g_dict, cond2_g_g_dict
        # --------------------------------------------------------------------------------

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
        #      amount of SC haplotypes
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

        freq_SC_haps = num_SC_haps / Model.popu_size / 2

        # -------------------------------------------------------------------------
        #     count the number of original ancestor RNases left in the population
        # -------------------------------------------------------------------------
        ancestor_rnases_arr = np.zeros(len(haplotypes_dict), dtype=int)
        for hap_ii, [hap_ind, hap] in enumerate(haplotypes_dict.items()):
            ancestor_rnases_arr[hap_ii] = hap.RNasesAncestralIndices[0]
        if len(np.unique(ancestor_rnases_arr)) == 1:
            more_than_one_original_ancestor = False

        # ---------------------------------------------------
        #   write haplotype/Diploid freqs/iters dicts
        # ---------------------------------------------------

        if not more_than_one_original_ancestor:
            first_iter_of_one_ancestor = iter_ind
            d_iters = 1
            for current_hap in haplotypes_dict:
                parents_haps[current_hap] = haplotypes_dict[current_hap].PreviousHaplotype
            if iter_ind % d_iters == 0:
                with open(path2 + '/parents_haps.pkl', 'wb') as output:
                    pickle.dump(parents_haps, output, pickle.HIGHEST_PROTOCOL)

        if iter_ind % d_iters == 0:
            # print("start generation # %d" % iter_ind)
            Logic.save_files(path2, AA_slfs_dict, AA_rnases_dict, E_Ri_Fj_dict)

            with open(path2 + '/diploid_count_dict_' + str(iter_ind) + '.pkl', 'wb') as output:
                pickle.dump(diploid_count_dict, output, pickle.HIGHEST_PROTOCOL)

            with open(path2 + '/haplotypes_count_dict_' + str(iter_ind) + '.pkl', 'wb') as output:
                pickle.dump(haplotypes_count_dict, output, pickle.HIGHEST_PROTOCOL)

            with open(path2 + '/haplotypes_dict_' + str(iter_ind) + '.pkl', 'wb') as output:
                pickle.dump(haplotypes_dict, output, pickle.HIGHEST_PROTOCOL)

        RNasesInHaplotypesDict = {}
        for hap_ind, hap in haplotypes_dict.items():
            RNases_indices = hap.RNasesIndices

            for rnase_ii, rnase_ind in enumerate(RNases_indices):
                if rnase_ind not in RNasesInHaplotypesDict:
                    RNasesInHaplotypesDict[rnase_ind] = []
                RNasesInHaplotypesDict[rnase_ind].append(hap_ind)

        if rnase_ind in RNasesInHaplotypesDict:
            rnase_still_alive = True
        else:
            rnase_still_alive = False
        iter_ind += 1

    Logic.save_files(path2, AA_slfs_dict, AA_rnases_dict, E_Ri_Fj_dict)

    with open(path2 + '/diploid_count_dict_' + str(iter_ind) + '.pkl', 'wb') as output:
        pickle.dump(diploid_count_dict, output, pickle.HIGHEST_PROTOCOL)

    with open(path2 + '/haplotypes_count_dict_' + str(iter_ind) + '.pkl', 'wb') as output:
        pickle.dump(haplotypes_count_dict, output, pickle.HIGHEST_PROTOCOL)

    with open(path2 + '/haplotypes_dict_' + str(iter_ind) + '.pkl', 'wb') as output:
        pickle.dump(haplotypes_dict, output, pickle.HIGHEST_PROTOCOL)