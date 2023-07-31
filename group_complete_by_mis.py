import numpy as np
import pickle
import os
import calc_weighted_num_partners_of_allele as weighted_num_partners
import simulationVars as Model

interaction_thresh = -6
d_iters = 10
N = Model.popu_size*2
os.chdir("..")
os.chdir("FlowersEvolutionsSimulation_ver13_data")
directory = "output_data_10haplotypes_k5_18AAs_E_thresh-6_p10-4"  # output_data_10haplotypes_k5_" + str(L) + "AAs_E_thresh-6_p2_2_10-5_mut_2" #test2" #test" #

subdirectories = next(os.walk(directory))[1]

num_files = len(subdirectories)
for subdirectory in subdirectories:
    output_dir = directory + "/" + subdirectory
    analyze_the_file = False
    if (not (os.path.exists(output_dir + '/mis_groups_hap_all_iters_dict.pkl'))):
        unique_RNases_all_iters_dict = {}
        mis_groups_hap_all_iters_dict = {}
        mis_groups_rnase_all_iters_dict = {}
        with open(output_dir + '/AA_rnases_dict.pkl', 'rb') as input:
            AA_rnases_dict = pickle.load(input)
        with open(output_dir + '/AA_slfs_dict.pkl', 'rb') as input:
            AA_slfs_dict = pickle.load(input)
        with open(output_dir + '/E_Ri_Fj_dict.pkl', 'rb') as input:
            E_Ri_Fj_dict = pickle.load(input)
        E_Ri_Fj_last_dict = {}

        # analyze the data from the moment of only one ancestor
        # (the population is composed of one ancestor -  the data was saved every 10 iters
        # before that - the population was composed of more than one ancestor
        # - the data was saved every 500 iters)
        # -------------------------------------------
        files_lst = os.listdir(output_dir)
        iters_lst = []
        for file_str in files_lst:
            if 'diploid_count_dict_' in file_str:
                iters_lst.append(int(file_str[19:-4]))
        if len(iters_lst) > 0:
            iters_arr = np.array(iters_lst)
            last_iter = np.max(iters_arr)
            first_iter = np.min(iters_arr)
            iters_lst = []
            for iter_ind in range(last_iter, first_iter, -d_iters):
                if (os.path.exists(output_dir + '/diploid_count_dict_' + str(iter_ind) + '.pkl')):
                    iters_lst.append(iter_ind)
                else:
                    break
            if len (iters_lst) >1:
                analyze_the_file = True
            last_iter = np.max(iters_lst)
            first_iter = np.min(iters_lst)
            # -------------------------------------------
            if analyze_the_file:
                print(os.path.join(subdirectory))
                hap_lst = []
                rnase_lst = []

                RNasesInHaplotypes_all_iters_Dict = {}
                for ii_iter, iter_ind in enumerate(range(first_iter, last_iter, d_iters)):
                    if not (os.path.isfile(output_dir + '/haplotypes_count_dict_' + str(iter_ind) + '.pkl') \
                        and os.path.isfile(output_dir + '/haplotypes_dict_' + str(iter_ind) + '.pkl')):
                        break
                    with open(output_dir + '/haplotypes_count_dict_' + str(iter_ind) + '.pkl', 'rb') as input:
                        haplotypes_count_dict = pickle.load(input)

                    with open(output_dir + '/haplotypes_dict_' + str(iter_ind) + '.pkl', 'rb') as input:
                        haplotypes_dict = pickle.load(input)

                    #a dict of RNases and their haplotypes
                    RNasesInHaplotypesDict = {}
                    for hap_ind, hap in haplotypes_dict.items():
                        RNases_indices = hap.RNasesIndices
                        for rnase_ii, rnase_ind in enumerate(RNases_indices):
                            if rnase_ind not in RNasesInHaplotypesDict:
                                RNasesInHaplotypesDict[rnase_ind] = []
                            RNasesInHaplotypesDict[rnase_ind].append(hap_ind)

                    RNasesInHaplotypes_all_iters_Dict[iter_ind] = RNasesInHaplotypesDict

                    haps_lst = list(haplotypes_count_dict.keys())
                    order = np.argsort(list(haplotypes_count_dict.values()))[::-1] #np.argsort(rnases_count)[::-1]
                    haps_arr = np.array(haps_lst)  # /np.sum(list(dt_new_any_RNase_dict.values())))
                    haps_arr_sorted = haps_arr[order]
                    comp_set = []
                    comp_set_by_RNases = []
                    #add first haplotype to first group (still empty)
                    # first check if SI:
                    ii = 0
                    SC_cond = True
                    while SC_cond:
                        hap0 = haps_arr_sorted[ii]
                        SC_cond = False
                        rnase_ind = haplotypes_dict[hap0].RNasesIndices[0]
                        rnase_AA = AA_rnases_dict[rnase_ind]
                        mis_groups_hap_dict = {}
                        mis_groups_rnase_dict = {}
                        for slf_ind in haplotypes_dict[hap0].SLFsIndices:
                            slf_AA = AA_slfs_dict[slf_ind]

                            if tuple([rnase_ind, slf_ind]) not in E_Ri_Fj_dict:
                                E = np.sum(Model.eij[rnase_AA, slf_AA], axis=0)
                            else:
                                E = E_Ri_Fj_dict[rnase_ind, slf_ind]
                            if E < interaction_thresh:
                                # interacting
                                SC_cond = True
                                break
                        ii+=1
                    #if not SC_cond: # SI
                    mis_groups_rnase_dict[0] = []
                    mis_groups_rnase_dict[0].append(haplotypes_dict[hap0].RNasesIndices[0])
                    mis_groups_hap_dict[0] = []
                    mis_groups_hap_dict[0].append(hap0)
                    for hap_partner in haps_arr_sorted[ii:]:
                        #check if hap_partner is SI:
                        SC_cond = False
                        rnase_ind = haplotypes_dict[hap_partner].RNasesIndices[0]
                        rnase_AA = AA_rnases_dict[rnase_ind]
                        for slf_ind in haplotypes_dict[hap_partner].SLFsIndices:
                            slf_AA = AA_slfs_dict[slf_ind]
                            if tuple([rnase_ind, slf_ind]) not in E_Ri_Fj_dict:
                                E = np.sum(Model.eij[rnase_AA, slf_AA], axis=0)
                            else:
                                E = E_Ri_Fj_dict[rnase_ind, slf_ind]
                            if E < interaction_thresh:
                                # interacting
                                SC_cond = True
                                break
                        if not SC_cond: #only SI haplotypes can be labeled to one of the groups
                            hap_ii = 0
                            comp_check_dist = {}
                            comp_def = []
                            for group_ii, hap_lst in mis_groups_hap_dict.items():
                                comp_check_dist[group_ii] = np.zeros(len(hap_lst), dtype = int)
                                for hap_ii, hap in enumerate(hap_lst):
                                    comp_cond1 = False
                                    rnase_ind = haplotypes_dict[hap].RNasesIndices[0]
                                    rnase_AA = AA_rnases_dict[rnase_ind]
                                    for slf_ind in haplotypes_dict[hap_partner].SLFsIndices:
                                        slf_AA = AA_slfs_dict[slf_ind]
                                        if tuple([rnase_ind, slf_ind]) not in E_Ri_Fj_dict:
                                            E = np.sum(Model.eij[rnase_AA, slf_AA], axis=0)
                                        else:
                                            E = E_Ri_Fj_dict[rnase_ind, slf_ind]
                                        if E < interaction_thresh:
                                            # interacting
                                            comp_cond1 = True
                                            break
                                    comp_cond2 = False
                                    rnase_ind = haplotypes_dict[hap_partner].RNasesIndices[0]
                                    rnase_AA = AA_rnases_dict[rnase_ind]
                                    for slf_ind in haplotypes_dict[hap].SLFsIndices:
                                        slf_AA = AA_slfs_dict[slf_ind]
                                        if tuple([rnase_ind, slf_ind]) not in E_Ri_Fj_dict:
                                            E = np.sum(Model.eij[rnase_AA, slf_AA], axis=0)
                                        else:
                                            E = E_Ri_Fj_dict[rnase_ind, slf_ind]
                                        if E < interaction_thresh:
                                            # interacting
                                            comp_cond2 = True
                                            break
                                    if comp_cond1 and comp_cond2:
                                        # symetric compatibility
                                        comp_check_dist[group_ii][hap_ii] = 2
                                    elif not comp_cond1 and not comp_cond2:
                                        comp_check_dist[group_ii][hap_ii] = 0
                                    else:
                                        comp_check_dist[group_ii][hap_ii] = 1

                            for group_ii in comp_check_dist.keys():

                                if all(comp_check_dist[group_ii] == 2):
                                    comp_def.append(2) #comp
                                elif all (comp_check_dist[group_ii] == 0):
                                    comp_def.append(0) #non_comp
                                else:
                                    comp_def.append(1)  # partial

                            if all(np.array(comp_def) == 2):
                                ind = len(mis_groups_hap_dict)
                                mis_groups_hap_dict[ind] = []
                                mis_groups_hap_dict[ind].append(hap_partner)
                                mis_groups_rnase_dict[ind] = []
                                mis_groups_rnase_dict[ind].append(haplotypes_dict[hap_partner].RNasesIndices[0])
                            elif len(np.where(np.array(comp_def)==0)[0] == 1) and not any(np.array(comp_def) == 1):
                                ind = np.where(np.array(comp_def)==0)[0][0]
                                mis_groups_hap_dict[ind].append(hap_partner)
                                mis_groups_rnase_dict[ind].append(haplotypes_dict[hap_partner].RNasesIndices[0])

                    mis_groups_hap_all_iters_dict[iter_ind] = mis_groups_hap_dict
                    mis_groups_rnase_all_iters_dict[iter_ind] = mis_groups_rnase_dict

                with open(output_dir + '/mis_groups_hap_all_iters_dict.pkl', 'wb') as output:
                    pickle.dump(mis_groups_hap_all_iters_dict, output, pickle.HIGHEST_PROTOCOL)

                with open(output_dir + '/mis_groups_rnase_all_iters_dict.pkl', 'wb') as output:
                    pickle.dump(mis_groups_rnase_all_iters_dict, output, pickle.HIGHEST_PROTOCOL)

                with open(output_dir + '/RNasesInHaplotypes_all_iters_Dict.pkl', 'wb') as output:
                    pickle.dump(RNasesInHaplotypes_all_iters_Dict, output, pickle.HIGHEST_PROTOCOL)