import numpy as np
import pickle
import os
import matplotlib.pyplot as plt
import simulationVars as Model
min_duration = 10
d_iters = 10
thresh_num_hap_in_group = 10
interaction_thresh = -6
os.chdir("..")
os.chdir("FlowersEvolutionsSimulation_ver13_data")
directory = "output_data_10haplotypes_k5_18AAs_E_thresh-6_p10-4"
subdirectories = next(os.walk(directory))[1]

num_files = len(subdirectories)
perc_extinct_rnase_mut = np.zeros(len(subdirectories))
perc_extinct_rnase_mut[:] = np.nan
perc_extinct_mut_out_orig_in = np.zeros(len(subdirectories))
perc_extinct_mut_out_orig_in[:] = np.nan
perc_mut_out_orig_in_that_cause_extinct = np.zeros(len(subdirectories))
perc_mut_out_orig_in_that_cause_extinct[:] = np.nan
dt_mut_out_orig_in_extinct_dict = {}
for ii_dir, subdirectory in enumerate(subdirectories):
    output_dir = directory + "/" + subdirectory
    if not os.path.isfile(output_dir + '/analyze_extinct_group_all_iters_' + str(min_duration * 10) + '_iters.pkl'):
        if os.path.isfile(output_dir + '/extinct_group_all_iters_' + str(min_duration*10) + '_iters_dict.pkl') \
                and os.path.isfile(output_dir + '/mut_rnases_' + str(min_duration * 10) + '_iters_dict.pkl'):
            print(os.path.join(subdirectory))
            with open(output_dir + '/AA_rnases_dict.pkl', 'rb') as input:
                AA_rnases_dict = pickle.load(input)
            with open(output_dir + '/AA_slfs_dict.pkl', 'rb') as input:
                AA_slfs_dict = pickle.load(input)
            with open(output_dir + '/extinct_group_all_iters_' + str(min_duration*10) + '_iters_dict.pkl', 'rb') as input:
                extinct_group_all_iters_dict = pickle.load(input)
            with open(output_dir + '/reorder_mis_strict_groups_hap_all_iters_dict.pkl', 'rb') as input:
                reorder_mis_strict_groups_hap_all_iters_dict = pickle.load(input)
            with open(output_dir + '/reorder_mis_rnase_strict_groups_all_iters_dict.pkl', 'rb') as input:
                reorder_mis_rnase_strict_groups_all_iters_dict = pickle.load(input)
            with open(output_dir + '/reorder_mis_strict_group_extension_hap_all_iters_dict.pkl', 'rb') as input:
                reorder_mis_strict_group_extension_hap_all_iters_dict = pickle.load(input)
            with open(output_dir + '/mut_rnases_' + str(min_duration * 10) + '_iters_dict.pkl', 'rb') as input:
                mut_rnases_dict = pickle.load(input)

            analyze_extinct_group_all_iters_dict = {}
            small_ext_groups_all_iters_dict = {}
            iter_ext = np.array(list(extinct_group_all_iters_dict.keys())) #\
            mut_iters_lst = list(mut_rnases_dict.keys())
            mut_iters_lst.sort(reverse=True)
            mut_iters_arr = np.array(mut_iters_lst)
            for iter_ind in iter_ext: #loop over all extinction events
                gro_lst = extinct_group_all_iters_dict[iter_ind]
                for gro_ind in gro_lst: #loop over all the extinct groups in this generation
                    mut_iters_arr_before_ext = mut_iters_arr[
                        mut_iters_arr <= iter_ind]  # RNase mutations that occured before that extinction
                    mut_iters_arr_before_ext_copy = []
                    for mut_iter_ii, mut_iter_ind in enumerate(mut_iters_arr_before_ext):
                        mut_rnase_cond1 = False
                        mut_rnase_cond2 = False
                        if gro_ind in reorder_mis_strict_group_extension_hap_all_iters_dict[mut_iter_ind]:
                            mut_rnase_cond1 = True
                        if gro_ind in reorder_mis_strict_groups_hap_all_iters_dict[mut_iter_ind]:
                            mut_rnase_cond2 = True
                        if mut_rnase_cond1 or mut_rnase_cond2:
                            mut_iters_arr_before_ext_copy.append(mut_iter_ind)
                    mut_iters_arr_before_ext = mut_iters_arr_before_ext_copy
                    if len(mut_iters_arr_before_ext) == 1:
                        if iter_ind not in small_ext_groups_all_iters_dict:
                            small_ext_groups_all_iters_dict[iter_ind] = []
                        small_ext_groups_all_iters_dict[iter_ind].append(gro_ind)
                    if len(mut_iters_arr_before_ext) == 0:
                        if iter_ind not in small_ext_groups_all_iters_dict:
                            small_ext_groups_all_iters_dict[iter_ind] = []
                        small_ext_groups_all_iters_dict[iter_ind].append(gro_ind)
                        continue # non of the RNase mutations occured before that extinction
                    ii = 0
                    iter_first = 0
                    rnase_mut_iter_ii = mut_iters_arr_before_ext[ii]
                    ii_iter = rnase_mut_iter_ii
                    # start at the RNase mutation iteration, and at the extinction iteration
                    haps_in_gro = []
                    haps_in_ext = []
                    haps_in_gro_ext1 = []
                    if gro_ind in list(reorder_mis_strict_groups_hap_all_iters_dict[ii_iter].keys()):
                        haps_in_gro = reorder_mis_strict_groups_hap_all_iters_dict[ii_iter][gro_ind]
                        haps_in_gro_ext1 = haps_in_gro
                    haps_in_gro_ext2 = []
                    if len(haps_in_gro_ext1)>0:
                        if gro_ind in list(reorder_mis_strict_group_extension_hap_all_iters_dict[ii_iter].keys()):
                            haps_in_ext = reorder_mis_strict_group_extension_hap_all_iters_dict[ii_iter][gro_ind]
                            haps_in_gro_ext2 = haps_in_ext
                        x = range(len(mut_rnases_dict[rnase_mut_iter_ii]))
                        for ii_mut_rnase in range(len(mut_rnases_dict[rnase_mut_iter_ii])):
                            rnase_orig = mut_rnases_dict[rnase_mut_iter_ii][ii_mut_rnase]['orig RNase']
                            rnase_mut = mut_rnases_dict[rnase_mut_iter_ii][ii_mut_rnase]['mutant RNase']
                            rnase_gro_lst = []
                            if rnase_orig not in reorder_mis_rnase_strict_groups_all_iters_dict[ii_iter] and \
                                    rnase_mut not in reorder_mis_rnase_strict_groups_all_iters_dict[ii_iter]:
                                continue
                            with open(output_dir + '/haplotypes_dict_' + str(ii_iter) + '.pkl', 'rb') as input:
                                haplotypes_dict = pickle.load(input)
                            with open(output_dir + '/haplotypes_count_dict_' + str(ii_iter) + '.pkl', 'rb') as input:
                                haplotypes_count_dict = pickle.load(input)
                            #check if male haplotypes (in generation 'ii_iter')
                            # in the group that is about to extinct pollinate the orig RNase
                            rnase_AA = AA_rnases_dict[rnase_orig]
                            is_orig_comp = np.zeros(len(haps_in_gro_ext1), dtype=bool)
                            for hap_ii, hap_ind in enumerate(haps_in_gro_ext1):
                                slfs_ind = haplotypes_dict[hap_ind].SLFsIndices
                                slf_cond = False
                                for slf_ii, slf_ind in enumerate(slfs_ind):
                                    slf_AA = AA_slfs_dict[slf_ind]
                                    E = np.sum(Model.eij[rnase_AA, slf_AA], axis=0)
                                    if E < interaction_thresh:
                                        # interacting SLF
                                        slf_cond = True
                                        is_orig_comp[hap_ii] = True
                                        break
                            # check if male haplotypes in generation ('ii_iter')
                            # in the group that is about to extinct pollinate the mut RNase
                            rnase_AA = AA_rnases_dict[rnase_mut]
                            is_mut_comp = np.zeros(len(haps_in_gro_ext1), dtype=bool)
                            n_comp_mut = 0
                            n_non_comp_mut = 0
                            for hap_ii, hap_ind in enumerate(haps_in_gro_ext1):
                                slfs_ind = haplotypes_dict[hap_ind].SLFsIndices
                                slf_cond = False
                                for slf_ii, slf_ind in enumerate(slfs_ind):
                                    slf_AA = AA_slfs_dict[slf_ind]
                                    E = np.sum(Model.eij[rnase_AA, slf_AA], axis=0)
                                    if E < interaction_thresh:
                                        # interacting SLF
                                        slf_cond = True
                                        is_mut_comp[hap_ii] = True
                                        n_comp_mut += haplotypes_count_dict[hap_ind]
                                        break
                                if not slf_cond:
                                    n_non_comp_mut += haplotypes_count_dict[hap_ind]

                            if n_non_comp_mut/(n_comp_mut + n_non_comp_mut)>0.75 and all(is_orig_comp):
                                if rnase_orig in reorder_mis_rnase_strict_groups_all_iters_dict[ii_iter]:
                                    rnase_gro_lst.append(reorder_mis_rnase_strict_groups_all_iters_dict[ii_iter][rnase_orig])
                                if rnase_mut in reorder_mis_rnase_strict_groups_all_iters_dict[ii_iter]:
                                    rnase_gro_lst.append(reorder_mis_rnase_strict_groups_all_iters_dict[ii_iter][rnase_mut])
                                unique_gro, count = np.unique(rnase_gro_lst, return_counts=True)
                                rnase_gro = unique_gro[np.argmax(count)]

                                tmp_dict = {}
                                tmp_dict['mut iter'] = rnase_mut_iter_ii
                                tmp_dict['extinct_group'] = gro_ind
                                tmp_dict['orig RNase'] = rnase_orig
                                tmp_dict['mut RNase'] = rnase_mut
                                tmp_dict['RNase group'] = rnase_gro
                                # tmp_dict['total_non_comp_by_mut'] = iters_lst
                                if iter_ind not in analyze_extinct_group_all_iters_dict:
                                    analyze_extinct_group_all_iters_dict[iter_ind] = []
                                analyze_extinct_group_all_iters_dict[iter_ind].append(tmp_dict)

                        if len(mut_iters_arr_before_ext) == 1:
                            continue
                    cond = len(mut_iters_arr_before_ext) >1
                    ii = 0
                    while cond:
                        ii += 1
                        mut_iters_arr_before_ext[ii]
                        rnase_mut_iter_ii = mut_iters_arr_before_ext[ii]
                        ii_iter = rnase_mut_iter_ii
                        # start at the RNase mutation iteration, and at the extinction iteration
                        haps_in_gro = []
                        haps_in_ext = []
                        # haps_in_gro_ext = []
                        if gro_ind in list(reorder_mis_strict_groups_hap_all_iters_dict[ii_iter].keys()):
                            haps_in_gro = reorder_mis_strict_groups_hap_all_iters_dict[ii_iter][gro_ind]
                        haps_in_gro_ext1 = haps_in_gro
                        if gro_ind in list(reorder_mis_strict_group_extension_hap_all_iters_dict[ii_iter].keys()):
                            haps_in_ext = reorder_mis_strict_group_extension_hap_all_iters_dict[ii_iter][gro_ind]
                        haps_in_gro_ext2 = haps_in_ext
                        if len(haps_in_gro_ext1) > 0:
                            x = range(len(mut_rnases_dict[rnase_mut_iter_ii]))
                            for ii_mut_rnase in range(len(mut_rnases_dict[rnase_mut_iter_ii])):
                                # for cases of more than one RNase mutant in a certain iteration
                                rnase_orig = mut_rnases_dict[rnase_mut_iter_ii][ii_mut_rnase]['orig RNase']
                                rnase_mut = mut_rnases_dict[rnase_mut_iter_ii][ii_mut_rnase]['mutant RNase']
                                iters_lst = []
                                rnase_gro_lst = []
                                rnase_append_gro_lst = []

                                if rnase_orig not in reorder_mis_rnase_strict_groups_all_iters_dict[ii_iter] and \
                                        rnase_mut not in reorder_mis_rnase_strict_groups_all_iters_dict[ii_iter]:
                                    cond = False
                                    break

                                with open(output_dir + '/haplotypes_dict_' + str(ii_iter) + '.pkl', 'rb') as input:
                                    haplotypes_dict = pickle.load(input)
                                with open(output_dir + '/haplotypes_count_dict_' + str(ii_iter) + '.pkl', 'rb') as input:
                                    haplotypes_count_dict = pickle.load(input)
                                # check if male haplotypes in the group that is about to extinct pollinate the orig RNase
                                rnase_AA = AA_rnases_dict[rnase_orig]
                                is_orig_comp = np.zeros(len(haps_in_gro_ext1), dtype=bool)
                                for hap_ii, hap_ind in enumerate(haps_in_gro_ext1):
                                    haplotypes_dict[hap_ind]
                                    slfs_ind = haplotypes_dict[hap_ind].SLFsIndices
                                    slf_cond = False
                                    for slf_ii, slf_ind in enumerate(slfs_ind):
                                        slf_AA = AA_slfs_dict[slf_ind]
                                        E = np.sum(Model.eij[rnase_AA, slf_AA], axis=0)
                                        if E < interaction_thresh:
                                            # interacting SLF
                                            slf_cond = True
                                            is_orig_comp[hap_ii] = True
                                            break
                                # check if male haplotypes in the group that is about to extinct pollinate the mut RNase
                                rnase_AA = AA_rnases_dict[rnase_mut]
                                is_mut_comp = np.zeros(len(haps_in_gro_ext1), dtype=bool)
                                n_comp_mut = 0
                                n_non_comp_mut = 0
                                for hap_ii, hap_ind in enumerate(haps_in_gro_ext1):
                                    slfs_ind = haplotypes_dict[hap_ind].SLFsIndices
                                    slf_cond = False
                                    for slf_ii, slf_ind in enumerate(slfs_ind):
                                        slf_AA = AA_slfs_dict[slf_ind]
                                        E = np.sum(Model.eij[rnase_AA, slf_AA], axis=0)
                                        if E < interaction_thresh:
                                            # interacting SLF
                                            slf_cond = True
                                            is_mut_comp[hap_ii] = True
                                            n_comp_mut += haplotypes_count_dict[hap_ind]
                                            break
                                    if not slf_cond:
                                        n_non_comp_mut += haplotypes_count_dict[hap_ind]

                                if n_non_comp_mut/(n_comp_mut + n_non_comp_mut) > 0.75 and all(is_orig_comp):
                                    if rnase_orig in reorder_mis_rnase_strict_groups_all_iters_dict[ii_iter]:
                                        rnase_gro_lst.append(reorder_mis_rnase_strict_groups_all_iters_dict[ii_iter][rnase_orig])
                                    else:
                                        for app_gro, haps_lst in reorder_mis_strict_group_extension_hap_all_iters_dict[ii_iter].items():
                                            for hap in haps_lst:
                                                if rnase_orig == haplotypes_dict[hap].RNasesIndices[0]:
                                                    rnase_append_gro_lst.append(app_gro)
                                    if rnase_mut in reorder_mis_rnase_strict_groups_all_iters_dict[ii_iter]:
                                        rnase_gro_lst.append(reorder_mis_rnase_strict_groups_all_iters_dict[ii_iter][rnase_mut])
                                    else:
                                        for app_gro, haps_lst in reorder_mis_strict_group_extension_hap_all_iters_dict[ii_iter].items():
                                            for hap in haps_lst:
                                                if rnase_mut == haplotypes_dict[hap].RNasesIndices[0]:
                                                    rnase_append_gro_lst.append(app_gro)
                                    unique_gro, count = np.unique(rnase_gro_lst, return_counts=True)
                                    rnase_gro = unique_gro[np.argmax(count)]
                                    tmp_dict = {}
                                    tmp_dict['mut iter'] = rnase_mut_iter_ii
                                    tmp_dict['extinct_group'] = gro_ind
                                    tmp_dict['orig RNase'] = rnase_orig
                                    tmp_dict['mut RNase'] = rnase_mut
                                    tmp_dict['RNase group'] = rnase_gro
                                    if iter_ind not in analyze_extinct_group_all_iters_dict:
                                        analyze_extinct_group_all_iters_dict[iter_ind] = []
                                    analyze_extinct_group_all_iters_dict[iter_ind].append(tmp_dict)

                                cond = ii < len(mut_iters_arr_before_ext) - 1  # true if before going through all the neutral rnase mutants

                                if not cond:
                                    break
                        else:
                            cond = ii < len(
                                mut_iters_arr_before_ext) - 1  # true if before going through all the neutral rnase mutants

                            if not cond:
                                break
            with open(output_dir + '/analyze_extinct_group_all_iters_' + str(min_duration * 10) + '_dict.pkl', 'wb') as output:
                pickle.dump(analyze_extinct_group_all_iters_dict, output, pickle.HIGHEST_PROTOCOL)
            with open(output_dir + '/small_ext_groups_all_iters_dict_' + str(min_duration * 10) + '_dict.pkl', 'wb') as output:
                pickle.dump(small_ext_groups_all_iters_dict, output, pickle.HIGHEST_PROTOCOL)