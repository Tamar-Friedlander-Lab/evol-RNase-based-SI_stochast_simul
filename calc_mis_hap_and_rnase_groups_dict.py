import pickle
import os

d_iters = 10
os.chdir("..")
os.chdir("FlowersEvolutionsSimulation_ver13_data")
directory = "output_data_10haplotypes_k5_18AAs_E_thresh-6_p10-4"
subdirectories = next(os.walk(directory))[1]
num_files = len(subdirectories)
for subdirectory in subdirectories:
    output_dir = directory + "/" + subdirectory
    if os.path.isfile(output_dir + '/mis_groups_hap_all_iters_dict.pkl'):
        print(os.path.join(subdirectory))
        
        with open(output_dir + '/mis_groups_hap_all_iters_dict.pkl', 'rb') as input:
            mis_groups_hap_all_iters_dict= pickle.load(input)
    
        with open(output_dir + '/mis_groups_rnase_all_iters_dict.pkl', 'rb') as input:
            mis_groups_rnase_all_iters_dict= pickle.load(input)
    
        with open(output_dir + '/RNasesInHaplotypes_all_iters_Dict.pkl', 'rb') as input:
            RNasesInHaplotypes_all_iters_Dict= pickle.load(input)

        mis_rnase_groups_all_iters_dict = {}
        mis_hap_groups_all_iters_dict = {}
        for iter_ind in mis_groups_hap_all_iters_dict.keys(): 
    
            mis_rnase_groups_dict = {}
            mis_hap_groups_dict = {}
            mis_groups_hap_dict = mis_groups_hap_all_iters_dict[iter_ind]
            mis_groups_rnase_dict = mis_groups_rnase_all_iters_dict[iter_ind]
            for group_ind, rnases_lst in mis_groups_rnase_dict.items():
                for rnase in rnases_lst:
                    mis_rnase_groups_dict[rnase] = group_ind
    
            for group_ind, haps_lst in mis_groups_hap_dict.items():
                for hap in haps_lst:
                    mis_hap_groups_dict[hap] = group_ind
    
            mis_hap_groups_all_iters_dict[iter_ind] = mis_hap_groups_dict
            mis_rnase_groups_all_iters_dict[iter_ind] = mis_rnase_groups_dict

        with open(output_dir + '/mis_hap_groups_all_iters_dict.pkl', 'wb') as output:
            pickle.dump(mis_hap_groups_all_iters_dict, output, pickle.HIGHEST_PROTOCOL)

        with open(output_dir + '/mis_rnase_groups_all_iters_dict.pkl', 'wb') as output:
            pickle.dump(mis_rnase_groups_all_iters_dict, output, pickle.HIGHEST_PROTOCOL)