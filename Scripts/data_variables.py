from os.path import join

# file paths and names
dataset_choice = 1  # [0 - Practice, 1 - Undergraduate, 2 - Graduate]
base_folder = "../../Data/Data - Homework 2"  # location of the data files
file_name = ['practice_E_1_chr_1.txt', 'hw2undergrad_E_2_chr_1.txt', 'hw2grad_M_1_chr_1.txt']  # base file names

reference_file_name = ['Practice/ref_practice_E_1_chr_1.txt', 'Undergraduate/ref_hw2undergrad_E_2_chr_1.txt', 'Graduate/ref_hw2grad_M_1_chr_1.txt']  # reference file
reads_file_name = ['Practice/reads_practice_E_1_chr_1.txt', 'Undergraduate/reads_hw2undergrad_E_2_chr_1.txt', 'Graduate/reads_hw2grad_M_1_chr_1.txt']  # reads file
reference_file_path = join(base_folder, reference_file_name[dataset_choice])
reads_file_path = join(base_folder, reads_file_name[dataset_choice])

hashed_file_path = join(base_folder, "Hashed/hashed_" + file_name[dataset_choice])  # hashed map location
SNP_file_path = join(base_folder, "SNPs/snps_" + file_name[dataset_choice])  # snps location
INDEL_file_path = join(base_folder, "INDELs/indels_" + file_name[dataset_choice])  # indel location
unmapped_read_file_path = join(base_folder, "Unmapped/unmapped_" + file_name[dataset_choice])  # unmapped reads location
