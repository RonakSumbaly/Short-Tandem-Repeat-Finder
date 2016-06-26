import re
import logging as logger
from itertools import groupby
from itertools import product
from operator import itemgetter
from collections import defaultdict

logger.basicConfig(level=logger.INFO,format='> %(message)s')

STR_THRESHOLD = 5  # threshold limit of STR length

def read_reads_file(reads_file_loc):
    """
    :param reads_file_loc: location of the reads file
    :return: list of all the paired-ends reads
    """
    input_reads = []
    omit_first_line = True

    read_file = open(reads_file_loc, 'r')
    while True:
        lines = read_file.readlines(100000)
        if not lines:
            break
        if omit_first_line:  # omit first line of file
            omit_first_line = False
            del lines[0]

        input_reads.extend(lines)

    return input_reads


def read_reference_file(reference_file_name):
    """
    :param reference_file_name: location of the reference file
    :return: reference genome
    """
    reference_seq = ''
    omit_first_line = True

    with open(reference_file_name, 'r') as f:
        for line in f.readlines():
            if omit_first_line:  # omit first line of file
                omit_first_line = False
                continue
            reference_seq += line.strip()

    return reference_seq


def output_to_file(output, file_name):
    """
    :param output: dictionary comprises of output variations
    :param file_name: output file name
    :return: creates output file with all variants
    """
    fopen = open(file_name, "wb")

    fopen.write(">" + file_name.split('.')[0] + "\n")
    for key in output:
        fopen.write(">{}".format(key) + "\n")
        if key == "SNP":
            for i in output[key]:
                fopen.write(i[0] + "," + i[1] + "," + str(i[2]) + "\n")
        elif key == "STR":
            for i in output[key]:
                fopen.write(i[0] + "," + str(i[1]) + "\n")
        elif key == "INS" or key == "DEL":
            for i in output[key]:
                fopen.write(i[0] + "," + str(i[1]) + "\n")

    tails = ('>' + x for x in ('CNV', 'ALU', 'INV'))
    fopen.write('\n'.join(tails))


def tandem_combos(repeat):
    """
    :param repeat: length of repeats
    :return: combination of all ACGT repeats of specified length
    """
    word_list = []
    for p in product('ACGT', repeat=repeat):
        word_list.append("".join(p))
    return word_list


def rotate(str, n):
    return str[n:] + str[:n]  # rotate the string


def preprocess_tandems(tandem_repeats):
    """
    :param tandem_repeats: STRs
    :return: processed STRs - remove ones which are close by and similar
    """
    logger.info("Preprocessing Tandem Repeats - Removing outliers")

    # expand short tandem repeats
    expand_tandem_repeats = []
    for STR in tandem_repeats:
        for lst in tandem_repeats[STR]:
            expand_tandem_repeats.append([STR * lst[1], lst[0]])

    # get indices of all short tandem repeats
    clean = list(k for k, _ in groupby(expand_tandem_repeats))
    indices = [expand_tandem_repeats[1] for expand_tandem_repeats in clean]
    unique_indices = []

    # combine indices that lie close to each other
    for expand_tandem_repeats, g in groupby(enumerate(sorted(indices)), lambda (i, x): i - x):
        unique_indices.append(map(itemgetter(1), g)[0])

    unique_indices = list(set(unique_indices))

    # only consider indices that occur first
    processed_tandem_repeats = []
    for i in clean:
        if i[1] in unique_indices:
            processed_tandem_repeats.append(i)
            del unique_indices[unique_indices.index(i[1])]

    return processed_tandem_repeats


def get_reference_tandem_repeats(genome):
    """
    :param genome: genome sequence with STRs
    :return: list of all STRs and their position
    """
    short_tandem_repeats = defaultdict(list)

    # loop through each length of STR
    for STR_LENGTH in xrange(2, 6):
        logger.info("Checking STR of Length : {}".format(STR_LENGTH))

        word_list = tandem_combos(STR_LENGTH)  # get all combos of STRs for specify length
        temp_tandem_repeats = defaultdict(list)

        for STR in word_list:  # consider each STR separately
            indices = []

            # append all occurrences of the STR in reference
            for m in re.finditer(STR * STR_THRESHOLD, genome):
                indices.append(m.start())

            threshold = STR_LENGTH * STR_THRESHOLD
            clean_indices = []

            for index in range(len(indices)):
                if index == 0:
                    clean_indices.append(indices[index])
                else:
                    if indices[index] - indices[index - 1] > threshold:
                        clean_indices.append(indices[index])

            # loop through each start-index of STR
            for index in clean_indices:
                cur_pos = index
                repeats = 0

                # calculate actual number of repeats in the reference
                while True:
                    check = genome[cur_pos:cur_pos + STR_LENGTH]
                    if check == STR:
                        repeats += 1
                    else:
                        break
                    cur_pos += STR_LENGTH

                # append temporary found tandem repeats
                temp_tandem_repeats[STR].append((index, repeats))

        # loop through each temporary found tandem repeat
        # delete complementary repeats whose index is greater
        for STR in temp_tandem_repeats.keys():
            current = temp_tandem_repeats[STR]
            for each_index in current:
                flag = False
                for com in range(1, len(STR)):
                    complement = temp_tandem_repeats[rotate(STR, com)]
                    other_index = 0
                    while other_index < len(complement):
                        if abs(each_index[0] - complement[other_index][0]) < 10:
                            if each_index[0] > complement[other_index][0]:
                                flag = False
                                short_tandem_repeats[rotate(STR, com)].append(complement[other_index])
                                break
                            else:
                                temp_tandem_repeats[rotate(STR, com)].remove(complement[other_index])
                                continue
                        other_index += 1
                if not flag:
                    short_tandem_repeats[STR].append(each_index)

    return short_tandem_repeats
