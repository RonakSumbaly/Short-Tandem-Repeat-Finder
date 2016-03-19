import time
from helpers import *
from os.path import exists
from data_variables import *
from improved_helpers import *

logger.basicConfig(level=logger.INFO,format='> %(message)s')

def map_reads_snps(reads, hashed_reference_map, reference):
    """
    :return: map reads to get SNPs
    """
    if exists(SNP_file_path):
        logger.info("Loading SNPs from {}".format(SNP_file_path))
        SNPs = cPickle.load(open(SNP_file_path, "rb"))  # get existing SNPs if present
    else:
        logger.info("Finding SNPs in the donor")
        unmapped_reads = []
        start = time.clock()

        for count, read in enumerate(reads):
            if count % 1000 == 0 and count != 0:
                logger.info(
                    "Mapping reads to reference. Completed {} %".format(str(100 * count / float(len(reads)))[:5]))
                logger.info("Reads remaining : {0}".format(str(len(reads) - count)))

            read = read.strip().split(',')

            if len(read[0]) != 50 or len(read[1]) != 50:
                continue

            check_read_1 = map_read_to_reference(read[0], hashed_reference_map, reference)  # map read
            reverse_read_1 = False
            if not check_read_1:
                check_read_1 = map_read_to_reference(read[0][::-1], hashed_reference_map, reference) # map reverse read
                reverse_read_1 = True

            check_read_2 = map_read_to_reference(read[1], hashed_reference_map, reference)
            reverse_read_2 = False

            if not check_read_2:
                check_read_2 = map_read_to_reference(read[1][::-1], hashed_reference_map, reference)
                reverse_read_2 = True

            # only one of them is matched
            # if reversed then other one would not be reversed
            if (bool(check_read_1) == False and bool(check_read_2) == True) or \
                    (bool(check_read_1) == True and bool(check_read_2) == False):
                if bool(check_read_1):
                    if reverse_read_1:        recheck_read = read[1][::-1]
                    else:                     recheck_read = read[1]
                    check_position = check_read_1
                else:
                    if reverse_read_2:        recheck_read = read[0][::-1]
                    else:                     recheck_read = read[0]
                    check_position = check_read_2

                unmapped_reads.append([recheck_read[::-1], check_position - 200, check_position + 200])

        SNPs = get_snps(reference)
        logger.info("Dumping Unmapped Reads at {}".format(unmapped_read_file_path))
        cPickle.dump(unmapped_reads, open(unmapped_read_file_path, "wb"))
        logger.info("Dumping SNPs at {}".format(SNP_file_path))
        cPickle.dump(SNPs, open(SNP_file_path, "wb"))

    return SNPs


def map_reads_indels(reference):
    """
    :return: check unmapped reads for indels
    """
    if exists(INDEL_file_path):
        logger.info("Loading INDELs from {}".format(INDEL_file_path))
        indels = cPickle.load(open(INDEL_file_path, "rb"))
        return indels

    if exists(unmapped_read_file_path):
        logger.info("Loading Unmapped Reads from {}".format(unmapped_read_file_path))
        unmapped_reads = cPickle.load(open(unmapped_read_file_path, "rb"))
        logger.info("Checking {} unmapped reads for insertions / deletions".format(len(unmapped_reads)))
        start = time.clock()

        for count, unmapped in enumerate(unmapped_reads):
            if count % 1000 == 0 and count != 0:
                logger.info("Mapping reads to reference. Completed {} %".format(
                    str(100 * count / float(len(unmapped_reads)))[:5]))
                logger.info("Reads remaining : {0}".format(str(len(unmapped_reads) - count)))

            check_for_indels(unmapped[0], reference, unmapped[1], unmapped[2])

        ins = get_insertions()
        DEL = get_deletions()

        logger.info("Dumping INDELs at {}".format(INDEL_file_path))
        cPickle.dump((ins, DEL), open(INDEL_file_path, "wb"))

        return ins, DEL
    else:
        logger.info("Run STR Code first to get unmapped reads")
        exit()


if __name__ == "__main__":
    start_time = time.clock()

    logger.info("Reading Reference File : {}".format(reference_file_name[dataset_choice]))
    reference = read_reference_file(reference_file_path)  # read the reference file
    logger.info("Length of Reference Sequence : {}".format(len(reference)))

    logger.info("Reading Reads File : {}".format(reads_file_name[dataset_choice]))
    reads = read_reads_file(reads_file_path) # read the read file
    logger.info("Number of Paired-End Reads : {}".format(len(reads)))

    hashed_reference_map = created_hashed_map(reference, hashed_file_path)  # create hash map of reference
    SNPs = map_reads_snps(reads, hashed_reference_map, reference)  # get SNPs by mapping reads
    SNPs = process_snps(SNPs)  # process the SNPs

    INDELs = map_reads_indels(reference)  # get INDELs by checking unmapped reads
    processed_INDELs = (process_indels(INDELs[0]), process_indels(INDELs[1]))  # process the INDELs

    donor,reference_location = create_donor_sequence(reference, SNPs, processed_INDELs)  # recreate donor sequence

    STRs = get_tandem_repeats(donor)  # get STRs from donor genome
    STRs = find_reference_position(STRs, reference_location)  # get corresponding reference location

    logger.info("Total Number of SNPs : {}".format(len(SNPs)))
    logger.info("Total Number of INS : {}".format(len(processed_INDELs[1])))
    logger.info("Total Number of DEL : {}".format(len(processed_INDELs[0])))
    logger.info("Total Number of STRs : {}".format(len(STRs)))

    # combine outputs into dictionary
    output = defaultdict()
    output['SNP'] = SNPs
    output['STR'] = STRs
    output['DEL'] = processed_INDELs[0]
    output['INS'] = processed_INDELs[1]

    output_to_file(output, "improved_" + file_name[dataset_choice])
    logger.info("Time to execute {} secs".format(time.clock() - start_time))
    logger.info("Process completed")
