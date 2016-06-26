import time
from helpers import *
import logging as logger
from data_variables import *
from collections import defaultdict

logger.basicConfig(level=logger.INFO,format='> %(message)s')

if __name__ == "__main__":
    start_time = time.clock()

    logger.info("Reading Reference File : {}".format(reference_file_name[dataset_choice]))
    reference = read_reference_file(reference_file_path)  # read the reference file
    logger.info("Length of Reference Sequence : {}".format(len(reference)))

    tandem_repeats = get_reference_tandem_repeats(reference)  # get short tandem repeats
    STRs = preprocess_tandems(tandem_repeats)  # pre-process tandem repeats

    logger.info("Total Number of STRs found : {}".format(len(STRs)))

    # combine outputs into dictionary
    output = defaultdict()
    output['STR'] = STRs
    output['SNP'] = []
    output['INS'] = []
    output['DEL'] = []

    output_to_file(output, "baseline_" + file_name[dataset_choice])  # print output to file

    logger.info("Time to execute {} secs".format(time.clock() - start_time))
    logger.info("Process completed")



