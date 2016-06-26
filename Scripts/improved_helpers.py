import os
import cPickle
import textwrap
import itertools
from helpers import *
from smith_waterman_align import *
from collections import defaultdict
import numpy as np

logger.basicConfig(level=logger.INFO,format='> %(message)s')

READ_LENGTH = 50
SNP_THRESHOLD = 5
KEY_LENGTH = 5
SNPS = []
INSERTIONS = []
DELETIONS = []
VARIATIONS = defaultdict()


def created_hashed_map(reference, hash_file_path):
    """
    :param reference: reference genome
    :param hash_file_path: location of hashtable
    :return: hashed reference genome
    """
    if os.path.exists(hash_file_path):
        logger.info("Loading hashed file from location {}".format(hash_file_path))
        hashed_reference = cPickle.load(open(hash_file_path, 'rb'))  # load from hard-disk
    else:
        logger.info("Creating hashed file from reference")
        hash_size = READ_LENGTH / KEY_LENGTH
        hashed_reference = defaultdict(list)

        for i in range(len(reference) - hash_size):
            if i % 100000 == 0:
                logger.info("Creating Hashed Reference Map. Completed {} %".format(
                    str(100 * i / float(len(reference) - hash_size))[:5]))
            hashed_reference[reference[i: i + hash_size]].append(i)

        logger.info("Dumping hashed file at {}".format(hash_file_path))
        cPickle.dump(hashed_reference, open(hash_file_path, "wb"))  # dump on hard-disk

    return hashed_reference


def match_read_with_reference(reference, read, start_ref_pos, end_ref_pos):
    """
    :return: get all variations between read and part of reference
    """
    small_variations = []
    check_reference = reference[start_ref_pos:end_ref_pos]
    for i in range(len(check_reference)):
        if check_reference[i] != read[i]:
            small_variations.append([start_ref_pos + i, read[i]])

    return small_variations


def map_read_to_reference(read, hashed_reference_map, reference):
    """
    :return: start location of read in the reference else false if cannot find location
    """
    read_hash = textwrap.wrap(read, len(read) / KEY_LENGTH)
    for i, part_read in enumerate(read_hash):
        if hashed_reference_map.has_key(part_read):
            for positions in hashed_reference_map[part_read]:
                start_ref_pos = positions - (i * (len(read) / KEY_LENGTH))
                end_ref_pos = min(len(reference), start_ref_pos + len(read))
                score = match_read_with_reference(reference, read, start_ref_pos, end_ref_pos)
                if len(score) <= KEY_LENGTH - 1:
                    add_variations(score)
                    return start_ref_pos

    return False


def add_variations(scores):
    """
    :return: add all variations at specific position
    """
    for score in scores:
        if VARIATIONS.has_key(score[0]):
            temp = VARIATIONS[score[0]]
            temp.append(score[1])
            VARIATIONS[score[0]] = temp
        else:
            VARIATIONS[score[0]] = [score[1]]


def get_snps(reference):
    """
    :return: all SNPS that have coverage greater than threshold
    """
    logger.info("Finding SNPs using variations")
    for key in VARIATIONS:
        snp = max(VARIATIONS[key], key=VARIATIONS[key].count)
        if VARIATIONS[key].count(snp) >= SNP_THRESHOLD:
            SNPS.append([reference[key], snp, key])
    return SNPS


def get_tandem_repeats(donor):
    """
    :return: STRs in the donor sequence
    """
    logger.info("Getting Short Tandem Repeats")
    tandem_repeats = get_reference_tandem_repeats(''.join(donor))  # get short tandem repeats
    STRs = preprocess_tandems(tandem_repeats)  # pre-process tandem repeats

    return STRs


def check_for_indels(recheck_read, reference, start_ref, end_ref):
    """
    :return: insertions and deletions between read and reference
    """
    align_seq_1, align_seq_2 = waterman_algorithm(reference[start_ref:end_ref], recheck_read)

    if bool(align_seq_1) and bool(align_seq_2):
        if "-" in align_seq_1 and '-' not in align_seq_2:
            check_for_insertions(align_seq_1, align_seq_2, start_ref, reference[start_ref:end_ref])
        elif '-' in align_seq_2 and '-' not in align_seq_1:
            check_for_deletions(align_seq_1, align_seq_2, start_ref, reference[start_ref:end_ref])


def create_donor_sequence(reference, SNPs, INDELs):
    """
    :return: donor sequence and corresponding location in reference
    """
    logger.info("Reassembling Donor Sequence")
    donor = list(reference)
    reference_location = np.arange(1, len(reference))

    # change SNPs in the donor sequence
    logger.info("Modifying SNPs")
    for s in SNPs:
        donor[int(s[2])] = s[1]

    # insertions
    for count, ins in enumerate(INDELs[1]):
        reference_location[int(ins[1]):] = np.add(reference_location[int(ins[1]):],
                                                  (len(reference_location) - int(ins[1])) * [len(ins[0])])
        donor.insert(int(ins[1]), ins[0])

    # deletions
    # reducing accuracy :'( by 5 points
    for dels in INDELs[0]:
        delete_location = reference_location[dels[1]] - 1
        reference_location[dels[1]:] = np.subtract(reference_location[dels[1]:],
                                                   (len(reference_location) - dels[1]) * [len(dels[0])])
        del donor[delete_location:delete_location + len(dels[0])]

    return donor, reference_location


def process_snps(snps):
    """
    :return: processed SNPs
    """
    processed_snps = []
    count = 0
    for snp in snps:
        if count == 0:
            processed_snps.append(snp)
            count += 1
        else:
            if abs(int(processed_snps[count - 1][2]) - int(snp[2])) <= 10:
                continue
            else:
                processed_snps.append(snp)
                count += 1

    return processed_snps


def process_indels(indel):
    """
    :return: processed indels
    """
    processed_indels = []
    count = 0
    for ins in indel:
        if count == 0:
            processed_indels.append(ins)
            count += 1
        else:
            if processed_indels[count - 1][0] == ins[0] or abs(int(processed_indels[count - 1][1]) - int(ins[1])) <= 10:
                continue
            else:
                processed_indels.append(ins)
                count += 1

    return processed_indels


def find_reference_position(STRs, reference_location):
    """
    :return: real STR position in reference
    """
    processed_STRs = []
    reference_location = reference_location.tolist()
    reference_location_set = set(reference_location)
    for STR in STRs:
        if STR[1] in reference_location_set:
            new_location = reference_location.index(STR[1])
        # messed up the index
        elif STR[1] + 1 in reference_location_set:
            new_location = reference_location.index(STR[1] + 1)
        # messed up the index
        elif STR[1] - 1 in reference_location_set:
            new_location = reference_location.index(STR[1] - 1)
        # just omit the STRs - better than finding their approximate location
        else:
            continue

        processed_STRs.append([STR[0],new_location])

    return processed_STRs


def add_insertions(start_pos, insert_string):
    INSERTIONS.append((insert_string, start_pos))


def add_deletions(start_pos, insert_string):
    DELETIONS.append((insert_string, start_pos))


def get_insertions():
    return list(k for k, _ in itertools.groupby(sorted(INSERTIONS, key=itemgetter(1))))


def get_deletions():
    return list(k for k, _ in itertools.groupby(sorted(DELETIONS, key=itemgetter(1))))
