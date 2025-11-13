import pandas as pd
import numpy as np
import os
import time
import glob
from covid_fasta_cleaner import CovidFastaProcessor, Sequence
from primer_loader import PrimerOperation

# ----------------------------
# Classes
# ----------------------------
class GenomeData:
    def __init__(self):
        self.country = ""
        self.date = ""
        self.sequence = ""
        self.length = 0
        self.id = ""
        self.virus_type = ""
        self.lineage = ""

# ----------------------------
# Global variables
# ----------------------------
primer_op = PrimerOperation()
ROW_COUNT = 1000

primer_groups, primers_id = primer_op.load_primers_from_file("data/raw/primers_info.xlsx")
primers_id = np.atleast_1d(primers_id)

stat_country = np.empty((ROW_COUNT, len(primer_groups) + 1, 2), dtype=object)
stat_date = np.empty((ROW_COUNT, len(primer_groups) + 1, 2), dtype=object)
stat_virus_type = np.empty((ROW_COUNT, len(primer_groups) + 1, 2), dtype=object)
stat_lineage = np.empty((ROW_COUNT, len(primer_groups) + 1, 2), dtype=object)
match_result = np.empty((ROW_COUNT, len(primer_groups) + 1, 2), dtype=object)
reasons_result = np.empty((22000, len(primers_id) + 1), dtype=object)

primer_group_names = ["Country"]  # will append groups later

# ----------------------------
# Initialization
# ----------------------------
def init_primer_groups():
    global primer_group_names
    primer_group_names += list(primer_groups.keys())

# ----------------------------
# Utility functions
# ----------------------------
def get_temperature(sequence):
    try:
        from Bio.SeqUtils import MeltingTemp
        return MeltingTemp.Tm_NN(sequence)
    except:
        return 100

def match_prob(seqA, seqB, max_mismatch=5):
    mismatches = sum(1 for a, b in zip(seqA, seqB) if a != b)
    return (mismatches <= max_mismatch, mismatches if mismatches <= max_mismatch else 500)

def match_non_prob(seqA, seqB):
    counter = 0
    mismatch = 0
    for a, b in zip(seqA, seqB):
        if counter < 5:
            if a == b:
                counter += 1
                continue
            else:
                return False, 1000
        else:
            if a != b:
                mismatch += 1
                if mismatch > 3:
                    return False, 600
    return True, mismatch

def get_align_score(seq1, seq2, primer_temp, primer_candidate):
    match_func = match_prob if "p" in primer_candidate.subgroup else match_non_prob
    match, score = match_func(seq1, seq2)

    temp1 = get_temperature(seq1)
    temp2 = get_temperature(seq2)
    diff_temp = abs(temp1 - temp2)

    if match and diff_temp <= 8:
        return True, "", score
    reason = f"Score {score}\n{seq1}\n{seq2}" if not match else f"Temp: {diff_temp}\n{seq1}\n{seq2}"
    return False, reason, score

def add_and_get_index(item, array):
    for i in range(array.shape[0]):
        if array[i, 0, 0] == item:
            return i
        if array[i, 0, 0] is None:
            array[i, 0, 0] = array[i, 0, 1] = item
            return i

def add_and_get_index_reason(item, array):
    for i in range(array.shape[0]):
        if array[i, 0] == item:
            return i
        if array[i, 0] is None:
            array[i, 0] = item
            return i

def update_cell(cell, value):
    return (0 if cell is None else cell) + value

def update_statistic(genome, col_idx, value):
    col_idx += 1
    country_idx = add_and_get_index(genome.country, stat_country)
    date_idx = add_and_get_index(genome.date, stat_date)
    virus_idx = add_and_get_index(genome.virus_type, stat_virus_type)
    lineage_idx = add_and_get_index(genome.lineage, stat_lineage)

    for array, idx in [(stat_country, country_idx), (stat_date, date_idx),
                       (stat_virus_type, virus_idx), (stat_lineage, lineage_idx)]:
        array[idx, col_idx, 0] = update_cell(array[idx, col_idx, 0], 1)
        array[idx, col_idx, 1] = update_cell(array[idx, col_idx, 1], value)

def update_match(genome, col_idx, value):
    col_idx += 1
    match_idx = add_and_get_index(genome.id, match_result)
    match_result[match_idx, col_idx, 0] = value
    match_result[match_idx, col_idx, 1] = value

def get_quarter(date_str):
    date = pd.to_datetime(date_str)
    return f"{date.year}.Q{date.quarter}"

def update_reject_list(genome):
    with open("reject_list.txt", "a") as f:
        f.write(f"{genome.id},{genome.date},{genome.virus_type},{genome.sequence}\n")

def update_reason(genome_id, primer_id, reason):
    col_idx_array = np.where(primers_id == primer_id)[0]
    if len(col_idx_array) == 0:
        return  # skip if primer not found
    col_idx = col_idx_array[0]
    row_idx = add_and_get_index_reason(genome_id, reasons_result)
    reasons_result[row_idx, col_idx + 1] = reason

# ----------------------------
# Alignment functions
# ----------------------------
def align_sequence(genome, primer_candidate):
    primer_match = False
    for seq, length, pid in zip(primer_candidate.primers, primer_candidate.lengths, primer_candidate.ids):
        start = primer_candidate.start_position - Sequence.VALID_RANGE
        end = primer_candidate.start_position + Sequence.VALID_RANGE
        best_score = 1000
        reason = ""

        for pos in range(start, end):
            if pos < 0:
                continue
            genome_seq = genome.sequence[pos: pos + length]
            match, reason_tmp, score = get_align_score(genome_seq, seq, primer_candidate.temperature, primer_candidate)
            if score < best_score:
                best_score = score
                reason = reason_tmp
                if match:
                    primer_match = True
                    break

        if not primer_match:
            update_reason(genome.id, pid, reason)
    return primer_match

def do_align_reason_report(genome):
    group_match = []
    for i, key in enumerate(primer_groups.keys()):
        primers = primer_groups[key]
        all_match = [align_sequence(genome, pc) for pc in primers]
        match = any(all_match) if len(primers) == 2 else all(all_match)

        update_statistic(genome, i, int(match))
        update_match(genome, i, int(match))

        group_match.append(match)

    if not any(group_match):
        update_reject_list(genome)
    return any(group_match)

# ----------------------------
# Output functions
# ----------------------------
def save_results():
    os.makedirs("results", exist_ok=True)
    arrays = [stat_country, stat_date, stat_virus_type, stat_lineage, match_result]
    sheet_names = ["Country_Count", "Date_Count", "VirusType", "Lineage", "Match"]
    with pd.ExcelWriter("results/result.xlsx", engine="xlsxwriter") as writer:
        for arr, name in zip(arrays, sheet_names):
            df_count = pd.DataFrame(arr[:, :, 0], columns=primer_group_names)
            df_match = pd.DataFrame(arr[:, :, 1], columns=primer_group_names)
            df_count.to_excel(writer, sheet_name=name, index=False)
            df_match.to_excel(writer, sheet_name=name + "_Match", index=False)

def save_reasons():
    os.makedirs("results", exist_ok=True)
    columns = ["GenomeID"] + list(primers_id)
    pd.DataFrame(reasons_result, columns=columns).to_excel("results/report_reason.xlsx", index=False)

# ----------------------------
# Main processing
# ----------------------------
def process_file(filename):
    with open(filename, "r") as f:
        f.readline()  # skip header
        for line in f:
            data = line.strip().split(",")
            genome = GenomeData()
            genome.country = data[0]
            genome.sequence = data[1]
            genome.id = data[4]
            genome.date = get_quarter(data[5])
            genome.lineage = data[12]
            genome.virus_type = data[11]

            try:
                do_align_reason_report(genome)
                print(f"{genome.id},{data[5]}")
            except Exception as e:
                print(e)

# ----------------------------
# Entry point
# ----------------------------
if __name__ == "__main__":
    init_primer_groups()
    input_files = glob.glob("data/clean/*.csv")

    start_time = time.time()
    for file in input_files:
        process_file(file)

    save_results()
    save_reasons()
    print("Processing complete")
    print("Time elapsed:", time.time() - start_time)
