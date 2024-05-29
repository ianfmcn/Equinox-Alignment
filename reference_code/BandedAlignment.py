from Bio import SeqIO
import argparse

def BandedAlignment(ref, read, match_reward, mismatch_penalty, indel_penalty, band_parameter):
    ref_len = len(ref) + 1
    read_len = len(read) + 1
    current_row = [0] * (read_len)
    previous_row = [0] * (read_len)

    max_score = 0
    max_i = 0
    max_j = 0
    

    for i in range(1, ref_len):
        for j in range(1, read_len):
            diff_down = abs(i - 1 - j)
            if diff_down < band_parameter:
                down = previous_row[j] + indel_penalty
            else:
                down = -10000

            diff_side = abs(i - (j - 1))
            if diff_side < band_parameter:
                side = current_row[j - 1] + indel_penalty
            else: 
                side = -10000
            if ref[i - 1] == read[j - 1]:
                diag = previous_row[j - 1] + match_reward
            else:
                diag = previous_row[j - 1] + mismatch_penalty

            current_row[j] = max(down, side, diag, 0)
            if current_row[j] > max_score:
                max_score = current_row[j]
                max_i = i
                max_j = j
        
        previous_row, current_row = current_row, [0] * (read_len + 1) 
    return max_score, max_i, max_j

ignore = """
def find_start_indices(ref, read, match_reward, mismatch_penalty, indel_penalty, band_parameter, end_i, end_j):
    ref_reversed = ref[:end_i][::-1]
    read_reversed = read[:end_j][::-1]

    _, start_i_reversed, start_j_reversed = BandedAlignment(ref_reversed, read_reversed, match_reward, mismatch_penalty, indel_penalty, band_parameter)

    start_i = end_i - start_i_reversed + 1
    start_j = end_j - start_j_reversed + 1

    return start_i, start_j
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('reference', type=str)
    parser.add_argument('reads', type=str)
    parser.add_argument('bandwidth', type=int)
    parser.add_argument('-m', '--match', type=float, default=1)
    parser.add_argument('-s', '--mismatch', type=float, default=-1)
    parser.add_argument('-d', '--indel', type=float, default=-1)

    args = parser.parse_args()

    # Read and parse input files
    reference_sequences = {}
    for record in SeqIO.parse(args.reference, "fasta"):
        reference_sequences[record.id] = str(record.seq)

    reads = []
    for record in SeqIO.parse(args.reads, "fastq"):
        reads.append(str(record.seq))

    # Perform alignment for each read
    for read in reads:
        best_score = float('-inf')
        best_i = ""
        best_j = ""
        best_ref = ""
        start_i = 0
        start_j = 0

        for ref_id, ref_seq in reference_sequences.items():
            max_score, max_i, max_j = BandedAlignment(ref_seq, read, args.match, args.mismatch, args.indel, args.bandwidth)
            if max_score > best_score:
                best_score = max_score
                best_i = max_i
                best_j = max_j
                best_ref = ref_id

        if best_score > float('-inf'):
            start_i, start_j = find_start_indices(reference_sequences[best_ref], read, args.match, args.mismatch, args.indel, args.bandwidth, best_i, best_j)

        print(f"Read: {read}")
        print(f"Reference sequence: {best_ref}")
        print(f"Score: {best_score}")
        print(f"Alignment starts at ref index: {start_i} and read index: {start_j}")
        print(f"Alignment ends at ref index: {best_i} and read index: {best_j}")
        print("\n")

if __name__ == "__main__":
    main() """