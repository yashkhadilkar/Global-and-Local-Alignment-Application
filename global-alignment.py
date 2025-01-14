import numpy as np

'''
The compute_alignment_matrix function computes the global alignment matrix for two nucleotide sequences 
using dynamic programming.

Input:
seq1: First nucleotide sequence
seq2: Second nucleotide sequence
match_score: Score for a match
mismatch_score: Score for a mismatch
gap_penalty: Penalty for introducing a gap

Output:
alignment_matrix: The matrix containing optimal alignment scores
'''
def compute_alignment_matrix(seq1, seq2, match_score, mismatch_score, gap_penalty):
    # Create the alignment matrix with dimensions based on input sequences
    alignment_matrix = np.zeros((len(seq2), len(seq1)), dtype=int)

    # Initialize the first row and column with gap penalties
    for i in range(len(seq2)):
        alignment_matrix[i][0] = i * gap_penalty
    for j in range(len(seq1)):
        alignment_matrix[0][j] = j * gap_penalty

    # Fill the matrix using dynamic programming
    for i in range(1, len(seq2)):
        for j in range(1, len(seq1)):
            if seq1[j] == seq2[i]:
                score_diagonal = alignment_matrix[i-1][j-1] + match_score
            else:
                score_diagonal = alignment_matrix[i-1][j-1] + mismatch_score

            score_up = alignment_matrix[i-1][j] + gap_penalty
            score_left = alignment_matrix[i][j-1] + gap_penalty

            # Assign the maximum score to the current cell
            alignment_matrix[i][j] = max(score_diagonal, score_up, score_left)

    return alignment_matrix

'''
The trace_alignment function traces back through the alignment matrix to find the optimal alignment
and build the aligned sequences.

Input:
alignment_matrix: The computed alignment matrix
seq1: First nucleotide sequence
seq2: Second nucleotide sequence

Output:
aligned_seq1: The aligned version of sequence 1
aligned_seq2: The aligned version of sequence 2
alignment_path: The sequence of scores from the optimal alignment path
'''
def trace_alignment(alignment_matrix, seq1, seq2):
    aligned_seq1 = []
    aligned_seq2 = []
    alignment_path = []

    i = len(seq2) - 1
    j = len(seq1) - 1

    # Trace the optimal alignment path starting from the bottom-right corner
    while i > 0 or j > 0:
        current_score = alignment_matrix[i][j]
        alignment_path.append(current_score)

        if i > 0 and j > 0 and (seq1[j] == seq2[i] or alignment_matrix[i][j] == alignment_matrix[i-1][j-1]):
            aligned_seq1.append(seq1[j])
            aligned_seq2.append(seq2[i])
            i -= 1
            j -= 1
        elif j > 0 and alignment_matrix[i][j] == alignment_matrix[i][j-1] + gap_penalty:
            aligned_seq1.append(seq1[j])
            aligned_seq2.append('-')
            j -= 1
        else:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[i])
            i -= 1

    # Reverse the aligned sequences and path
    aligned_seq1.reverse()
    aligned_seq2.reverse()
    alignment_path.reverse()

    return aligned_seq1, aligned_seq2, alignment_path

'''
The display_alignment function prints the alignment matrix, optimal path, and aligned sequences.

Input:
alignment_matrix: The computed alignment matrix
alignment_path: The optimal path through the alignment matrix
seq1: First nucleotide sequence
seq2: Second nucleotide sequence
aligned_seq1: The aligned version of sequence 1
aligned_seq2: The aligned version of sequence 2
'''
def display_alignment(alignment_matrix, alignment_path, seq1, seq2, aligned_seq1, aligned_seq2):
    print('Alignment Matrix:')
    print_matrix(alignment_matrix, seq1, seq2)
    print()

    print('Optimal Alignment Path:')
    print(alignment_path)
    print()

    print('Aligned Sequences:')
    print(' '.join(aligned_seq1))
    
    alignment_indicators = ''
    for a, b in zip(aligned_seq1, aligned_seq2):
        if a == b:
            alignment_indicators += '| '
        elif a == '-' or b == '-':
            alignment_indicators += '  '
        else:
            alignment_indicators += '. '
    print(alignment_indicators)
    
    print(' '.join(aligned_seq2))
    print()

def print_matrix(alignment_matrix, seq1, seq2):
    print('  |', ' | '.join(seq1), '|')
    print('--' + '---' * len(seq1) + '+')
    for i in range(len(seq2)):
        row = [str(alignment_matrix[i][j]).rjust(2) for j in range(len(seq1))]
        print(f"{seq2[i]} | {' | '.join(row)} |")
        print('--' + '---' * len(seq1) + '+')

def format_sequence(seq):
    formatted_seq = ['-'] + list(seq)
    return formatted_seq


