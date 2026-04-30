import sys

# Scoring Scheme
match = 1
mismatch = -1
gap = -2


#FASTA READER
def read_fasta(file):
    seq = ""
    with open(file) as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip()
    return seq
    
# ---- READ FILES INSTEAD OF INPUT ----
if len(sys.argv)!= 3:
    print("Usage: python script.py seq1.txt seq2.txt")
    sys.exit()

#Read sequence
seq1 = read_fasta(sys.argv[1])
seq2 = read_fasta(sys.argv[2])

def Needleman_Algo(seq1, seq2):
    row = len(seq1) + 1
    col = len(seq2) + 1

    dp = [[0]*col for _ in range(row)]

    for i in range(row):
        dp[i][0] = i * gap
    for j in range(col):
        dp[0][j] = j * gap

    for i in range(1, row):
        for j in range(1, col):
            if seq1[i-1] == seq2[j-1]:
                score = match
            else:
                score = mismatch

            dp[i][j] = max(
                dp[i-1][j-1] + score,
                dp[i-1][j] + gap,
                dp[i][j-1] + gap
            )

    # TRACEBACK
    align1 = ""
    align2 = ""
    i = row - 1
    j = col - 1

    while i > 0 and j > 0:
        if seq1[i-1] == seq2[j-1]:
            score = match
        else:
            score = mismatch

        if dp[i][j] == dp[i-1][j-1] + score:
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1

        elif dp[i][j] == dp[i-1][j] + gap:
            align1 = seq1[i-1] + align1
            align2 = "-" + align2
            i -= 1

        else:
            align1 = "-" + align1
            align2 = seq2[j-1] + align2
            j -= 1

    # (optional but better) handle remaining
    while i > 0:
        align1 = seq1[i-1] + align1
        align2 = "-" + align2
        i -= 1

    while j > 0:
        align1 = "-" + align1
        align2 = seq2[j-1] + align2
        j -= 1

    # Alignment line
    line = ""
    for l in range(len(align1)):
        if align1[l] == align2[l]:
            line += "|"
        elif align1[l] == "-" or align2[l] == "-":
            line += " "
        else:
            line += "."

    print("\nNeedleman-Wunsch Algorithm\n")
    print("Sequence 1:", seq1)
    print("Sequence 2:", seq2)
    print("\nAlignment:\n")
    print(align1)
    print(line)
    print(align2)
    print("\nFinal score:", dp[row-1][col-1])


Needleman_Algo(seq1, seq2)
