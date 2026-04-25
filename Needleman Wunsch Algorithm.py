# Needleman-Wunsch Global Pairwise DNA Alignment

#Scoring Scheme
match = 1
mismatch = -1
gap = -2
seq1 = input("Enter your 1st sequence here:")
seq2 = input("Enter your 2nd Sequence here:")
def Needleman_Algo(seq1,seq2):
    row = len(seq1) + 1
    col = len(seq2) + 1

    dp = [[0]* col for _ in range(row)]
    for i in range(row):
        dp[i][0] = i * gap
    for j in range(col):
        dp[0][j] = j * gap

    for i in range(1,row):
        for j in range(1,col):
            if seq1[i-1] == seq2[j-1]:
                score = match
            else:
                score = mismatch
                
            dp[i][j] = max(dp[i-1][j-1] + score,
                               dp[i-1][j] + gap,
                               dp[i][j-1] + gap)
        
    #TRACE BACK
    align1 = ""
    align2 = ""
    i = row - 1
    j = col- 1
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
#Building Alignment line
    line = ""

    for l in range(len(align1)):
         
        if align1[l] == align2[l]:
            line = line + "|"
        elif align1[l] == "-" or align2[l] == "-": #Base pair matched
            line = line + " " # Gap
        else:
            line = line + "." # Mismatched


    print("             Needleman-Wunsch Algorithm              ")
    print("Sequence 1 :",seq1)
    print("Sequence 2 :",seq2)
    print("Alignment of Sequences:")
    print(align1)
    print(line)
    print(align2)       
    print("Final score :",dp[row -1][col-1])
Needleman_Algo(seq1,seq2)
