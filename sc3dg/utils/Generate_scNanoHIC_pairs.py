import pandas as pd
import sys


def WriteCellPairs(ChromSize, SamFile, GenomeInfo, ResultDir):

    fp = open(ResultDir + "/" + str(SamFile).split("/")[-1].split(".")[0]+".pairs", "w")
    fp.write("## pairs format v1.0.0\n#shape: whole matrix\n#genome_assembly: " + GenomeInfo + "\n")
    for line in open(ChromSize):
        fp.write("#chromsize: " + line.split("\t")[0] + " " + line.split("\t")[1])
    

    Anchor_longreads = ''
    Anchor_longreads_pairsInfo = []
    flag=0
    for line in open(SamFile):

        ####Deal Sam file head
        if "@" == line[0]:
            fp.write("#samheader: " + line)

        else:
            if flag==0:
                fp.write("#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type mapq1 mapq2\n")
                flag+=1
            ReadComponents = line.split("\t")

            ####判断正负链
            if bool(int(ReadComponents[1]) & 0x1000):
                ReadComponents[1] = '+'
            else:
                ReadComponents[1] = '-'

            current_longreads = ReadComponents[0].split(":")[0]

            if current_longreads == Anchor_longreads:
                Anchor_longreads_pairsInfo.append([ReadComponents[0], ReadComponents[2], ReadComponents[3], ReadComponents[1], ReadComponents[4]])

            else:
                for i in range(len(Anchor_longreads_pairsInfo)):
                    for j in range(0,i):
                        fp.write(Anchor_longreads_pairsInfo[i][0] + "_" + Anchor_longreads_pairsInfo[j][0] + "\t" +
                                 Anchor_longreads_pairsInfo[i][1] + "\t" + Anchor_longreads_pairsInfo[i][2] + "\t" +
                                 Anchor_longreads_pairsInfo[j][1] + "\t" + Anchor_longreads_pairsInfo[j][2] + "\t" +
                                 Anchor_longreads_pairsInfo[i][3] + "\t" + Anchor_longreads_pairsInfo[j][3] + "\tUU\t" +
                                 Anchor_longreads_pairsInfo[i][4] + "\t" + Anchor_longreads_pairsInfo[j][4] + "\n")

                Anchor_longreads = current_longreads
                Anchor_longreads_pairsInfo = []
                Anchor_longreads_pairsInfo.append([ReadComponents[0], ReadComponents[2], ReadComponents[3], ReadComponents[1], ReadComponents[4]])

    for i in range(len(Anchor_longreads_pairsInfo)):
        for j in range(0,i):
            fp.write(Anchor_longreads_pairsInfo[i][0] + "_" + Anchor_longreads_pairsInfo[j][0] + "\t" +
                     Anchor_longreads_pairsInfo[i][1] + "\t" + Anchor_longreads_pairsInfo[i][2] + "\t" +
                     Anchor_longreads_pairsInfo[j][1] + "\t" + Anchor_longreads_pairsInfo[j][2] + "\t" +
                     Anchor_longreads_pairsInfo[i][3] + "\t" + Anchor_longreads_pairsInfo[j][3] + "\tUU\t" +
                     Anchor_longreads_pairsInfo[i][4] + "\t" + Anchor_longreads_pairsInfo[j][4] + "\n")

    fp.close()









