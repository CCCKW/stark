import pandas as pd

def generate_barcode(index, pcr_file, tn5_file,out):
    # 定义序列
    seq1 = "AATGATACGGCGACCACCGAGATCT"
    seq2 = "TCGTCGGCAGCGTC"
    seq3 = "AGATGTGTATAAGAGACAG"
    
    # 读取PCR和TN5的条形码文件
    PCR = pd.read_csv(pcr_file, sep='\t', header=None, index_col=0)
    TN5 = pd.read_csv(tn5_file, sep='\t', header=None, index_col=0)
    
    with open(out, 'w') as output:
        for line in index:
            a = line.strip().split()
            if a[0].isdigit():
                PCR_BC = int(a[0])
                BC_start = int(a[1]) if len(a) > 1 else 1
                BC_end = int(a[2]) if len(a) > 2 else 24
                PCR_Name = f'P{PCR_BC}'

                PCR_Seq = PCR.loc[PCR_BC, 1]
                for TN5_BC in range(BC_start, BC_end + 1):
                    TN5_Seq = TN5.loc[TN5_BC, 1]
                    Name = f'{PCR_Name}B{TN5_BC}'
                    BARCODE = f'{seq1}{PCR_Seq}{seq2}{TN5_Seq}{seq3}'
                    
                    # print(f"Name: {Name}")
                    # print(f"PCR_Seq: {PCR_Seq}")
                    # print(f"TN5_Seq: {TN5_Seq}")
                    # print("==============================================================")

                    output.write(f">{Name}\n")
                    output.write(f"{BARCODE}\n")
            else:
                print("Barcode is not an integer!")
                break