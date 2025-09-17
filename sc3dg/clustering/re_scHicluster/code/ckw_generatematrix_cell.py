import json
import os
import pandas as pd


def generatematrix_cell(
    config,
    infile,
    
    ):
    
    config = json.load(open(config))
    chrom = config['chrom_list']
    outdir = config['temp_dir'] + '/cell_matrix'
    chr1, pos1, chr2, pos2 = 1, 2, 3, 4
    res = int(config['res'])
    dist=int(config['dist'])
    
    cell = os.path.basename(infile).split('.')[0]
    
    
    
    
    chrom_split = chrom.copy() 
    for c in chrom_split:
       
        os.makedirs(f'{outdir}/{c}/', exist_ok=True)
    
    data = pd.read_csv(infile, sep='\t', header=None,comment='#')
    
    data.loc[data[pos1]>data[pos2], [pos1, pos2]] = data.loc[data[pos1]>data[pos2], [pos2, pos1]]
    data = data[(data[chr1]==data[chr2]) & (data[chr1].isin(chrom)) & (data[pos2] - data[pos1] > dist)]
    data[[pos1, pos2]] = data[[pos1,pos2]] // res
    data = data[(data[chr1]==data[chr2]) & (data[chr1].isin(chrom_split)) & (data[pos1]!=data[pos2])]
    # count contacts at each pixel
    data = data.groupby(by=[chr1, pos1, pos2])[chr2].count().reset_index()
    for c, tmp in data.groupby(by=chr1):
            if os.path.exists(f'{outdir}/{c}/{cell}_{c}.txt'):
                continue

            tmp[[pos1, pos2, chr2]].astype(int).to_csv(f'{outdir}/{c}/{cell}_{c}.txt', sep='\t', header=False, index=False)

            
    
    return
    
    
if __name__ == '__main__':
    generatematrix_cell(
        config='config.JSON',
        infile='GSM7682104_BJ.pairs.gz'
    )