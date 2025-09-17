import json
import os
import pandas as pd
import cooler

def generatematrix_cell_mcool(
    config,
    infile,
    ):
    
    config = json.load(open(config))
    chrom = config['chrom_list']
    outdir = config['temp_dir'] + '/cell_matrix'
    res = int(config['res'])
    dist = int(config['dist'])
    
    cell = os.path.basename(infile).split('.')[0]
    
    chrom_split = chrom.copy() 
    for c in chrom_split:
        os.makedirs(f'{outdir}/{c}/', exist_ok=True)
    
    # 检查文件类型并读取数据
    if infile.endswith(('.cool', '.mcool')):
        # 处理cool/mcool文件
        if infile.endswith('.mcool'):
            # 对于mcool文件，需要指定分辨率
            clr = cooler.Cooler(f'{infile}::/resolutions/{res}')
        else:
            clr = cooler.Cooler(infile)
        
        # 提取接触数据
        data_list = []
        for c in chrom_split:
            if c not in clr.chromnames:
                continue
            
            # 获取染色体内的接触矩阵
            matrix = clr.matrix(balance=False).fetch(c)
            bins = clr.bins().fetch(c)
            
            # 转换为三元组格式 (pos1, pos2, count)
            for i in range(len(bins)):
                for j in range(i+1, len(bins)):
                    if matrix[i, j] > 0:
                        pos1 = bins.iloc[i]['start'] // res
                        pos2 = bins.iloc[j]['start'] // res
                        if abs(pos2 - pos1) * res > dist:
                            data_list.append([c, pos1, pos2, matrix[i, j]])
        
        if data_list:
            data = pd.DataFrame(data_list, columns=[1, 2, 4, 'count'])
        else:
            return
            
    else:
        # 原有的pairs文件处理逻辑
        chr1, pos1, chr2, pos2 = 1, 2, 3, 4
        data = pd.read_csv(infile, sep='\t', header=None, comment='#')
        
        data.loc[data[pos1]>data[pos2], [pos1, pos2]] = data.loc[data[pos1]>data[pos2], [pos2, pos1]]
        data = data[(data[chr1]==data[chr2]) & (data[chr1].isin(chrom)) & (data[pos2] - data[pos1] > dist)]
        data[[pos1, pos2]] = data[[pos1,pos2]] // res
        data = data[(data[chr1]==data[chr2]) & (data[chr1].isin(chrom_split)) & (data[pos1]!=data[pos2])]
        # count contacts at each pixel
        data = data.groupby(by=[chr1, pos1, pos2])[chr2].count().reset_index()
    
    # 输出结果
    for c, tmp in data.groupby(by=1):
        if os.path.exists(f'{outdir}/{c}/{cell}_{c}.txt'):
            continue
        
        if infile.endswith(('.cool', '.mcool')):
            tmp[[2, 4, 'count']].astype(int).to_csv(f'{outdir}/{c}/{cell}_{c}.txt', sep='\t', header=False, index=False)
        else:
            tmp[[2, 4, 3]].astype(int).to_csv(f'{outdir}/{c}/{cell}_{c}.txt', sep='\t', header=False, index=False)
    
    return

if __name__ == '__main__':
    generatematrix_cell(
        config='config.JSON',
        infile='example.cool'  # 或 'example.mcool'
    )