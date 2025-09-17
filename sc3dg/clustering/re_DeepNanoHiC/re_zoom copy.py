
import cooler
import numpy as np
import pandas as pd
import sys
import os
import matplotlib.pyplot as plt
import sys

import shutil
from joblib import Parallel, delayed
path = sys.argv[1]
outpath = sys.argv[2]
os.makedirs(outpath, exist_ok=True)

def func(idx,full_path, name, outpath):
   

  
    # cooler.zoomify_cooler(base_uris=full_path +  '::resolutions/10000',
    #                 outfile= outpath + '/' + name + '.mcool',
    #                     resolutions=[10000,50000,100000,500000,1000000],
    #                         chunksize=1000000, nproc=32)
    try:
        cooler.zoomify_cooler(base_uris=full_path ,
                        outfile= outpath + '/' + name + '.mcool',
                            resolutions=[10000,50000,100000,500000,1000000],
                                chunksize=1000000, nproc=32)
    except:
        print('Error processing:', full_path)
        return
    print('Finished:', idx)
# scan the directory
name = []
dirs = []
for value in os.listdir(path):
    if value.endswith('.cool'):
        name.append(value.split('.')[0])
        dirs.append(path + '/' + value)

print(len(name))
task = [[i, dirs[i], name[i], outpath] for i in range(len(name))]
os.makedirs(outpath, exist_ok=True)
# Parallelize the function
Parallel(n_jobs=10)(delayed(func)(task[i][0], task[i][1], task[i][2],task[i][3]) for i in range(len(task)))
print('all done')





