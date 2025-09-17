import os
from joblib import Parallel, delayed
import logging
import sys

path = sys.argv[1]
to_path = sys.argv[2]
os.makedirs(to_path, exist_ok=True)
# hg38 or mm10
org = sys.argv[3]
resolution = sys.argv[4]
 
print(org, resolution)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logging.info('Starting the cload process...')
to_do = []
for pair in os.listdir(path):
    if pair.endswith('.pairs.gz'):
        to_do.append(pair)

logging.info(f'Found {len(to_do)} pairs to process.')

def cload(pair):
    if os.path.exists(to_path + '/' + pair.replace('.pairs.gz', '.cool')):
        logging.info(f'Skipping {pair}, already processed.')
        return
    cmd = 'cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 '
    
    cmd += f'{org}:' + str(resolution) + ' ' 

    cmd += ' ' + path + '/' + pair
    cmd += ' ' + to_path + '/' + pair.replace('.pairs.gz', '.cool')
    logging.info(f'Running command: {cmd}')
    print(f'Processing {pair}...')
    os.system(cmd)
    
Parallel(n_jobs=10)(delayed(cload)(pair) for pair in to_do)