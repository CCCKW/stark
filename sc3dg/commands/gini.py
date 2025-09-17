import click
from sc3dg.analysis.api import MMCooler
import pandas as pd


@click.command('gini', short_help='calculate gini')
@click.option('--mcool',  help='mcool path',required=True, default=None)
@click.option('--output',  help='output path',required=True, default=None)
@click.option('--resolution', help='resolution',required=True,  type=str)
@click.option('--n_jobs', help='n_jobs',default=1,  type=int)
def gini(mcool, output, resolution,n_jobs):

    resolutions = list(resolution.split(','))
    resolutions = [int(x) for x in resolutions]
    print(mcool, output, resolutions, n_jobs)
    
    ob1 = MMCooler(directory=mcool,
                   resolution=resolution.split(','),
                   describe = 'gini',
                   merge=False
                   )
    gini = ob1.get_giniQC(n_jobs=n_jobs)

    gini.to_csv(output, sep='\t')