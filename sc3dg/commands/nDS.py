import click
from sc3dg.analysis.api import MMCooler
import pandas as pd

@click.command('nDS', short_help='nDS calculation')
@click.option('--mcool',  help='mcool path',required=True, default=None)
@click.option('--output',  help='output path',required=True, default=None)
@click.option('--resolution', help='resolution',required=True,  type=int)
@click.option('--mode', help='mode',required=True,type=click.Choice(['terr', 'tad', 'ab', 'cen']))
@click.option('--top_num', help='top number',required=False, default=None, type=int)
@click.option('--describe', help='describe',required=False, default='nDS')
@click.option('--epoch', help='epoch',required=False, default=10, type=int)
@click.option('--genome', help='fa genome',required=False, default='mm10')

def nDS(mcool, output, resolution, mode, top_num, describe, epoch, genome):

    '''

    Parameters
    ----------
    mcool:str
        mcool file path
    output:
        output file path
    resolution
    mode:str
        mode = ['terr', 'tad', 'ab', 'cen']
    top_num:int
        top number
    describe:str
        describe, for save file name
    epoch:int
        epoch
    Returns
    -------

    '''


    ob1 = MMCooler(directory=mcool,
                   resolution=[resolution],
                   describe = mode if not describe else describe,
                   merge=False
                   )

    if top_num:
        ob1.sorted_cell(num=top_num)

    

    if mode == 'terr':
        terr = ob1.get_nDS_chrom_terr(epoch=epoch,res=resolution)[0]
        idx = ob1.mcool_names
        terr.index = idx
        terr.to_csv(output, sep='\t')
    elif mode == 'tad':
        tad = ob1.get_TAD_nDS(epoch=epoch, res=resolution)
        df = pd.DataFrame({'tad':tad})
        idx = ob1.mcool_names
        df.index = idx
        df.to_csv(output, sep='\t')
    elif mode == 'ab':
        if not genome:
            raise ValueError('genome is required')
        ab = ob1.get_AB_compartment(res=resolution, epoch =epoch,genome=genome)
        idx = ob1.mcool_names
       
        df = pd.DataFrame({'ab':ab})
        df.index = idx
        df.to_csv(output, sep='\t')
    elif mode == 'cen':
        cen = ob1.get_Centromere_nDs(epoch=epoch, res=resolution)
        idx = ob1.mcool_names
        cen.index = idx
        cen.to_csv(output, sep='\t')





