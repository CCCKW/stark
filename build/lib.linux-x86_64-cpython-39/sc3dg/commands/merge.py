from sc3dg.analysis.api import MMCooler
import click

import cooler
import os

'''
    Traversing the directory, and mergeing the mcool files in the directory
'''


@click.command('merge', short_help='merge mcool files')
@click.option('--mcool', help='mcool directory path', required=True, default=None)
@click.option('--output', help='output path', required=True, default=None)
@click.option('--resolution', help='resolution', required=True, type=str)
@click.option('--num', help='top number', required=False, default=None, type=int)
@click.option('--sort', help='sort', required=False, default=True, type=bool)
@click.option('--n_jobs', help='n_jobs', required=False, default=10, type=int)

def  merge(mcool,
                           output,
                           resolution,
                           num,
                           sort=True,
                           n_jobs=10):
    '''
    Traversing the directory, and mergeing the mcool files in the directory,
    if num is larger than the number of the mcool files, then merge all the mcool files
    Parameters
    ----------
    mcool: str
        The directory path of the mcool files, traverse the directory and merge the mcool files
    output: str
        The output file path, which should be a file name
    resolution: list, like [10000, 40000, 100000, 500000, 1000000]
        The resolution of the mcool files

    num: int
        The number of the top level mcool files
    sort: bool
        Whether to sort the mcool files, default is True, merge Top level mcool files
        if false, merge the bottom level mcool files
    Returns
    -------
        None
    '''
    data = MMCooler(directory=mcool,
                    resolution=resolution,
                    merge=False,
                    n_jobs=n_jobs
                    )

    data.sorted_cell(num=num, reverse=sort)

    # save
    cooler.merge_coolers(input_uris=[x.uri for x in data.clr_obj[0]],
                         output_uri=output + '.cool',
                         mergebuf=1000000)

    cooler.zoomify_cooler(base_uris=output, outfile=output,
                          resolutions=resolution,
                          chunksize=1000000, nproc=16)
    os.remove(output + '.cool')








