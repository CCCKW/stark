import click
import pandas as pd


@click.command('clustering', short_help='clustering')
@click.option('--method',  help='mcool path',required=True, default='higashi', type=click.Choice(['higashi', 'fasthigashi', 'deepNanoHiC','schicluster'], case_sensitive=False))
@click.option('--config',  help='config',required=True, default=None)


def clustering(method, config):

    if method not in ['higashi', 'fasthigashi', 'deepNanoHiC']:
        raise ValueError("Method must be one of 'higashi', 'fasthigashi', or 'deepNanoHiC'")

    if method == 'higashi':
        from sc3dg.clustering.re_higashi.higashi_pipe import workflow as higashi_clustering
        
        higashi_clustering(config)
    elif method == 'fasthigashi':
        from sc3dg.clustering.fasthigashi.fasthigashi_pipe import fasthigashi_clustering
        fasthigashi_clustering(config, output)
    elif method == 'deepNanoHiC':
        from sc3dg.clustering.deepNanoHiC.deepNanoHiC_pipe import deepNanoHiC_clustering
        deepNanoHiC_clustering(config, output)
    elif method == 'schicluster':
        from sc3dg.clustering.schiCluster.schiCluster_pipe import schiCluster_clustering
        schiCluster_clustering(config, output)