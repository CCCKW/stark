import click
import time
from sc3dg.utils.download import make_index as download_make_index

@click.command('index', short_help='create the index for the count')
@click.option('-g', '--genome', type=click.Choice(['hg38', 'hg19', 'mm10', 'mm9']), required=True,
              help='Choose from hg38, hg19, mm10, mm9')
@click.option('-a', '--aligner', type=click.Choice(['bowtie2', 'bwa']), required=True,
              help='Choose from bowtie2, bwa')
@click.option('-p', '--path', required=True, help='Path to save index')
def index(genome, aligner, path):
    """Create an index for the specified genome and aligner."""
    click.echo(f"Creating index for {genome} using {aligner}...")
    download_make_index(path, genome, aligner)
    click.echo("Index creation completed.")

if __name__ == '__main__':
    index()