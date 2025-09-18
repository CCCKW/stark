import click
import importlib

class LazyGroup(click.Group):
    def get_command(self, ctx, cmd_name):
        try:
            module = importlib.import_module(f'sc3dg.commands.{cmd_name}')
            return getattr(module, cmd_name)
        except (ImportError, AttributeError) as e:
            click.echo(f"Error loading command {cmd_name}: {e}")
            return None

@click.group(cls=LazyGroup)
def cli():
    '''
        stark is a software package for single-cell three-dimensional genome sequencing data(sc3dg-seq) preprocessing and analysis.

        Usage: stark [OPTIONS] COMMAND [ARGS]...

        Options:
            --help  Show this message and exit.

        Commands: 
        
            count           the main pipeline of Hi-C \n
            model           3D genome model ,generate the pdb file \n 
            index           generate the index \n
            emptycells      emptycells for QC \n
            impute          impute the Hi-C data \n
            ssce            cal Spatial Structure Capture Efficiency \n
            merge           merge the Hi-C data \n
            accum           accumulative analysis \n
            nDS             calculate the nDS \n
            gini            calculate the gini \n
            loop            detect the loop \n
            clustring       using higashi,fasthigshi,deepNanoHIiC to clustring the single cell \n
            

    '''
    pass

commands = [
    'count', 'model', 'impute', 'ssce', 'merge', 'accum',
    'nDS', 'gini', 'loop', 'index', 'emptycells', 'clustring'
]

if __name__ == '__main__':
    cli()
