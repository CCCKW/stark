from sc3dg.model.nuc_dynamics import calc_genome_structure
from time import time
import os
import click
def generate_3d(pair_file_path:str,
                output:str,
                num_models:int,
                iter_steps:int,
                iter_res:list,
                support_chrom:list):

    # Number of alternative conformations to generate from repeat calculations
    # with different random starting coordinates
    num_models = num_models

    # Parameters to setup restraints and starting coords
    general_calc_params = {'dist_power_law': -0.33,
                           'contact_dist_lower': 0.8, 'contact_dist_upper': 1.2,
                           'backbone_dist_lower': 0.1, 'backbone_dist_upper': 1.1,
                           'random_seed': int(time()), 'random_radius': 10.0}

    # Annealing & dyamics parameters: the same for all stages
    # (this is cautious, but not an absolute requirement)
    anneal_params = {'temp_start': 5000.0, 'temp_end': 10.0, 'temp_steps': 100,
                     'dynamics_steps': iter_steps, 'time_step': 0.001}
                     # ynamics_step 精细度
    # Hierarchical scale protocol: calculations will initially use 8 Mb particles
    # but deminish to 100 kb at the end. The whole annealing protocol (hot to cold)
    # will be run at each size, but subsequent stages will start from the previous
    # structure
    # particle_sizes = [8e6, 4e6, 2e6, 4e5, 2e5, 1e5]
    particle_sizes =  iter_res
    # Contacts must be clustered with another within this separation threshold
    # (at both ends) to be considered supported, i.e. not isolated
    # This removes noise contacts
    isolation_threshold = 2e6

    Target_Chrom = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                    'chr13', 'chr14', 'chr15', 'chr16',
                    'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']
    if os.path.exists(output):
        import sys
        sys.exit("Output file already exists: %s" %output)



    # Actually run the calculation with the specified parameters and input
    # The below function will automatically create the appropriate distance restraints
    calc_genome_structure(pair_file_path, output, general_calc_params, anneal_params,
                          particle_sizes, num_models, isolation_threshold, support_chrom=support_chrom, out_format='pdb')


@click.command('model', short_help='3d model and generate pdb file.')
@click.option('--pair', required=True, type=click.Path(exists=True), help='Path to the pair file.')
@click.option('--output', required=True, type=click.Path(), help='Path to the output pdb file.')
@click.option('--num_models', default=5, type=int, help='Number of models to generate.')
@click.option('--iter_steps', default=10, type=int, help='Number of iteration steps.')
@click.option('--iter_res', default="8e6,4e6,2e6,4e5,2e5,1e5", help='Comma-separated list of iteration resolutions.')
@click.option('--support_chrom', default=None, help='Comma-separated list of supported chromosomes.')
def model(pair, output, num_models, iter_steps, iter_res, support_chrom):
    """Model command to run the specified model with the config file."""
    iter_res_list = [float(res) for res in iter_res.split(',')]
    generate_3d(pair, output, num_models, iter_steps, iter_res_list, support_chrom)