Clustering
=========================================

STARK has built-in clustering algorithms such as Higashi, Fast-higashi, deepnanoHi-C, and schicluster. These algorithms are all based on clustering from pairs.gz files.


To run clustering, use the following command:

.. code-block:: shell

    stark clustering --method <method_name> \
        --config <config_file> \


example:

.. code-block:: shell

    stark clustering --method higashi \
        --config config.json

The configuration file should be in JSON format and include the following fields:

Attention:

if you want to running Higashi, Fast-Higashi or DeepnanoHi-C, you need to install torch and cuda first.

You find the installation instructions here: https://pytorch.org/get-started/locally/




Higashi:

.. code-block:: json
    
    {
        "data_dir": "/path/to/final_dir",
        "structured": true,
        "input_format": "higashi_v2",
        "header_included": true,
        "temp_dir": "/path/to/final_dir",
        "genome_reference_path": "g38.fa.chrom.sizes",
        "cytoband_path": "cytoBand_hg38.txt",
        "chrom_list": ["chr1", "chr2", "chr3", "chr4", "chr5",
                        "chr6", "chr7", "chr8", "chr9", "chr10",
                        "chr11", "chr12", "chr13", "chr14", "chr15",
                        "chr16", "chr17", "chr18", "chr19", "chr20",
                        "chr21", "chr22"],
        "resolution": 50000,
        "resolution_cell": 50000,
        "resolution_fh": [50000],
        "embedding_name": "test",
        "minimum_distance": 50000,
        "maximum_distance": -1,
        "local_transfer_range": 0,
        "loss_mode": "zinb",
        "dimensions": 128,
        "impute_list": ["chr1", "chr2", "chr3", "chr4", "chr5",
                        "chr6", "chr7", "chr8", "chr9", "chr10",
                        "chr11", "chr12", "chr13", "chr14", "chr15",
                        "chr16", "chr17", "chr18", "chr19", "chr20",
                        "chr21", "chr22"],
        "neighbor_num": 5,
        "cpu_num": 10,
        "gpu_num": 8,
        "embedding_epoch": 60,
        "correct_be_impute": true
        }

Fast-higashi:

.. code-block:: json
    
   {
  "data_dir": "/path/to/final_dir",
  "structured": true,
  "input_format": "higashi_v2",
  "header_included": true,
  "temp_dir":  "/path/to/final_dir",
  "genome_reference_path": "hg38.fa.chrom.sizes",
  "cytoband_path": "cytoBand_hg38.txt",
  "chrom_list": ["chr1", "chr2", "chr3", "chr4", "chr5",
                 "chr6", "chr7", "chr8", "chr9", "chr10",
                 "chr11", "chr12", "chr13", "chr14", "chr15",
                 "chr16", "chr17", "chr18", "chr19", "chr20",
                 "chr21", "chr22"],
  "resolution": 50000,
  "resolution_cell": 50000,
  "resolution_fh": [50000],
  "embedding_name": "test",
  "minimum_distance": 50000,
  "maximum_distance": -1,
  "local_transfer_range": 0,
  "loss_mode": "zinb",
  "dimensions": 128,
  "impute_list": ["chr1", "chr2", "chr3", "chr4", "chr5",
                  "chr6", "chr7", "chr8", "chr9", "chr10",
                  "chr11", "chr12", "chr13", "chr14", "chr15",
                  "chr16", "chr17", "chr18", "chr19", "chr20",
                  "chr21", "chr22"],
  "neighbor_num": 5,
  "cpu_num": 10,
  "gpu_num": 8,
  "embedding_epoch": 60,
  "correct_be_impute": true
}


DeepnanoHi-C:

.. code-block:: json
    
    {
    "data_dir": "/path/to/final_dir",
    "temp_dir": ""/path/to/final_dir"",
    "structured": true,
    "input_format": "higashi_v2",
    "header_included": true,
    "genome_reference_path": "hg38.fa.chrom.sizes",
    "cytoband_path": "cytoBand_hg38.txt",
    "chrom_list": ["chr1", "chr2", "chr3", "chr4", "chr5",
                   "chr6", "chr7", "chr8", "chr9", "chr10",
                   "chr11", "chr12", "chr13", "chr14", "chr15",
                   "chr16", "chr17", "chr18", "chr19", "chr20",
                   "chr21", "chr22"],
    "resolution": 500000,
    "resolution_cell": 500000,
    "resolution_fh": [500000],
    "embedding_name": "test",
    "minimum_distance": 50000,
    "maximum_distance": -1,
    "local_transfer_range": 0,
    "loss_mode": "zinb",
    "dimensions": 128,
    "impute_list": ["chr1", "chr2", "chr3", "chr4", "chr5",
                    "chr6", "chr7", "chr8", "chr9", "chr10",
                    "chr11", "chr12", "chr13", "chr14", "chr15",
                    "chr16", "chr17", "chr18", "chr19", "chr20",
                    "chr21", "chr22"],
    "neighbor_num": 5,
    "cpu_num": 10,
    "gpu_num": 8,
    "embedding_epoch": 60,
    "correct_be_impute": true
  }


Schicluster:


.. code-block:: json
    
  {
    "data_dir": "/path/to/final_dir",
    "temp_dir": "/path/to/final_dir",
    "filelist": "/path/to/final_dir/filelist.txt",
    "genome_reference_path": "hg38.fa.chrom.sizes",
    "chrom_list": ["chr1", "chr2", "chr3", "chr4", "chr5",
                   "chr6", "chr7", "chr8", "chr9", "chr10",
                   "chr11", "chr12", "chr13", "chr14", "chr15",
                   "chr16", "chr17", "chr18", "chr19", "chr20",
                   "chr21", "chr22"],
    "impute_list": ["chr1", "chr2", "chr3", "chr4", "chr5",
                    "chr6", "chr7", "chr8", "chr9", "chr10",
                    "chr11", "chr12", "chr13", "chr14", "chr15",
                    "chr16", "chr17", "chr18", "chr19", "chr20",
                    "chr21", "chr22"],
    "dist": 2500,
    "res": 50000,
    "n_jobs": 10
  }


You need to give the filelist.txt file in the target data_dir directory, which contains the pairs.gz files of all the cells you need to cluster, like:

.. code-block:: text

    /path/to/final_dir/[Even2Bo10][Odd2Bo69][DPM6bot1].pairs.gz
    /path/to/final_dir/[Even2Bo11][Odd2Bo19][DPM6bot31].pairs.gz
    /path/to/final_dir/[Even2Bo11][Odd2Bo1][DPM6bot75].pairs.gz