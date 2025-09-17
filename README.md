# STARK tutorial

STARK is a  software for processing various types of single-cell Hi-C sequencing data, which currently includes the  single-cell Hi-C technologies as follows:

- [scHiC](https://www.nature.com/articles/nature12593)
- [scHi-C+](https://www.nature.com/articles/nature23001)
- [Dip-C](https://www.science.org/doi/10.1126/science.aat5641)
- [HiRES](https://www.science.org/doi/10.1126/science.adg3797)
- [sn-m3C](https://www.nature.com/articles/s41592-019-0547-z)
- [scSPRITE](https://www.nature.com/articles/s41587-021-00998-1)
- [sciHi-C](https://www.nature.com/articles/nmeth.4155)
- [snHi-C](https://www.nature.com/articles/nature21711)
- [snHi-C+](https://data.4dnucleome.org/experiment-set-replicates/4DNESF829JOW/)
- [scNanoHi-C](https://www.nature.com/articles/s41592-023-01978-w)
- [scMethyl](https://www.nature.com/articles/s41592-019-0502-z)
- [LiMAC](https://www.nature.com/articles/s41592-024-02239-0)
- [GAGE-seq](https://www.nature.com/articles/s41588-024-01745-3)
- [Droplet](https://www.nature.com/articles/s41587-024-02447-1)
- [Paired](https://www.nature.com/articles/s41587-024-02447-1)

At the same time, we also provide a generalized way to process the sc3DG for details please see the tutorial.

To get started：
- Install STARK
- Read the documentation and see the Jupyter Notebook walkthrough.
- test data can be download form [Zenodo](https://zenodo.org/records/12598215/files/sc3dg.tar.gz)
- Many more single-cell three dimension genome sequencing data are available on our website

# Installation


We suggest creating a new python environment to install STARK before installation and STARK is based on python≥3.9.

stark can only be installed on the linux.

```shell

    conda create -n sc3dg python=3.9
    conda activate sc3dg
    git clone https://github.com/CCCKW/stark.git
    cd stark-main
    pip install -r requirements.txt
    python setup.py install

```

Attention: if you have any problem with installation of pysam, try to install the pysam as the follow command before run the 'pip install -r requirements.txt':

```shell

    conda install -c bioconda pysam

```

For usage and more detailed information of stark , you can visit the [tutorial](https://sc3dg.readthedocs.io/)


Also, you can visit the [scNucleome](http://scnucleome.com/) for more detailed information.

