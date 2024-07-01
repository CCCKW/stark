# STARK tutorial

STARK is a  software for processing various types of single-cell Hi-C sequencing data, which currently includes the  single-cell Hi-C technologies as follows:

- scHiC
- scHi-C+
- Dip-C
- HiRES
- sn-m3C
- scSPRITE
- sciHi-C
- snHi-C
- snHi-C
- scNanoHi-C
- scMethyl
- scCARE-seq

At the same time, we also provide a generalized way to process the sc3DG for details please see the tutorial.

To get started：
- Install STARK
- Read the documentation and see the Jupyter Notebook walkthrough.
- test data can be download form [Zenodo](https://zenodo.org/records/12598215/files/sc3dg.tar.gz)
- Many more single-cell Hi-C data are available on our website

# Installation


We suggest creating a new python environment to install STARK before installation and STARK is based on python≥3.9

```shell

    conda create -n sc3dg python>=3.9
    conda activate sc3dg
    git clone https://github.com/CCCKW/stark.git
    cd stark-main
    pip install -r requirements.txt
    python setup.py install

```

For usage and more detailed information of stark , you can visit the [tutorial](https://sc3dg.readthedocs.io/en/latest/overview.html)


Also, you can visit the [scNucleome](https://www.baidu.com) for more detailed information.

