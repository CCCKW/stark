install the STARK
========================

We suggest creating a new Python environment to install STARK before installation. STARK is based on Python ≥ 3.9.
There are two ways to install stark:

Command
++++++++++++++++++++
This installation method is based on github.

step1.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

create a new environment

.. code-block:: bash

    conda create -n sc3dg python=3.9
    conda activate sc3dg
    git clone https://github.com/CCCKW/stark.git
    cd stark
    pip install -r requirements.txt
    python setup.py install

It 's recommended to install the requirements.txt first, then install the STARK package.

Attention: it take a little time to install the requirements and the STARK package. We encapsulate most of the packages with complex dependencies. 
Also, you can use the nohup command to install the package in the background.
Besides, if you have any problems with the installation, especially "pip install -r requirements.txt", the installation of pysam may fail. 
You can try to install it with the following command:

.. code-block:: bash

    conda install -c bioconda pysam

After the installation of pysam, you can continue to install the requirements.



step2.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

prepare index and enzyme bed

STARK has its own format for reading reference genomes and enzyme bed files. For convenience, we provide an automated program for building indexes and enzyme bed files:

.. code-block:: shell

    stark index -g hg38 -a bwa -p ./save/path1 -e MboI,DpnII,BglII,HpyCH4V
    stark index -g hg38 -a bowtie2 -p ./save/path2 -e MboI,DpnII,BglII,HpyCH4V

Parameter description:

- **-g**: The genome version you want to assemble, currently supports hg38, hg19, mm10, mm9.
- **-a**: The software used for assembly, different software produces different index files.
- **-p**: Path where the results are stored.
- **-e**: The required enzyme bed file, can be single or multiple inputs (different genomes correspond to different enzyme files).
