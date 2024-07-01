calculate the GiniQC
=========================

Single-cell Hi-C (scHi-C) allows the study of cell-to-cell variability in chromatin structure and dynamics. However, the high level of noise inherent in current scHi-C protocols necessitates careful assessment of data quality before biological conclusions can be drawn. Here, we present GiniQC, which quantifies unevenness in the distribution of inter-chromosomal reads in the scHi-C contact matrix to measure the level of noise. Our examples show the utility of GiniQC in assessing the quality of scHi-C data as a complement to existing quality control measures.

To calculate the GiniQC, you need to provide the path to the mcool file and the output path to save the results.

Usage

.. code-block:: shell

    stark gini --mcool/absolute/path/to/data/gini \
        --output/absolute/path/to/result/gini/test_gini.csv \
        --resolution 10000,40000,100000 