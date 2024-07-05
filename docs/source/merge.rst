Merge the topN/bottomN cell
=================================
.. raw:: html

   <style type="text/css">
   body {
       text-align: justify;
   }
   </style>

Given the sparsity of sc3DG-seq data, an integrated analysis of data from multiple
cells is often beneficial.
STARK facilitates this by enabling the aggregation of data from multiple single cells
for subsequent analysis. The aggregation is mathematically represented by the equation:

:math:`A\left(i,j\right)=\ \sum_{k=1}^{n}{a_k(i,j)}`

Matrix A denotes the integrated 3D genome data from multiple cells, and a_k  represents the single-cell 3D genome data from the k-th cell. The indices i and j indicate the genomic positions of interactions. This aggregation function serves to generate a pseudo-bulk dataset that represents cell types or clusters of interest. The creation of pseudo-bulk data allows for the enhancement of signal-to-noise ratio and can reveal patterns that may not be evident when analyzing single cells in isolation. This approach is particularly useful for downstream analyses such as identifying consensus chromosomal structures, studying cell type–specific interactions, and exploring the organization of the genome across different cellular states.



If you want to integrate the sc3DG data in the entire directory, you can use the following command:

.. code-block:: shell

    stark --mcool ./data/mcool \
        --output ./output/merge.mcool \
        --resolution 10000,400000,1000000 \
        --num 100 \
        --sort True \
        --n_jobs 32

else if you want to integrate the sc3DG data which you prefer, you can use the following command:

.. code-block:: shell

    stark --mcool ./data/mcool \
        --output ./output/merge.mcool \
        --resolution 10000,400000,1000000 \
        --num 100 \
        --sort True \
        --n_jobs 32 \
        --select ./data/select.txt

The select.txt file is a file that contains the names of the sc3DG data you want to integrate.
The file format is as follows:

.. code-block:: shell

    cell_1.mcool
    cell_2.mcool
    cell_3.mcool


