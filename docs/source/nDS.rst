Calculate the normlized detection score
===========================================


Now, STARK support to calculate 4 kind of detection score, as follows:

* chromosome territory detection score
* chromosome TADs detection score
* chromosome A/B compartment detection score

It is worth noting that nDS calculations often take a long time.

The data needed for the test can be downloaded from `here <www.baidu.com>`_


introduction of normalized detection score
++++++++++++++++++++++++++++++++++++++++++++
Normalized Detection Scores (nDS) are utilized to evaluate the clarity of 3D genome structures, including chromosome territories, A/B compartments, and TADs. We have refined the approach introduced in scSPRITE, leveraging the high-throughput capabilities of single-cell sequencing while bypassing the requirement for pre-computed genome-wide 3D structural data. The Detection Scores (DS) are computed using the native contact matrix, which defines the frequency of interactions among various genomic bins. A higher DS for chromatin compartments in a cell, for instance, suggests a more pronounced preference for intra-chromosomal interactions over inter-chromosomal interactions.

Detection scores were calculated for each structure in each cell as follows:

Chromosome Territories. The unnormalized Detection Score (DS) for chromatin territories is calculated using the following formula:
:math:`DS_{terr} = \frac{O_{chr}^{intra}}{E_{chr}^{intra}} - \frac{O_{chr}^{inter}}{E_{chr}^{inter}}`

Here, :math:`O_{chr}^{intra}` and :math:`E_{chr}^{intra}` represent the observed and expected intra-chromosomal contact frequencies, respectively, while :math:`O_{chr}^{inter}` and :math:`E_{chr}^{inter}` denote the observed and expected inter-chromosomal contact frequencies. The DS of chromatin territory is calculated in different chromatin interactions, such as chr1-chr2, and does not include the same chromosome interaction, such as chr1-chr1, with a contact map resolution of 100 Kb. Thus, taking detection score of the chr1 and chr2 as an example, :math:`E_{chr}^{intra}` is calculated as :math:`N_1^c \times N_1^c + N_2^c \times N_2^c`, where :math:`N_1^c` and :math:`N_2^c` are the number of bins of chr1 and chr2. :math:`E_{chr}^{inter}` is calculated as :math:`N_1^c \times N_2^c \times 2`. Finally, :math:`DS_{terr}` of a single cell is the mean of the different chromatin interactions’ DS.

A/B compartments. In contrast to scSPRITE, which requires a prior bulk A/B compartment structural data, we leverage the high-throughput capacity of sc3DG-seq data. By merging the top 300 cells with the most contacts from each dataset or the entire dataset, we generate a pseudobulk data. This data is used to identify A/B compartments, creating a prior compartment annotation for subsequent single-cell DS calculation. This approach avoids the need for prior bulk data, which may not always be available. The unnormalized DS of A/B compartments is calculated as follows:
:math:`DS_{AB} = \frac{O_{AB}^{intra}}{E_{AB}^{intra}} - \frac{O_{AB}^{inter}}{E_{AB}^{inter}}`

Here, :math:`O_{AB}^{intra}` and :math:`E_{AB}^{intra}` represent the observed and expected contact frequencies within compartments, respectively, while :math:`O_{AB}^{inter}` and :math:`E_{AB}^{inter}` denote the observed and expected contact frequencies between compartments. The :math:`DS_{AB}` of a cell is obtained by calculating the A/B compartment switches, which are transitions like ‘A-B-A’ or ‘B-A-B’, derived from the pseudobulk. Thus, taking a ‘A_1-B_1-A_2’ as an example, :math:`E_{AB}^{intra}` is calculated as :math:`N_1^A \times N_1^A + N_1^B \times N_1^B + N_2^A \times N_2^A`, where :math:`N_1^A`, :math:`N_1^B` and :math:`N_2^A` are the number of bins of A_1 compartment, B_1 compartment and A_2 compartment. :math:`E_{chr}^{inter}` is calculated as :math:`N_1^B \times (N_1^A + N_2^A) \times 2`. Finally, :math:`DS_{AB}` of a single cell is the mean of all transitions from pseudobulk. The calculation is performed on contact maps with a resolution of 100 Kb.



calculate the chromosome territory detection score
++++++++++++++++++++++++++++++++++++++++++++++++++++
.. code-block:: shell

    stark nDS --mcool /absolute/path//data/nDS \
        --output /absolute/path//result/nDS/terr.txt\
        --resolution 1000000 \
        --mode terr \
        --top_num 300 \
        --epoch 50 \
        --genome /absolute/path//hg38/hg38.fa

If top_num is -1, the nDS of all cells is calculated; otherwise, topN is calculated.
And The recommended setting for epoch is 10 or more

And the output like this:

+------------+------------+------------+
|            | chr1_chr2  | chr1_chr3  |
+============+============+============+
| Cell 1     | 0.0549252  | 0.0496     |
+------------+------------+------------+
| Cell 2     | 0.0549252  | 0.05123    |
+------------+------------+------------+

Row is the cell, and column is the chromosome interaction,
You can average the values in the same row to get the nDS of the cell



calculate the chromosome TADs detection score
+++++++++++++++++++++++++++++++++++++++++++++++
.. code-block:: shell

    stark nDS --mcool ./data/mcool \
        --output ./output/tads.nds.txt \
        --resolution 40000 \
        --mode tads \
        --top_num -1 \
        --describe test \
        --epoch 10 \

It is worth noting that nDS of TADs calculation takes a long time


calculate the chromosome A/B compartment detection score
++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.. code-block:: shell

    stark nDS --mcool ./data/mcool \
        --output ./output/compartment.nds.txt \
        --resolution 1000000 \
        --mode compartment \
        --top_num -1 \
        --describe test \
        --epoch 10 \
        --genome ./data/genome.fa

If the genome file is not provided, the nDS of A/B compartment will not be calculated


And the output like this:

+------------+------------+
|            | ab         |
+============+============+
| Cell 1     | .0197576   |
+------------+------------+
| Cell 2     | 0.01627818 |
+------------+------------+


