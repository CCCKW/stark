Spatial Structure Capture Efficiency
=========================================

Spatial  Structure Capture Efficiency (SSCE) is introduced as a quantitative metric designed to evaluate the efficacy of single-cell sequencing in capturing the chromotin structural information. It is crucial to recognize that while the total number of sequenced contacts can vary significantly among cells, a lower count does not necessarily indicate inferior sequencing quality. Such a result may simply reflect the capture of a specific subset of the chromotin architecture, which can still have substantial biological implications. Conversely, a high number of detected contacts does not guarantee the biological relevance of these interactions or confirm an accurate representation of the cell's chromotin spatial organization. With this in mind, we introduce a novel indicator, termed Spatial Structure Capture Efficiency (SSCE), which quantifies the proportion of intra-cellular contacts contributing to the formation of the cell's chromotin architecture. This metric is essential for a more nuanced understanding of the structural fidelity of single-cell sequencing data.The SSCE are calculated as follows:
:math:`SSCE=\frac{T+A}{I+C}`,
where  :math:`T`, :math:`A`, :math:`I`, and :math:`C` are derived from fitting the following regression:

:math:`C_i～B+I*inter_i+T*tad_i+A*ab_i`

In this equation, :math:`C_i` denotes the aggregate count of contacts associated with chromosome :math:`i`, :math:`inter_i` represents the number of interactions that chromosome :math:`i` engages in with other chromosomes in the genome. :math:`tad_i`  corresponds to the count of TAD boundaries identified on chromosome :math:`i` at a resolution of 40 kilobases, which is pivotal for elucidating chromatin domain organization. :math:`ab_i` signifies the count of A/B compartment switches for chromosome :math:`i` determined at a resolution of 100 kilobases, providing insights into the chromatin's compartmentalization. The constant :math:`B` in the equation serves two critical roles. Firstly, it accounts for contacts that cannot be attributed to any specific chromosomal structure feature. Secondly, it includes contacts that may result from technical or systematic errors during the sequencing process. To ensure the accuracy of the model, any chromosomes where all variables are nullified are systematically excluded from the regression analysis.
Considering the inherent collinearity between the counts of :math:`ab_i` and :math:`tad_i`, we employ Ridge regression as our fitting method. The coefficients derived from the Ridge regression model are then utilized to compute the SSCE. It is logical to infer that higher values of structural features :math:`ab_i` and :math:`tad_i`, along with lower values of non-structural features :math:`B` and :math:`I`, contribute to higher SSCE scores. An elevated SSCE score, therefore, implies a more thorough and precise depiction of the single-cell chomotin structure.

To calculate the SSCE, you can use the following command:

.. code-block:: shell

    stark ssce --mcool /absolute/path/to/data/test_ssce \
    --genome /absolute/path/to/mm10/mm10.fa \
    --output /absolute/path/to/result/ssce \
    --nproc 2

