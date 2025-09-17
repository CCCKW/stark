EmptyCells
==========

sc3DG-seq technologies
----------------------

sc3DG-seq technologies are broadly categorized into two approaches:

1. Low-throughput non-barcoded methods
2. High-throughput barcoded methods

Non-barcoded methods:

- Involve the precise isolation of individual cells
- Subjected to library construction and sequencing

Barcoded methods:

- Allow for the high-throughput sequencing of multiple cells simultaneously by incorporating cell barcodes
- Enable the sequencing of a larger number of cells per experiment
- Can result in a substantial variation in the number of contacts captured from each cell

Importance of Quality Control
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The effective filtering of low-quality cells is critical.

Common quality control approaches:

- Rely on total contact number thresholds to filter cells/barcodes
- Have limitations:
   - Library preparation efficiencies can vary among barcodes
   - Stringent thresholds may inadvertently exclude cells with less compact chromatin structures
   - Could be particularly problematic when analyzing cell-state transitions (Nagano et al., 2017)

EmptyCells Method
-----------------

Inspired by EmptyDrops (Lun et al., 2019) and Cellranger (10 x Genomics, version 7.0.1), we propose a two-step method for filtering low-quality cells in sc3DG-seq data:

1. Preliminary filtering:
   - Setting a low threshold based on the total number of contacts
   - Remove low-quality cells with insufficient sequencing or capture efficiency

2. Further refinement:
   - Modeling the profile of these low-quality cells
   - Enhance the filtering process

Hypothesis
^^^^^^^^^^

High-quality cells exhibit a diverse array of interaction types.

Mathematical Representation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- All barcodes of a sample can be represented as :math:`A= \{y_{bg}\in \mathbb{R}^{N\times\mathcal{g}}\}`
- :math:`y_{bg}` denotes the number of contacts for interaction type g in barcode b
- N denotes the total number of barcodes
- :math:`\mathcal{g}` denotes the interaction type

Note: We exclude interaction types that are not present in all barcodes from our analysis to maintain consistency and accuracy.

Step 1: Threshold Establishment
-------------------------------

We establish a threshold T based on the total number of contacts per barcode:

.. math::

   T = \max(500, \alpha f(A_b))

Where:

- :math:`\alpha` is a hyperparameter, which we set to 0.01
- :math:`f(\cdot)` calculates the median value of the barcode contact vector :math:`A_b`
- :math:`A_b = \sum_{g\in\mathcal{g}} y_{gb}`

Barcodes that fall below this threshold T are discarded, while the remaining barcodes proceed to the next step of the filtering process.

Step 2: Low-Quality Barcode Profile Construction
------------------------------------------------

1. Select a subset of barcodes characterized by a low number of contacts (bottom 5%) to form a set G
2. Calculate A as the sum of contacts across all barcodes in G:

   .. math::

      A_g = \sum_{b\in G}{y_{gb},}

3. Yielding a contact vector :math:`A=(A_1,\ldots,A_\mathcal{g})`
4. Apply the Simple Good-Turing algorithm to obtain the posterior expectation :math:`\hat{p}_g`

Barcode Quality Assessment
--------------------------

For a given barcode b:

- The total number of contacts is denoted as :math:`t_b`
- Modeled using a Dirichlet multinomial distribution
- The likelihood is defined as:

  .. math::

     L_b = \frac{t_b!\Gamma(\alpha)}{\Gamma(t_b+\alpha)}\prod_{g=1}^{N}\frac{\Gamma(y_{bg}+\alpha_g)}{y_{bg}!\Gamma(\alpha_g)}

  where :math:`\alpha_g = \beta\hat{p}_g`, :math:`\beta` is a scaling factor

P-value Calculation
-------------------

We employ the Monte Carlo simulation method to calculate the p-value for each barcode:

1. For each iteration i of a given barcode, a new contact number vector is randomly sampled
2. The likelihood :math:`{L^\prime}_{bi}` is calculated based on :math:`\hat{p}_g`, :math:`t_b` and :math:`\beta`
3. The count :math:`R_b` of likelihoods from the R round simulation that are less than :math:`L_b` is utilized
4. P-value calculation (Phipson and Smyth, 2010):

   .. math::

      P_b = \frac{R_b+1}{R+1}

This approach allows us to quantitatively assess the quality of each barcode and to filter out those that do not meet the predetermined quality thresholds, thereby enhancing the reliability of downstream analyses.



To conduct the EmptyCells method, use the command below:

.. code-block:: bash

   stark emptycells --mcool-path /absolute/path/to/mcool \
        --output /absolute/path/to/output


The results includes two csv file and one rank-plot png file:
    1. res.csv: in this file, the column "ambient" is False means the barcode is a high-quality barcode, and True means the barcode is a low-quality barcode.
    2. matrix.csv: it records the data for calculating the p-value.
    3. rank-plot.png: the rank-plot of the cells.