Get imputed contact map
===========================

The RWR algorithm is a variant of the traditional random walk technique. It is designed to explore the graph by randomly walking through its nodes, with the added feature that at each step, there is a probability the walk will "restart" and return to a predefined starting node. We implement this algorithm, as referenced in previous research, for the interpolation of single-cell 3D genome data.

First, we construct a probability transition matrix by normalizing the single-cell contact matrix (A) such that each row sums to 1. The walk restarts at the starting node with a probability :math:`𝑟`, and with the probability :math:`1−𝑟`, it proceeds to a neighboring node. The goal is to find the steady-state π which represents the matrix given the restart probability :math:`𝑟`. The RWR is solved iteratively through the following steps:

1.	Initialization: Start with an initial probability vector :math:`𝜋0 = A`.

2.	Iteration: Update 𝜋 using the formula: :math:`𝜋𝑡+1 = (1−𝑟)𝑃 𝜋𝑡 + 𝑟A`. :math:`P` is XXX.

3.	Convergence: Repeat the iteration until π converges, when :math:`|| 𝜋𝑡+1 - 𝜋𝑡 ||< 10^{-6}`.

Upon completion of the iterations, we obtain 𝜋𝑡+1 which serves as the imputated single-cell contact matrix.




.. code-block:: shell

    stark impute --mcool-path/absolute/path/to/data/gini/dipC_GSE117874_hum \
    --resolution 1000000 \
    --output/absolute/path/to/result/impute \
    --nthread 2



The result like this:

.. image:: ./_static/impute.png
    :width: 400
    :align: center

