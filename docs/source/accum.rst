Accumlated Analysis
========================

For a comprehensive quality analysis of single-cell data,
we implemented an iterative strategy to cumulatively aggregate individual cells based on either
the top 300 cells with the highest contact counts or the entire dataset. T
his approach allows us to evaluate the consistency of 3D genome structure detection across cells
within a sc3DG-seq dataset.
To measure the similarity between the aggregated cells at each iteration and the final aggregate,
we calculate the Intersection over Union (IoU) of the TAD boundaries.
The IoU is calculated using the following formula:

.. math::
   IOU_i = \frac{A_i \cap T}{A_i \cup T}

Here, :math:`A_i` represents the set of TAD boundaries identified in the aggregated single-cell
data at the :math:`i^{th}` iteration, and T represents the set of TAD boundaries
in the final aggregated cell. The IoU is a metric that quantifies the overlap
between two sets of TAD boundaries. A higher IoU value signifies a greater
degree of overlap, indicating higher consistency in the detected 3D structures.

.. code-block:: shell

    stark accum   --mcool /cluster/home/Gaoruixiang/test/example/data/nDS \
    --output /cluster/home/Gaoruixiang/test/example/result/accum \
    --resolution 40000 \
    --top_num 30 \
