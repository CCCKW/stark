What's the STARK?
===================


.. image:: _static/stark_logo.png
    :width: 300px

sc3DG-seq has emerged as an essential tool for understanding the variability in 3D chromatin structure among individual cells. The process encompasses cell sorting, chromosome cross-linking, digestion, ligation, and sequencing, culminating in a map of interactions across the genome. Despite the development of over ten distinct sc3DG-seq techniques, each with unique processing requirements, a unified computational framework has been lacking. This gap poses a significant challenge to the research community, hindering the effective utilization of these extensive and valuable data.

To address this challenge, we propose STARK, a unified framework that conducts preprocessing, quality control and downstream analysis for all current sc3DG-seq data types (Fig S1). STARK encompasses three modules: Preprocess, Cell QC (Quality Control) and Downstream Analysis.


STARK is a  software for processing various types of sc3DG-seq sequencing data,
which currently includes the sc3DG-seq technologies as follows:

* scHi-C
* scHi-C+
* Dip-C
* HiRES
* sn-m3C
* scSPRITE
* sciHi-C
* snHi-C
* snHi-C
* scNanoHi-C
* scMethyl
* scCARE-seq