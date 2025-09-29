# MethXY

**1. Specialized for sex chromosome analysis, with full autosome support**
**MethXY** is a flexible and comprehensive R package designed for **DNA methylation analysis**, with a particular focus on **sex chromosome methylation**.
It is especially suitable for studies exploring methylation differences on the X and Y chromosomes while also supporting the simultaneous analysis of **sex chromosomes and autosomes.**

**2.Flexible input handling with unified QC and preprocessing**
Public DNA methylation datasets are typically provided in **two common format**s:
**IDAT** files generated directly by Illumina methylation arrays.
Processed **intensity data**, often shared in repositories such as GEO.

MethXY accepts both formats and performs standardized quality control (QC) and preprocessing, ensuring consistent downstream analysis.
This includes essential steps such as sex prediction, age estimation, and cell composition deconvolution, which are critical for adjusting biological and technical confounders.

**3. Comprehensive methylation analysis: DMP and VMP**
MethXY provides a full suite of methylation analysis functions, allowing researchers to explore both **differentially methylated positions (DMPs)** and **variably methylated positions (VMPs)**:
