# Introgression_mapping

This pipeline identifies and visualises candidate introgressions in wheat and wild relative genomes by analysing genome-wide coverage deviation. It processes raw sequencing data, maps reads to a combined reference, computes normalised coverage profiles, and detects introgressed regions based on characteristic signal shifts.

This pipeline is combined approach based on the [alien_detection](https://github.com/benedictcoombes/alien_detection) pipeline by Coombes et al. (2023)and the [SkimSeq_Method](https://github.com/sandeshsth/SkimSeq_Method) method by Adhikari et al. (2022).

These resources laid the groundwork for coverage-based introgression detection in low-coverage resequencing data.

The core structure and logic are inherited from these pipelines, but this implementation uses a custom coverage deviation script:
cov_deviation_new.py, which introduces:

Asymmetric normalisation for wheat and wild relative chromosomes

Dynamic scaling for wild relative chromosomes based on the top 10% coverage bins

Output formatted for compatibility with high-resolution R plotting

This enhancement improves detection sensitivity for subtle introgressions in ultra-low coverage datasets.
