
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Overview

infercna aims to provide functions for inferring CNA values from
scRNA-seq data and related queries.

  - `infercna()` to infer copy-number alterations from single-cell
    RNA-seq data
  - `refCorrect()` to convert relative CNA values to absolute values +
    computed in `infercna()` if reference cells are provided
  - `cnaPlot()` to plot a heatmap of CNA values
  - `cnaScatterPlot()` to visualise malignant and non-malignant cell
    subsets
  - `cnaCor()` a parameter to identify cells with high CNAs + computed
    in `cnaScatterPlot()`
  - `cnaSignal()` a second parameter to identify cells with high CNAs +
    computed in `cnaScatterPlot()`
  - `findMalignant()` to find malignant subsets of cells
  - `findClones()` to identify genetic subclones
  - `fitBimodal()` to fit a bimodal gaussian distribution + used in
    `findMalignant()` + used in `findClones()`
  - `filterGenes()` to filter genes by their genome features
  - `splitGenes()` to split genes by their genome features
  - `orderGenes()` to order genes by their genomic position
  - `useGenome()` to change the default genome configured with infercna
  - `addGenome()` to configure infercna with a new genome specified by
    the user

*See Reference tab for a full list and documentation pages.*

## Installation

To install `infercna`:

``` r
# install.packages("devtools")
devtools::install_github("jlaffy/infercna")
```

## References

The methodology behind infercna has been tried and tested in several
high-impact publications. It was actually in the earliest of these
papers (last listed) that the idea to infer CNAs from single-cell
RNA-sequencing data was first formulated.

  - [An Integrative Model of Cellular States, Plasticity, and Genetics
    for Glioblastoma (Neftel, Laffy et al., 2019,
    *Cell*)](https://doi.org/10.1016/j.cell.2019.06.024)

  - [Developmental and oncogenic programs in H3K27M gliomas dissected by
    single-cell RNA-seq (Filbin, Tirosh et al., 2018,
    *Science*)](https://doi.org/10.1126/science.aao4750)

  - [Decoupling genetics, lineages, and microenvironment in IDH-mutant
    gliomas by single-cell RNA-seq (Venteicher, Tirosh et al., 2017,
    *Science*)](https://doi.org/10.1126/science.aai8478)

  - [Single-cell RNA-seq supports a developmental hierarchy in human
    oligodendroglioma (Tirosh, Venteicher et al., 2016,
    *Nature*)](https://doi.org/10.1038/nature20123)

  - [Single-cell RNA-seq highlights intratumoral heterogeneity in
    primary glioblastoma (Patel, Tirosh et al., 2014,
    *Science*)](https://doi.org/10.1126/science.1254257)

## Data requirements

The bare minimum for use in infercna is:

  - a single-cell expression matrix of genes by cells + *not* centered +
    normalised for sequencing depth and gene length (e.g. one of TPM,
    RPKM, CPM, etc). + *optionally* in log space. e.g. `log2(TPM/10
    + 1)` + Note: also see `infercna::TPM` and `infercna::logTPM`

If you would like to compute *absolute*, rather than *relative*, CNA
values, you should additionally provide:

  - a list of length two or more containing reference cell IDs of normal
    cells. For example list(macrophages, oligodendrocytes). + see
    example reference `infercna::refCells`

Finally, if your genome is not available in the current implementation
of infercna, you should additionally provide:

  - a genome dataframe, containing the columns: `symbol`,
    `chromosome_name`, `start_position`, `arm`.

## Example data

infercna is built with two example datasets of scRNA-seq data from two
patients with Glioblastoma, `infercna::bt771` and `infercna::mgh125`,
along with two normal reference groups, `infercna::refCells`. The
matrices are stored as sparse matrices and you can use
`infercna::useData()` to load them as normal matrices. These patients
are taken from a much larger cohort of 28 Glioblastoma samples. You can
look at the complete study
[here](Neftel*,%20Laffy*%20et%20al.%202019,%20Cell) and can download the
complete dataset via the [Single Cell
Portal](https://portals.broadinstitute.org/single_cell/study/SCP393/single-cell-rna-seq-of-adult-and-pediatric-glioblastoma).

## Future implementations

Future implementations will include:

  - more default genomes to choose from
  - option to correct CNA values (to absolute) when just *one* reference
    is available.
  - more stuff…
