# RIMA: RIgorous Matching of single-cell transcriptomics Atlases

## Overview

RIMA (RIgorous Matching of single-cell transcriptomics Atlases) is an R package for matching single-cell transcriptomics data across two datasets at the level of **cell neighbourhoods**. It is particularly useful in scenarios where data integration is difficult, such as cross-species analyses, or comparisons across different experimental conditions.

Unlike most integration-based approaches, RIMA identifies an explicit mapping between cell neighbourhoods (small groups of similar cells) and retains quantitative gene expression information. This mapping enables a comparison between atlases on a "like for like" basis and allows to:

- **Compare cell states** across different datasets in a rigorous, statistically sound manner, with minimal preprocessing and assumptions
- **Conserve the resolution offered by single-cell data** by avoiding cell aggregations (e.g. clustering, pseudo-bulking) or subsetting to known gene sets
- **Retain actual gene expression values** for interpretable downstream analysis
- **Analyse associations between neighbourhood's metadata**, for example to align cell trajectories across conditions

## Key Features

- **Neighbourhood-based matching**: Works at the level of cell neighbourhoods rather than individual cells or clusters, avoiding averaging out biological signal while improving robustness against sequencing noise
- **Statistical significance testing**: Derives p-values for neighbourhood pairs similarities by creating a baseline with scrambled cell identity
- **Gene conservation scoring**: Identifies genes with conserved expression (CoPE score) across matched neighbourhoods
- **Visualise matching and metadata association**: Includes functions for creating informative heatmaps and embedding visualizations of RIMA's matching
- **Flexible gene mapping**: Supports feature mapping between datasets (e.g., for cross-species comparisons using orthologs)

## Installation

You can install RIMA directly from GitHub using the `devtools` package:

```r
# Install devtools if you don't have it
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

# Install RIMA from GitHub
devtools::install_github("majpark21/RIMA")
```

### Requirements

RIMA requires R ≥ 4.0 and the following packages:
- `data.table`
- `miloR` (Bioconductor)
- `SingleCellExperiment` (Bioconductor)
- `Matrix`
- `igraph`
- `ggplot2`
- `viridis`
- `rdist`

## Quick Start

The best way to get started with RIMA is to go through the usage example in the [vignettes directory](vignettes/).  It can also be accessed by running `browseVignettes("RIMA")`.

Here's a minimal example to get you started:

```r
library(RIMA)
library(miloR)
library(SingleCellExperiment)

# Assuming you have two SingleCellExperiment objects, containing sc data
# Here we load the built-in example datasets of mouse and rabbit gastrulation
sce_mouse <- RIMA::sce_mouse_gastrulation
sce_rabbit <- RIMA::sce_rabbit_gastrulation

# Step 0: Define the neighbourhoods (here with Milo's implementation, but could use others, e.g. metacells)
define_neighbourhoods <- function(sce, prop_seeds, knn=10, reduced.dim="PCA"){
  n_components <- ncol(reducedDim(sce, reduced.dim))  # use all available PCs
  mi <- Milo(sce)
  mi <- miloR::buildGraph(mi, k = knn, d = n_components, reduced.dim = reduced.dim)
  mi <- miloR::makeNhoods(mi, prop = prop_seeds, k = knn, d=n_components, refined = TRUE, reduced_dims = reduced.dim)
  return(mi)
}
mi_mouse <- define_neighbourhoods(sce_mouse, prop_seeds = 0.02)
mi_rabbit <- define_neighbourhoods(sce_rabbit, prop_seeds = 0.02)

# Step 1: Preprocess the Milo objects
milos <- preprocess_milos(mi_mouse, mi_rabbit)

# Step 2: Calculate neighbourhood similarities
dt_sims <- calculate_similarities(milos, method = "spearman")

# Step 3: Assess statistical significance of nhood-nhood similarity
dt_sims_sig <- calculate_nhoodnhood_significance(
  milos, dt_sims,
  n_scrambles = 10,
  col_scramble_label = "celltype",
  direction = "b"
)

# Step 4: Match significant nhood-nhood connections
dt_match <- match_nhoods(dt_sims_sig[is_significant == TRUE])

# Step 5: Visualize and analyze results
plot_matches_embed(milos, dt_match, cols_color = c("celltype", "celltype"), dimred="PCA")
plot_matches_map(milos, dt_match, cols_label = c("celltype", "stage"))

# Example downstream analysis: Find the 3 genes with the most conserved expression across matches
dt_cope <- calculate_cope(milos, dt_match, genes = NULL)
dt_cope <- dt_cope[order(dt_cope$cope, na.last = FALSE), ]
plot_paired_expression(milos, dt_match, genes = tail(dt_cope$gene, 3))
```

## Citation

If you use RIMA in your research, please cite:

```
Jacques, M.-A., et al. (2025). RIMA: Rigorous Matching of single-cell transcriptomics Atlases. 
GitHub: https://github.com/ma-jacques/RIMA
```

## Example Data

RIMA includes example data from mouse and rabbit gastrulation atlases:
- `sce_mouse_gastrulation`: Mouse gastrulation SingleCellExperiment
- `sce_rabbit_gastrulation`: Rabbit gastrulation SingleCellExperiment

Load them with:
```r
RIMA::sce_mouse_gastrulation
RIMA::sce_rabbit_gastrulation
```

Both example datasets are small subsets of large gastrulation atlases, see:
- Pijuan-Sala, Blanca, et al. "A single-cell molecular map of mouse gastrulation and early organogenesis." Nature 566.7745 (2019): 490-495.
- Imaz-Rosshandler, Ivan, et al. "Tracking early mammalian organogenesis–prediction and validation of differentiation trajectories at whole organism scale." Development 151.3 (2024): dev201867.
- Ton, Mai-Linh Nu, et al. "An atlas of rabbit development as a model for single-cell comparative genomics." Nature cell biology 25.7 (2023): 1061-1072.

## License

RIMA is licensed under the GPL-3 License. See [LICENSE](LICENSE) for details.

## Contact

For questions or feedback, please contact:
- **Marc-Antoine Jacques** - jacques@ebi.ac.uk
