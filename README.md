# Genomic Prediction Pipeline

A simulation-based genomic prediction framework comparing statistical methods (GBLUP, BayesA, BayesR) with machine learning approaches (Random Forest, XGBoost, SVR).

## Overview

This pipeline consists of two main components:
1. **Simulation** (`simulation.R`): Generates synthetic genomic data using AlphaSimR
2. **Prediction** (`prediction_pipeline.R`): Trains and evaluates genomic prediction models

## Requirements

### R Packages

```r
# Core packages
install.packages(c("AlphaSimR", "yaml", "jsonlite", "tidyverse"))

# Statistical genomic prediction
install.packages(c("rrBLUP", "hibayes"))

# Machine learning
install.packages(c("randomForest", "xgboost", "e1071"))

# Visualization
install.packages("gridExtra")
```

## Data and pipeline download

All genotype and phenotype data can be downloaded via [huggingface.co](https://huggingface.co/datasets/bowo1745/Genomic_prediction_simulation_data/tree/main) using `huggingface_hub` using the following command, but ensure you have installed conda or python3

```bash
conda install huggingface_hub # via conda or

#pip install huggingface_hub # via pip

mkdir -p output/data_gen

huggingface-cli download bowo1745/Genomic_prediction_simulation_data \
  --repo-type dataset \
  --local-dir output/data_gen \
  --include "run_scenario*/*"
```

## Configuration

Edit `config.yaml` to customize simulation parameters:

```yaml
# Key parameters
random_seed: 12345

population:
  n_individuals: 2000          # Population size
  n_chromosomes: 25            # Number of chromosomes
  n_snps_per_chrom: 1000       # SNPs per chromosome
  effective_pop_size: 300      # Ne - controls genetic diversity
  heritability: 0.2            # Target h²

breeding:
  burn_in_generations: 10      # Equilibrium phase (discarded)
  reference_generations: [14, 15]  # Training data
  validation_generation: 16    # Hyperparameter tuning
  test_generations: [16, 17, 18]   # Model evaluation
  selection_intensity: 0.8     # Proportion selected

genetic_architecture:
  n_qtl: 500                   # Number of QTLs
  heritability: 0.2            # Narrow-sense h²
  qtl_effect_distribution: "normal"  # or "mixture"

output:
  output_dir: "sim_output/run_scenario1"
```

### Genetic Architecture Options

| Parameter | Polygenic | Oligogenic |
|-----------|-----------|------------|
| n_qtl | 500-1000 | 50-100 |
| heritability | 0.2-0.4 | 0.4-0.6 |
| qtl_effect_distribution | "normal" | "mixture" |

## Running the Pipeline

### Step 1: Run Simulation

```bash
Rscript simulation.R
```

**Outputs** (saved to `output_dir`):
- `genotype_reference.csv` - Training genotypes
- `phenotype_reference.csv` - Training phenotypes
- `genotype_test_gen*.csv` - Test genotypes
- `phenotype_test_gen*.csv` - Test phenotypes
- `simulation_parameters.json` - Reproducibility info

### Step 2: Run Prediction

```bash
Rscript prediction_pipeline.R
```

**Outputs** (saved to `prediction_output/`):
- `model_performance_metrics.csv` - All evaluation metrics
- `predictions_gen*.csv` - Individual predictions
- `*.png` - Visualization plots

## Models Implemented

### Statistical Methods
| Model | Package | Description |
|-------|---------|-------------|
| GBLUP | rrBLUP | Genomic BLUP with GRM |
| BayesA | hibayes | t-distributed priors |
| BayesR | hibayes | Mixture model with spike |

### Machine Learning
| Model | Package | Key Parameters |
|-------|---------|----------------|
| Random Forest | randomForest | ntree=500, mtry=200 |
| XGBoost | xgboost | eta=0.1, max_depth=4 |
| SVR | e1071 | kernel=radial, cost=1 |

## Evaluation Metrics

- **Correlation**: Pearson correlation with true breeding values
- **R²**: Coefficient of determination
- **RMSE**: Root mean squared error
- **Bias slope**: Regression slope (ideal = 1.0)
- **Rank correlation**: Spearman correlation for selection

## Project Structure

```
├── config.yaml              # Simulation parameters
├── simulation.R             # AlphaSimR simulation
├── prediction_pipeline.R    # Model training & evaluation
├── sim_output/              # Simulation outputs
│   └── run_scenario1/
└── prediction_output/       # Prediction results
    └── run_scenario1/
```

## Citation

If using this pipeline, please cite:
- AlphaSimR: Gaynor et al. (2021)
- hibayes: Yin et al. (2023)
- rrBLUP: Endelman (2011)
