#!/usr/bin/env Rscript

#===============================================================================
# AlphaSimR Simulation Script
# 
# This script simulates a complete genetic population with realistic structure
# using AlphaSimR, including genetic drift, recombination, and selection.
#
# Output:
# - genotype_data.csv: Genotype matrix (n_individuals × n_snps)
# - phenotype_data.csv: Phenotypes and True Breeding Values
# - simulation_parameters.json: Simulation parameters for reproducibility
#
# Author: Agus Wibowo
# Date: 2025-01-26
#===============================================================================

#install.packages(c("AlphaSimR", "yaml", "jsonlite", 
#                   "tidyverse"))

# Load required libraries
suppressPackageStartupMessages({
  library(AlphaSimR)
  library(yaml)
  library(jsonlite)
  library(tidyverse)
})

#' Load configuration from YAML
#' 
#' @param config_file Path to configuration file
#' @return List containing configuration parameters
load_config <- function(config_file = "config.yaml") {
  if (!file.exists(config_file)) {
    stop(sprintf("Configuration file not found: %s", config_file))
  }
  
  config <- yaml::read_yaml(config_file)
  cat(sprintf("✓ Configuration loaded from: %s\n", config_file))
  
  return(config)
}

#' 
#' @param config 
print_config_summary <- function(config) {
  cat("\n=== Configuration Summary ===\n")
  cat(sprintf("Population: %d individuals, %d chromosomes, %d SNPs per chromosome\n",
              config$population$n_individuals,
              config$population$n_chromosomes,
              config$population$n_snps_per_chrom))
  cat(sprintf("Breeding: %d generations, selection intensity %.2f\n",
              config$breeding$generations,
              config$breeding$selection_intensity))
  cat(sprintf("Genetic Architecture: %d QTL, h2 = %.2f\n",
              config$genetic_architecture$n_qtl,
              config$genetic_architecture$heritability))
  cat("\n")
}

#' @param config 
#' @return AlphaSimR founder population object
initialize_founder_population <- function(config) {
  
  cat("=== STEP 1: Initialize Founder Population ===\n")
  
  # Parameter dari config
  seed <- as.integer(config$random_seed)
  n_ind <- as.integer(config$population$n_individuals)
  n_chr <- as.integer(config$population$n_chromosomes)
  n_snps <- as.integer(config$population$n_snps_per_chrom)
  chr_length_bp <- as.numeric(config$population$chr_length_bp)
  gen_len <- as.numeric(config$population$genome_length)
  ne <- as.numeric(config$population$effective_pop_size)
  hist_gen <- as.numeric(config$population$historical_generations)
  hist_ne <- as.numeric(config$population$historical_Ne)
  mutRate <- as.numeric(config$population$mutation_rate)

  # calculate SNP density
  snp_density <- (n_snps * n_chr) / ((chr_length_bp * n_chr) / 1000000)  # Total SNPs per total Mb
  
  cat(sprintf("Creating founder population with:\n"))
  cat(sprintf("  - %d chromosomes\n", n_chr))
  cat(sprintf("  - %d SNPs per chromosome (total: %d SNPs)\n", n_snps, n_chr * n_snps))
  cat(sprintf("  - Effective population size (Ne): %d\n", ne))
  cat(sprintf("  - Historical generations: %d\n", hist_gen))
  cat(sprintf("SNP density: %.2f SNPs per Mb\n", snp_density))
  
  set.seed(seed)
  # Create founder population using coalescent simulation
  founderPop <- runMacs2(
    nInd = n_ind,
    nChr = n_chr,
    segSites = n_snps,
    bp = chr_length_bp,
    genLen = gen_len,
    Ne = ne,
    mutRate = mutRate,
    histNe = hist_ne,
    histGen = hist_gen
  )
  
  # Validation
  actual_snps <- sum(sapply(founderPop@genMap, length))
  cat(sprintf("✓ Founder population created with %d total SNPs\n", actual_snps))
  cat(sprintf("  Average SNPs per chromosome: %.0f\n", actual_snps / n_chr))
  
  return(founderPop)
}


#' 
#' @param founderPop AlphaSimR founder population
#' @param config
#' @return Final population after breeding
run_breeding_program <- function(founderPop, config) {
  
  cat("=== STEP 2: Run Breeding Program ===\n")
  
  # Initialize simulation parameters
  SP <- SimParam$new(founderPop)
  assign("SP", SP, envir = .GlobalEnv)
  SP$setSexes("no")  # Hermaphrodite 
  SP$setTrackPed(TRUE)
  
  # Add traits with QTL
  n_qtl <- as.integer(config$genetic_architecture$n_qtl)
  h2 <- as.numeric(config$genetic_architecture$heritability)
  
  cat(sprintf("Adding trait with:\n"))
  cat(sprintf("  - %d QTL\n", n_qtl))
  cat(sprintf("  - Target heritability: %.2f\n", h2))

  # Setup QTL effect distribution
  n_qtl_per_chr <- as.integer(n_qtl %/% config$population$n_chromosomes)
  total_qtl <- n_qtl_per_chr * as.integer(config$population$n_chromosomes)

  if (config$genetic_architecture$qtl_effect_distribution == "mixture") {
    SP$addTraitAG(
      nQtlPerChr = n_qtl_per_chr,
      mean = 0,
      var = 1,
      gamma = TRUE,
      shape = 0.4,  # Controls mixture: many small + few large effects
      name = "trait1"
    )
    cat(sprintf("✓ QTL added: %d (gamma distribution for realistic mixture)\n", total_qtl))
    
  } else {
    # Default: normal distribution
    SP$addTraitA(
      nQtlPerChr = n_qtl_per_chr,
      mean = 0,
      var = 1,
      name = "trait1"
    )
    cat(sprintf("✓ QTL added: %d (normal distribution)\n", total_qtl))
  }
  
  # Rescale to target heritability
  SP$setVarE(h2 = h2)

  # SNP access strategy: Use all segregating sites directly
  # Use ALL segregating sites as SNP markers
  cat(sprintf("\nDefining SNP chip markers:\n"))
  SP$addSnpChip(nSnpPerChr = SP$segSites)  # Use all segregating sites
  total_seg_sites <- sum(SP$segSites)
  cat(sprintf("✓ SNP chip defined: %d total SNPs\n", total_seg_sites))
  cat(sprintf("  Average per chromosome: %.1f SNPs\n", mean(SP$segSites)))
  cat(sprintf("  Range: %d - %d SNPs per chromosome\n", 
              min(SP$segSites), max(SP$segSites)))
  cat(sprintf("  This ensures QTLs are within genotyped markers\n\n"))
  
  # Create initial population
  pop <- newPop(founderPop, simParam = SP)
  pop <- setPheno(pop, simParam = SP)

  # VALIDATE realized heritability
  varG <- varG(pop)
  varP <- varP(pop)
  realized_h2 <- varG / varP

  cat(sprintf("✓ Realized heritability: %.3f (target: %.3f)\n", realized_h2, h2))
  varG_initial <- varG(pop)
  cat(sprintf("  Initial genetic variance (VarG): %.3f\n", varG_initial))
  cat(sprintf("  Initial phenotypic variance (VarP): %.3f\n", varP(pop)))
  cat(sprintf("  Initial heritability: %.3f\n", varG_initial / varP(pop)))

  if (abs(realized_h2 - h2) > 0.05) {
    warning(sprintf(
      "Realized h² (%.3f) deviates from target (%.3f) by >0.05",
      realized_h2, h2
    ))
  }
  
  # Store for multi-generational populations
  pop_history <- list()
  
  # Extract configuration
  burn_in <- as.integer(config$breeding$burn_in_generations)
  ref_gens <- as.integer(config$breeding$reference_generations)
  val_gen <- as.integer(config$breeding$validation_generation)
  test_gens <- as.integer(config$breeding$test_generations)
  total_gens <- max(c(ref_gens, val_gen, test_gens))
  
  sel_intensity <- as.numeric(config$breeding$selection_intensity)
  
  cat(sprintf("\nRunning breeding program:\n"))
  cat(sprintf("  - Burn-in: %d generations (discarded)\n", burn_in))
  cat(sprintf("  - Reference: generations %s (training)\n", 
              paste(ref_gens, collapse=", ")))
  cat(sprintf("  - Validation: generation %d (tuning)\n", val_gen))
  cat(sprintf("  - Test: generations %s (accuracy evaluation)\n", 
              paste(test_gens, collapse=", ")))
  cat(sprintf("  - Selection intensity: %.2f\n", sel_intensity))
  cat(sprintf("  - Mating design: %s\n", config$breeding$mating_design))
  
  # PHASE 1: Burn-in generations (random mating, no selection, no storage)
  if (burn_in > 0) {
    cat(sprintf("\n=== Phase 1: Burn-in (%d generations) ===\n", burn_in))
    cat("Purpose: Establish mutation-drift equilibrium and LD structure\n")
    
    for (gen in 1:burn_in) {
      if (gen %% 5 == 0 || gen == burn_in) {
        cat(sprintf("  Burn-in generation %d/%d...\n", gen, burn_in))
      }
      
      # Random mating only (no selection pressure)
      pop <- randCross(pop, nCrosses = config$breeding$population_size, simParam = SP)
      pop <- setPheno(pop, simParam = SP)

      if (gen %% 5 == 0 || gen == burn_in) {
        varG_gen <- varG(pop)
        varA_gen <- varA(pop)
        cat(sprintf("    VarG=%.3f, VarA=%.3f\n", varG_gen, varA_gen))
      }
    }
    
    cat("✓ Burn-in complete. Population at mutation-drift equilibrium.\n\n")
  }
  
  # PHASE 2: Breeding program with selection
  cat("=== Phase 2: Breeding Program with Selection ===\n")

  for (gen in 1:total_gens) {
    if (gen %% 5 == 0 || gen == total_gens) {
      cat(sprintf("  Generation %d/%d (absolute gen: %d)...\n", gen, total_gens, gen + burn_in))
    }

    # Apply truncation selection for all generations
    n_elite <- round(nInd(pop) * sel_intensity)
    elite_pop <- selectInd(pop, nInd = n_elite, use = "pheno", simParam = SP)
    
    # Mating from elite parents
    if (config$breeding$mating_design == "random") {
      pop <- randCross(elite_pop, 
                       nCrosses = config$breeding$population_size, 
                       simParam = SP)
    } else if (config$breeding$mating_design == "assortative") {
      pop <- selectCross(elite_pop, 
                         nCrosses = config$breeding$population_size,
                         use = "pheno", simParam = SP)
    }
    
    pop <- setPheno(pop, simParam = SP)

    # === VARIANCE TRACKING ===
    varG_gen <- varG(pop)
    varA_gen <- varA(pop)
    realized_h2 <- varG_gen / varP(pop)
    
    if (gen %% 5 == 0 || gen == total_gens) {
      cat(sprintf("    Gen %d: VarG=%.3f, VarA=%.3f, h²=%.3f\n",
                  gen, varG_gen, varA_gen, realized_h2))
    }
    
    # WARNING if variance drops too much
    if (gen > 1 && varG_gen < 0.5 * varG_initial) {
      warning(sprintf(
        "Generation %d: Genetic variance dropped >50%% (%.3f → %.3f). Selection too intense!",
        gen, varG_initial, varG_gen
      ))
    }
    
    # save population for reference or validation
    if (gen %in% ref_gens) {
      pop_history[[paste0("reference_gen", gen)]] <- pop
      cat(sprintf("    ✓ Saved reference generation %d\n", gen))
    }
    if (gen == val_gen) {
      pop_history[["validation"]] <- pop
      cat(sprintf("    ✓ Saved validation generation %d\n", gen))
    }
    if (gen %in% test_gens) {
      pop_history[[paste0("test_gen", gen)]] <- pop
      cat(sprintf("    ✓ Saved test generation %d\n", gen))
    }
  }
  
    cat("✓ Continuous breeding program completed\n\n")
  
    # === POPULATION DIVERGENCE DIAGNOSTICS ===
    cat("=== Population Divergence Analysis ===\n")
    
    # Get reference and test populations
    ref_pops <- pop_history[grep("reference_gen", names(pop_history))]
    test_pops <- pop_history[grep("test_gen", names(pop_history))]
    
    if (length(ref_pops) > 0 && length(test_pops) > 0) {
      # Use last reference generation
      ref_pop <- ref_pops[[length(ref_pops)]]
      
      for (test_name in names(test_pops)) {
        test_pop <- test_pops[[test_name]]
        test_gen <- as.integer(sub("test_gen", "", test_name))
        ref_last_gen <- max(ref_gens)
        
        cat(sprintf("\nReference gen %d vs Test gen %d:\n", ref_last_gen, test_gen))
        
        # 1. Allele frequency correlation
        ref_geno <- pullSegSiteGeno(ref_pop, simParam = SP)
        test_geno <- pullSegSiteGeno(test_pop, simParam = SP)
        
        ref_freq <- colMeans(ref_geno, na.rm = TRUE) / 2
        test_freq <- colMeans(test_geno, na.rm = TRUE) / 2
        
        freq_cor <- cor(ref_freq, test_freq, use = "complete.obs")
        cat(sprintf("  1. Allele frequency correlation: %.4f\n", freq_cor))
        
        # 2. Mean allele frequency difference
        mean_freq_diff <- mean(abs(ref_freq - test_freq), na.rm = TRUE)
        cat(sprintf("  2. Mean |Δ allele freq|: %.4f\n", mean_freq_diff))
        
        # 3. Number of fixed SNPs in test
        ref_fixed <- sum(ref_freq == 0 | ref_freq == 1, na.rm = TRUE)
        test_fixed <- sum(test_freq == 0 | test_freq == 1, na.rm = TRUE)
        cat(sprintf("  3. Fixed SNPs - Ref: %d, Test: %d (Δ=%d)\n", 
                    ref_fixed, test_fixed, test_fixed - ref_fixed))
        
        # 4. Genetic variance comparison
        ref_varG <- varG(ref_pop)
        test_varG <- varG(test_pop)
        varG_change <- (test_varG - ref_varG) / ref_varG * 100
        cat(sprintf("  4. VarG change: %.3f → %.3f (%.1f%%)\n", 
                    ref_varG, test_varG, varG_change))
        
        # 5. Mean breeding value shift
        ref_mean_bv <- mean(gv(ref_pop))
        test_mean_bv <- mean(gv(test_pop))
        bv_shift <- test_mean_bv - ref_mean_bv
        cat(sprintf("  5. Mean BV shift: %.3f → %.3f (Δ=%.3f)\n", 
                    ref_mean_bv, test_mean_bv, bv_shift))
        
        # === DIVERGENCE WARNING ===
        warning_issued <- FALSE
        critical_issue <- FALSE
        
        # CRITICAL warnings (major issues)
        if (freq_cor < 0.95) {
          cat(sprintf("  ⚠️  CRITICAL: Low allele freq correlation (%.3f < 0.95)\n", freq_cor))
          warning_issued <- TRUE
          critical_issue <- TRUE
        }
        
        if (mean_freq_diff > 0.10) {
          cat(sprintf("  ⚠️  CRITICAL: High mean freq difference (%.3f > 0.10)\n", mean_freq_diff))
          warning_issued <- TRUE
          critical_issue <- TRUE
        }
        
        if (abs(varG_change) > 25) {
          cat(sprintf("  ⚠️  CRITICAL: Large VarG change (%.1f%% > 25%%)\n", abs(varG_change)))
          warning_issued <- TRUE
          critical_issue <- TRUE
        }
        
        # MODERATE warnings (expected with selection, but notable)
        if (test_fixed - ref_fixed > 500) {
          cat(sprintf("  ⚠️  MODERATE: Many new fixed SNPs (%d > 500)\n", test_fixed - ref_fixed))
          warning_issued <- TRUE
        } else if (test_fixed - ref_fixed > 200) {
          cat(sprintf("  ℹ️  INFO: Moderate SNP fixation (%d SNPs) - expected with continuous selection\n", 
                      test_fixed - ref_fixed))
        }
        
        if (abs(varG_change) > 15 && abs(varG_change) <= 25) {
          cat(sprintf("  ⚠️  MODERATE: Notable VarG change (%.1f%% > 15%%)\n", abs(varG_change)))
          warning_issued <- TRUE
        }
        
        # Final assessment
        if (critical_issue) {
          cat("\n  ❌ CRITICAL: Populations show SEVERE divergence\n")
          cat("     → Genomic prediction will likely fail\n")
          cat("     → MUST reduce selection intensity or shorten gap\n")
        } else if (warning_issued) {
          cat("\n  ⚠️  MODERATE: Populations show moderate divergence\n")
          cat("     → Genomic prediction accuracy may be reduced\n")
          cat("     → Consider evaluating robustness across methods\n")
        } else {
          cat("\n  ✓ EXCELLENT: Populations are genetically similar\n")
          cat("     → Ideal for genomic prediction validation\n")
        }
      }
    }
    
    cat("\n")
    
    return(list(
      pop_history = pop_history,
      SP = SP,
      config_generations = list(
        reference = ref_gens,
        validation = val_gen,
        test = test_gens
      )
    ))
  }

#' Extract pedigree information from AlphaSimR populations
#' 
#' @param pop_history List of populations from breeding program
#' @param config_generations List with reference, validation, test generations
#' @return Data frame with parent-offspring relationships
extract_pedigree_data <- function(pop_history, config_generations) {
  
  cat("=== Extracting Pedigree Information ===\n")
  
  pedigree_list <- list()
  
  # Process all populations
  all_pops <- names(pop_history)
  
  for (pop_name in all_pops) {
    pop <- pop_history[[pop_name]]
    
    # Extract generation number
    if (grepl("reference_gen", pop_name)) {
      gen_num <- as.integer(sub("reference_gen", "", pop_name))
      pop_type <- "reference"
    } else if (pop_name == "validation") {
      gen_num <- config_generations$validation
      pop_type <- "validation"
    } else if (grepl("test_gen", pop_name)) {
      gen_num <- as.integer(sub("test_gen", "", pop_name))
      pop_type <- "test"
    } else {
      next
    }
    
    # Extract pedigree from AlphaSimR object
    n_ind <- nInd(pop)
    
    for (i in 1:n_ind) {
      pedigree_list[[length(pedigree_list) + 1]] <- data.frame(
        individual_id = paste0(pop_type, "_gen", gen_num, "_ind", i),
        generation = gen_num,
        population_type = pop_type,
        mother_id = ifelse(is.null(pop@mother[i]) || is.na(pop@mother[i]), 
                          NA, pop@mother[i]),
        father_id = ifelse(is.null(pop@father[i]) || is.na(pop@father[i]), 
                          NA, pop@father[i]),
        stringsAsFactors = FALSE
      )
    }
  }
  
  pedigree_df <- do.call(rbind, pedigree_list)
  
  cat(sprintf("✓ Pedigree extracted for %d individuals\n", nrow(pedigree_df)))
  cat(sprintf("  - Reference: %d\n", sum(pedigree_df$population_type == "reference")))
  cat(sprintf("  - Validation: %d\n", sum(pedigree_df$population_type == "validation")))
  cat(sprintf("  - Test: %d\n\n", sum(pedigree_df$population_type == "test")))
  
  return(pedigree_df)
}

#' Extract genotype data from population
#' 
#' @param pop AlphaSimR population object
#' @param config
#' @return Data frame with genotype matrix
simulate_genotype_data <- function(pop, SP, config, haplo_raw = NULL) {
  
  cat("=== STEP 3: Simulate Genotype Data ===\n")
  
  # Pull SNP genotypes
  if (is.null(haplo_raw)) {
    haplo_raw <- pullSnpHaplo(pop, simParam = SP) # pullSegSiteHaplo(pop, simParam = SP) 
    cat("  Using SNP chip markers (includes QTLs)\n")
  } else {
    cat("  Using pre-extracted haplotype data for consistency\n")
  }

  # Convert to genotype (sum of two haplotypes)
  n_ind <- nInd(pop)
  geno <- haplo_raw[seq(1, 2*n_ind, by=2), ] + haplo_raw[seq(2, 2*n_ind, by=2), ]
  n_ind <- nrow(geno)
  n_snps <- ncol(geno)
  
  # Validation
  if (n_ind != nInd(pop)) {
    stop(sprintf("ERROR: Expected %d individuals but got %d rows", nInd(pop), n_ind))
  }

  cat(sprintf("Genotype matrix extracted:\n"))
  cat(sprintf("  - Dimensions: %d individuals × %d SNPs\n", n_ind, n_snps))
  cat(sprintf("  - Values: {0, 1, 2} (minor allele count)\n"))
  cat(sprintf("  - Validation: %d expected individuals = %d actual rows ✓\n", nInd(pop), n_ind))
  
  # Add missing data jika dikonfigurasi
  missing_rate <- as.numeric(config$quality_control$missing_data_rate)
  if (missing_rate > 0) {
    n_missing <- round(n_ind * n_snps * missing_rate)
    missing_idx <- sample(length(geno), n_missing)
    geno[missing_idx] <- NA
    
    cat(sprintf("  - Missing data: %.2f%% (%d values)\n", 
                missing_rate * 100, n_missing))
  }
  
  # Convert ke data frame
  geno_df <- as.data.frame(geno)
  colnames(geno_df) <- paste0("snp_", 1:n_snps)
  geno_df <- cbind(individual = paste0("ind_", 1:n_ind), geno_df)
  
  cat("✓ Genotype data created\n\n")
  
  return(list(geno_df = geno_df, haplo_raw = haplo_raw))
}

#' Generate true QTL effects at haplotype level
#' 
#' @param pop AlphaSimR population object
#' @param SP SimParam object
#' @param config 
#' @return Data frame with true QTL effects
simulate_qtl_haplotype_effects <- function(pop, SP, config) {
  
  cat("=== STEP 3B: Simulate True QTL Haplotype Effects ===\n")
  
  # Extract true QTL effects from AlphaSimR
  true_qtl_effects <- SP$traits[[1]]@addEff  # Effects pada additive scale
  
  n_qtl <- length(true_qtl_effects)
  
  cat(sprintf("True QTL effects extracted:\n"))
  cat(sprintf("  - Number of QTL: %d\n", n_qtl))
  cat(sprintf("  - Effect range: [%.4f, %.4f]\n", 
              min(true_qtl_effects), max(true_qtl_effects)))
  cat(sprintf("  - Mean absolute effect: %.4f\n", mean(abs(true_qtl_effects))))
  
  # For the simulation, we assume the haplotype effect is half the genotype effect
  # because in an additive model, the genotype effect = the effect of haplotype A + the effect of haplotype B
  haplotype_effects <- true_qtl_effects / 2
  
  # Create QTL map data frame
  qtl_map <- data.frame(
    qtl_id = paste0("qtl_", 1:n_qtl),
    chromosome = rep(1:config$population$n_chromosomes, 
                    length.out = n_qtl),
    position = sample(1:config$population$chr_length_bp, n_qtl, replace = TRUE),
    true_effect = haplotype_effects,
    effect_size_category = cut(abs(haplotype_effects),
                              breaks = c(0, 0.1, 0.5, 2.0),
                              labels = c("small", "medium", "large"))
  )
  
  cat(sprintf("✓ True QTL haplotype effects simulated\n"))
  cat(sprintf("  - Small effects: %d (|β| < 0.1)\n", 
              sum(qtl_map$effect_size_category == "small")))
  cat(sprintf("  - Medium effects: %d (0.1 ≤ |β| < 0.5)\n", 
              sum(qtl_map$effect_size_category == "medium")))
  cat(sprintf("  - Large effects: %d (|β| ≥ 0.5)\n", 
              sum(qtl_map$effect_size_category == "large")))
  cat("\n")
  
  return(qtl_map)
}


#' Generate phenotype data with QTL effects
#' 
#' @param pop AlphaSimR population object
#' @param SP SimParam object
#' @param config 
#' @return Data frame with phenotypes dan TBV
generate_phenotype_data <- function(pop, SP, config) {
  
  cat("=== STEP 4: Generate QTL Effects & Phenotypes ===\n")
  
  # Extract genetic values (True Breeding Values)
  tbv <- gv(pop)
  
  # Extract phenotypic values
  pheno_vals <- pop@pheno[, 1]
  
  # Calculate realized heritability
  var_g <- var(tbv)
  var_p <- var(pheno_vals)
  realized_h2 <- var_g / var_p
  
  cat(sprintf("Phenotype statistics:\n"))
  cat(sprintf("  - Target h2: %.3f\n", config$genetic_architecture$heritability))
  cat(sprintf("  - Realized h2: %.3f\n", realized_h2))
  cat(sprintf("  - Var(G): %.3f\n", var_g))
  cat(sprintf("  - Var(P): %.3f\n", var_p))
  cat(sprintf("  - Var(E): %.3f\n", var_p - var_g))
  
  # Create data frame
  pheno_df <- data.frame(
    individual = paste0("ind_", 1:nInd(pop)),
    tbv = as.vector(tbv),
    phenotype = as.vector(pheno_vals)
  )
  
  cat("✓ Phenotype data created\n\n")
  
  return(list(
    pheno_df = pheno_df,
    realized_h2 = realized_h2,
    var_g = var_g,
    var_e = var_p - var_g
  ))
}

extract_snp_map <- function(SP, config) {
  
  cat("=== Extracting SNP Map (Genetic-based) ===\n")
  
  map_list <- list()
  global_snp_index <- 1  # Global counter across all chromosomes
  
  for(chr in 1:config$population$n_chromosomes) {
    
    genetic_pos_morgan <- SP$genMap[[chr]]
    genetic_pos_cM <- genetic_pos_morgan * 100
    n_snps <- length(genetic_pos_morgan)
    
    chr_length_bp <- config$population$chr_length_bp
    max_genetic_cM <- max(genetic_pos_cM)
    
    if (max_genetic_cM > 0) {
      physical_pos_bp <- round((genetic_pos_cM / max_genetic_cM) * chr_length_bp)
    } else {
      physical_pos_bp <- round(seq(0, chr_length_bp, length.out = n_snps))
    }
    
    physical_pos_bp <- pmin(physical_pos_bp, chr_length_bp)
    physical_pos_bp <- pmax(physical_pos_bp, 0)
    physical_pos_bp <- as.integer(physical_pos_bp)
    
    snp_ids <- paste0("snp_", global_snp_index:(global_snp_index + n_snps - 1))
    global_snp_index <- global_snp_index + n_snps
    
    chr_map <- data.frame(
      SNPID = snp_ids,
      Chr = chr,
      Position = physical_pos_bp,
      GeneticPos_cM = round(genetic_pos_cM, 4),
      stringsAsFactors = FALSE
    )
    
    map_list[[chr]] <- chr_map
    
    cat(sprintf("  Chr %d: %d SNPs (genetic: 0-%.2f cM, physical: 0-%d bp)\n", 
                chr, n_snps, max(genetic_pos_cM), max(physical_pos_bp)))
  }
  
  map_df <- do.call(rbind, map_list)
  rownames(map_df) <- NULL
  
  cat(sprintf("✓ Total SNPs: %d\n", nrow(map_df)))
  cat("  Note: SNP IDs match genotype CSV column names (snp_1, snp_2, ...)\n\n")
  
  return(map_df)
}

#' Save simulation parameters to JSON
#' 
#' @param config 
#' @param realized_h2 Realized heritability
#' @param output_dir Output directory
save_simulation_parameters <- function(config, realized_h2, output_dir) {
  
  cat("=== STEP 5: Save Simulation Parameters ===\n")
  
  # Compile parameters
  params <- list(
    population = list(
      n_individuals = config$population$n_individuals,
      n_chromosomes = config$population$n_chromosomes,
      n_snps = config$population$n_chromosomes * config$population$n_snps_per_chrom,
      effective_pop_size = config$population$effective_pop_size
    ),
    breeding = list(
      burn_in_generations = config$breeding$burn_in_generations,
      reference_generations = config$breeding$reference_generations,
      validation_generation = config$breeding$validation_generation,
      test_generations = config$breeding$test_generations,
      selection_intensity = config$breeding$selection_intensity,
      mating_design = config$breeding$mating_design
    ),
    genetic_architecture = list(
      n_qtl = config$genetic_architecture$n_qtl,
      heritability = config$genetic_architecture$heritability,
      realized_heritability = realized_h2,
      qtl_effect_distribution = config$genetic_architecture$qtl_effect_distribution,
      qtl_effect_sizes = config$genetic_architecture$qtl_effect_sizes,
      qtl_proportions = config$genetic_architecture$qtl_proportions
    ),
    quality_control = config$quality_control,
    timestamp = as.character(Sys.time())
  )
  
  # Save to JSON
  json_file <- file.path(output_dir, "simulation_parameters.json")
  write_json(params, json_file, auto_unbox = TRUE, pretty = TRUE)
  
  cat(sprintf("✓ Parameters saved to: %s\n\n", json_file))
  
  return(params)
}

#===============================================================================
# MAIN EXECUTION
#===============================================================================

main <- function() {
  
  cat("\n")
  cat("========================================\n")
  cat("  AlphaSimR Genetic Simulation\n")
  cat("========================================\n")
  cat("\n")
  
  # Timestamp awal
  start_time <- Sys.time()
  
  # Load configuration
  config <- load_config("config.yaml")

  # validate parameters
  if (config$population$n_individuals != config$breeding$population_size) {
    warning(sprintf(
      "Inconsistency: founder population size (%d) != final population size (%d). Adjusting...",
      config$population$n_individuals, config$breeding$population_size
    ))
    config$breeding$population_size <- config$population$n_individuals
  }
  
  # Validation of n_qtl is divisible by n_chromosomes
  if (config$genetic_architecture$n_qtl %% config$population$n_chromosomes != 0) {
    adjusted_qtl <- config$genetic_architecture$n_qtl - 
      (config$genetic_architecture$n_qtl %% config$population$n_chromosomes)
    warning(sprintf(
      "n_qtl (%d) not divisible by n_chromosomes (%d). Using %d QTL.",
      config$genetic_architecture$n_qtl, config$population$n_chromosomes, adjusted_qtl
    ))
    config$genetic_architecture$n_qtl <- adjusted_qtl
  }

  print_config_summary(config)
  
  # Create output directory
  output_dir <- config$output$output_dir
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Run simulation pipeline
  founderPop <- initialize_founder_population(config)
  
  breeding_result <- run_breeding_program(founderPop, config)
  pop_history <- breeding_result$pop_history
  SP <- breeding_result$SP
  config_gens <- breeding_result$config_generations

  # Extract pedigree
  pedigree_df <- extract_pedigree_data(pop_history, config_gens)
  
  # Simulate QTL effects
  first_pop <- pop_history[[1]]
  qtl_map <- simulate_qtl_haplotype_effects(first_pop, SP, config)
  
  # Save output files
  cat("=== Saving Output Files ===\n")

  # Extract and save SNP map
  snp_map <- extract_snp_map(SP, config)
  old_scipen <- getOption("scipen", default = 0)
  options(scipen = 999)
  write.table(snp_map, 
              file.path(output_dir, "map.txt"), 
              row.names = FALSE, 
              col.names = TRUE,
              quote = FALSE,
              sep = "\t") 
  
  options(scipen = old_scipen)

  cat(sprintf("✓ SNP map saved: %d SNPs across %d chromosomes\n", 
              nrow(snp_map), config$population$n_chromosomes))

  # Initialize combined dataframes
  geno_ref_combined <- NULL
  pheno_ref_combined <- NULL
  
  # Process reference populations
  cat("\n--- Reference Populations (Training) ---\n")
  for (gen in config_gens$reference) {
    pop_name <- paste0("reference_gen", gen)
    pop <- pop_history[[pop_name]]
    
    cat(sprintf("Processing generation %d...\n", gen))
    
    # Extract haplotypes once, reuse for genotype
    geno_result <- simulate_genotype_data(pop, SP, config, haplo_raw = NULL)
    geno_df <- geno_result$geno_df
    
    geno_df$generation <- gen
    geno_df$population_type <- "reference"
    geno_df$individual_id <- paste0("reference_gen", gen, "_ind", 1:nrow(geno_df))
    
    pheno_result <- generate_phenotype_data(pop, SP, config)
    pheno_df <- pheno_result$pheno_df
    pheno_df$generation <- gen
    pheno_df$population_type <- "reference"
    pheno_df$individual_id <- paste0("reference_gen", gen, "_ind", 1:nrow(pheno_df))
    
    # Combine
    geno_ref_combined <- rbind(geno_ref_combined, geno_df)
    pheno_ref_combined <- rbind(pheno_ref_combined, pheno_df)
  }
  
  # Save combined reference data
  write.csv(geno_ref_combined, file.path(output_dir, "genotype_reference.csv"), 
            row.names = FALSE)
  write.csv(pheno_ref_combined, file.path(output_dir, "phenotype_reference.csv"), 
            row.names = FALSE)
  
  cat(sprintf("✓ Reference data saved (combined %d generations, %d individuals)\n", 
              length(config_gens$reference), nrow(geno_ref_combined)))
  
  # Process validation population
  cat("\n--- Validation Population (Hyperparameter Tuning) ---\n")
  pop_val <- pop_history[["validation"]]
  
  # Extract haplotypes ONCE, reuse for genotype
  geno_val_result <- simulate_genotype_data(pop_val, SP, config, haplo_raw = NULL)
  geno_val <- geno_val_result$geno_df
  
  geno_val$generation <- config_gens$validation
  geno_val$population_type <- "validation"
  geno_val$individual_id <- paste0("validation_gen", config_gens$validation, 
                                   "_ind", 1:nrow(geno_val))
  
  pheno_val_result <- generate_phenotype_data(pop_val, SP, config)
  pheno_val <- pheno_val_result$pheno_df
  pheno_val$generation <- config_gens$validation
  pheno_val$population_type <- "validation"
  pheno_val$individual_id <- paste0("validation_gen", config_gens$validation, 
                                    "_ind", 1:nrow(pheno_val))
  
  write.csv(geno_val, file.path(output_dir, "genotype_validation.csv"), 
            row.names = FALSE)
  write.csv(pheno_val, file.path(output_dir, "phenotype_validation.csv"), 
            row.names = FALSE)
  
  cat(sprintf("✓ Validation data saved (generation %d, %d individuals)\n", 
              config_gens$validation, nrow(geno_val)))
  
  # Process test populations
  cat("\n--- Test Populations (Accuracy Evaluation) ---\n")
  for (gen in config_gens$test) {
    pop_name <- paste0("test_gen", gen)
    pop_test <- pop_history[[pop_name]]
    
    cat(sprintf("Processing test generation %d...\n", gen))
    
    # CRITICAL: Extract haplotypes ONCE, reuse for genotype
    geno_test_result <- simulate_genotype_data(pop_test, SP, config, haplo_raw = NULL)
    geno_test <- geno_test_result$geno_df
    
    geno_test$generation <- gen
    geno_test$population_type <- "test"
    geno_test$individual_id <- paste0("test_gen", gen, "_ind", 1:nrow(geno_test))
    
    pheno_test_result <- generate_phenotype_data(pop_test, SP, config)
    pheno_test <- pheno_test_result$pheno_df
    pheno_test$generation <- gen
    pheno_test$population_type <- "test"
    pheno_test$individual_id <- paste0("test_gen", gen, "_ind", 1:nrow(pheno_test))
    
    write.csv(geno_test, 
              file.path(output_dir, paste0("genotype_test_gen", gen, ".csv")), 
              row.names = FALSE)
    write.csv(pheno_test, 
              file.path(output_dir, paste0("phenotype_test_gen", gen, ".csv")), 
              row.names = FALSE)
    
    cat(sprintf("  ✓ Test generation %d saved (%d individuals)\n", 
                gen, nrow(geno_test)))
  }
  
  # Save QTL map and pedigree
  write.csv(qtl_map, file.path(output_dir, "true_qtl_effects.csv"), 
            row.names = FALSE)
  write.csv(pedigree_df, file.path(output_dir, "pedigree.csv"), 
            row.names = FALSE)
  
  cat("\n✓ QTL effects saved\n")
  cat("✓ Pedigree data saved\n")
  
  # Save parameters dengan realized h2 dari reference
  realized_h2_ref <- var(pheno_ref_combined$tbv) / var(pheno_ref_combined$phenotype)
  save_simulation_parameters(config, realized_h2_ref, output_dir)
  cat(sprintf("  Columns: individual, tbv, phenotype, generation\n"))
  
  # Completion message
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "secs")
  
  cat("\n")
  cat("========================================\n")
  cat("  Simulation Completed Successfully!\n")
  cat("========================================\n")
  cat(sprintf("\nElapsed time: %.1f seconds\n", elapsed))
  cat(sprintf("Output directory: %s/\n", output_dir))
  cat("\nOutput files:\n")
  cat("  - genotype_reference.csv (training, multi-generational)\n")
  cat("  - phenotype_reference.csv (training)\n")
  cat("  - genotype_validation.csv (hyperparameter tuning)\n")
  cat("  - phenotype_validation.csv\n")
  cat("  - genotype_test.csv\n")
  cat("  - phenotype_test\n")
  cat("  - pedigree.csv (parent-offspring relationships)\n")
  cat("  - true_qtl_effects.csv\n")
  cat("  - simulation_parameters.json\n")
}

# Run main function
if (!interactive()) {
  main()
}
