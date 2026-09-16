# eccmr: Error-Correcting Code Inspired Robust Mendelian Randomization

`eccmr` implements **ECC-MR**, an error-correcting code inspired framework
for robust Mendelian randomization with GWAS summary statistics.

**Core idea.** Horizontal pleiotropy is treated as *sparse, correctable
corruption* of SNP-level causal signals, and linkage disequilibrium (LD)
among instruments is exploited as *structured redundancy* rather than a
nuisance:

- A **parity-check matrix** `H` is built from the LD graph so that the
  causal signal lies in its null space: `H %*% beta_X = 0` by construction.
- The **syndrome** `s = H %*% beta_Y = H alpha + H epsilon` isolates the
  pleiotropic errors — the causal component is annihilated.
- **Sparse syndrome decoding** (inverse-variance weighted least squares +
  L1 sparsity + syndrome-consistency penalty, jointly convex in
  `(theta, alpha)`) estimates the causal effect and the pleiotropic effects
  *jointly*, correcting — not discarding — pleiotropic instruments.
- A **syndrome-based global test of pleiotropy** (`T_syn`, chi-squared)
  generalizes the MR-PRESSO global test to correlated instruments.
- Tuning parameters are selected by a **BIC criterion** along the penalty
  path (LD-block cross-validation is available as an option), and the
  selected support is **refitted without the L1 penalty** (relaxed refit)
  to remove shrinkage bias from the causal estimate.
- Standard errors and confidence intervals come from an **LD-block
  bootstrap**.

Because correlated instruments are retained, ECC-MR is more efficient than
methods that require aggressive LD pruning, and its correction capacity
grows with the amount of LD redundancy.

## Installation

From a local copy of this directory:

```r
# install.packages("remotes")
remotes::install_local("ECCMR")          # path to the package directory
```

or, once the repository is on GitHub:

```r
remotes::install_github("DrWoodWood/ECCMR")
```

The only hard dependency besides base R is the recommended package
`Matrix`.

## Quick start

```r
library(eccmr)

## Simulated example (scenario A: directional pleiotropy)
sim <- simulate_eccmr(n_snps = 300, block_size = 10, theta = 0.3,
                      scenario = "A", seed = 42)

fit <- eccmr(sim$beta_X, sim$beta_Y, sim$se_Y,
             blocks = sim$blocks,        # or ld_mat = <LD correlation matrix>
             n_boot = 500, seed = 1)

print(fit)      # causal estimate, SE, 95% CI, syndrome test
summary(fit)    # + top pleiotropic instruments
plot(fit)       # standardized syndrome diagnostics
```

With real data (e.g. harmonized `TwoSampleMR` output), supply:

- `beta_X`, `beta_Y`: exposure and outcome association estimates;
- `se_Y`: standard errors of the outcome associations;
- `blocks`: LD-block membership (e.g. Berisa & Pickrell blocks), **or**
  `ld_mat`: a pairwise LD correlation matrix from a matched-ancestry
  reference panel — correlated instruments are welcome, no LD pruning is
  needed.

## Benchmarking against mainstream MR methods

The package ships reference implementations of the standard summary-data MR
estimators and a simulation engine that benchmarks ECC-MR against them.
Conventional methods are applied to one LD-pruned representative per block
(standard clumping practice); ECC-MR keeps all instruments — which is
exactly where its efficiency advantage comes from.

```r
## One data set, five methods
sim <- simulate_eccmr(n_snps = 300, block_size = 10, theta = 0.3,
                      scenario = "A", directional = 0.05, seed = 42)
compare_mr_methods(sim$beta_X, sim$beta_Y, sim$se_Y, blocks = sim$blocks,
                   n_boot = 200, seed = 1)

## Full simulation benchmark (bias / RMSE / coverage / rejection)
bench <- run_simulation(nsim = 50, scenario = "A", directional = 0.05,
                        seed = 1)          # cores = 4 for parallel replicates
print(bench)    # comparison table, lowest RMSE marked with *
plot(bench)     # bias and RMSE panels

## Scenario C (LD-clustered, weak directional pleiotropy) is where the
## ECC correction capacity is most visible:
run_simulation(nsim = 50, scenario = "C", pi_block = 0.3, seed = 1)
```

For production-grade versions of the comparator methods we recommend
`TwoSampleMR`, `MendelianRandomization` and `MRPRESSO`; the built-in
implementations exist to keep the benchmark dependency-free.

## Package layout

| Function | Purpose |
|---|---|
| `eccmr()` | Main entry point: decode, tune, infer |
| `partition_ld_blocks()` | Derive LD blocks from an LD matrix |
| `build_parity_check()` | Construct the parity-check matrix `H` |
| `syndrome_test()` | Syndrome + global pleiotropy test |
| `ecc_fit()` | Core convex solver (fixed tuning parameters) |
| `simulate_eccmr()` | Synthetic data (manuscript scenarios A/B/C) |
| `compare_mr_methods()` | ECC-MR vs mainstream methods, one data set |
| `run_simulation()` | Replicated simulation benchmark + summary/plot |
| `mr_ivw()`, `mr_egger()`, `mr_weighted_median()`, `mr_lasso()` | Comparator methods |

## Method reference

Jiang J, et al. *ECC-MR: an error-correcting-code framework for robust Mendelian randomization and pleiotropy decoding on correlated instruments.* Manuscript.

## License

MIT © Mu-BioDig Group
