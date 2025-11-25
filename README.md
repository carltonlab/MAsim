# Neo-XY Population Simulation (`gen3_co_SynAdj_population.py`)

This directory hosts the current neo-sex–chromosome simulator used to explore recombination, mutation, and fixation dynamics in a bottlenecked population of worms.

## Requirements

```bash
python -m pip install numpy scipy plotly Pillow
```

`scipy` is only needed if you edit the map generator; `plotly` and `Pillow` are for the HTML/PNG outputs.

## Running the simulation

```bash
cd MAsim
python gen3_co_SynAdj_population.py
```

Key parameters (edit inside `main()`):

- `iterations`: independent runs to average (default 100).
- `generations`: pedigree length per run (default 100).
- `chromosome_size`: loci per chromosome (default 10 000).
- `num_pairs`: number of breeding males/females per generation (default 5).
- `mutations_per_chromosome`: new mutations per chromosome per generation (default 3).
- `male_map`, `female_map`: choose crossover landscapes (`ivsqrt.txt`, `ivs.txt`, or `flat.txt`).

## What the script does

1. Initializes `num_pairs` males and females with labeled neo-Y/neo-X chromosomes.
2. Each generation, adds unique mutations to every chromosome, draws crossover positions from the specified maps, and resamples parents with replacement to form the next generation.
3. Records haplotype histories, cumulative Y–X SNP divergence, and whether any neo-Y locus is fixed for a derived allele across all sampled males.
4. Repeats steps 1–3 for `iterations` independent runs and aggregates the outputs.

## Outputs (written to `MAsim/`)

| File | Description |
| --- | --- |
| `adj_snpimg.html` | Heatmap of haplotypes over all generations (Y, X1, X2, X3 stacked). |
| `adj_snpplot.html` | Scatter of the final generation’s four chromatids. |
| `adj_snpaccum.html` | Per-locus SNP divergence (counts Y≠X with ≥1 derived allele). |
| `adj_snpaccum_image.html` | Heatmap of the divergence counts over generations. |
| `adj_snp_accum.txt` | Raw divergence values (`generations × chromosome_size`). |
| `A_plus_B_rgb.png` | RGB overlay showing where Y (red) and X1 (green) carry alleles >10. |
| `adj_population.vcf` | Haploid VCF of the final population (all male/female chromosomes). |
| `neoY_fixed_mutation_frequency.txt` | Fraction of runs where each locus was fixed for the same derived allele across all neo-Y copies. |

## Notes

- Mutation IDs start at 4 because founders use labels 0–3.
- SNP divergence currently counts loci where Y and X differ and at least one allele is derived; tweak `variant_mask` in `simulate_population()` if you want a different definition.
- `neoY_fixed_mutation_frequency.txt` checks both “all derived” and “all identical,” so only true fixation events contribute.

Feel free to adjust map files or parameters and rerun; the outputs listed above will always reflect the most recent invocation.
