# Neo-XY Population Simulation (`gen3_co_SynAdj_population_ondemand.py`)

This directory hosts the neo-sex–chromosome simulator used to explore recombination, mutation, and fixation dynamics in a bottlenecked population of worms. The active script performs on-demand meiosis so every inherited chromosome triggers its own crossover (sons always retain the neo-Y left arm; daughters inherit the male X plus a female X).

## Requirements

```bash
python -m pip install numpy scipy plotly Pillow
```

`scipy` is only needed if you edit the map generator; `plotly` and `Pillow` are for the HTML/PNG outputs.

## Running the simulation

```bash
python gen3_co_SynAdj_population_ondemand.py
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
2. Each generation, adds unique mutations to every chromosome, performs on-demand crossovers using the specified maps, and resamples parents with replacement to form the next generation.
3. Records haplotype histories, cumulative Y–X SNP divergence, and whether any neo-Y locus is fixed for a derived allele across all sampled males.
4. Repeats steps 1–3 for `iterations` independent runs and aggregates the outputs.

## Outputs (written to `MAsim/`)

| File | Description |
| --- | --- |
| `adj_snpimg.html` | Heatmap of haplotypes over all generations (Y, X1, X2, X3 stacked with spacers). |
| `adj_snpplot.html` | Scatter plot of the terminal-generation chromatids (offset by +0…+3). |
| `adj_snpaccum.html` | Scatter of per-locus SNP divergence (counts loci where Y≠X and at least one allele is derived). |
| `adj_snpaccum_image.html` | Heatmap of divergence counts across generations. |
| `adj_snp_accum.txt` | Raw divergence matrix (`generations × chromosome_size`). |
| `A_plus_B_rgb.png` | RGB overlay highlighting loci where Y (red) and X1 (green) exceed value 10. |
| `adj_population.vcf` | Haploid VCF of the final population (male Y, male X, female X2, female X3). |
| `neoY_fixed_mutation_frequency.txt` | Fraction of iterations where each neo-Y locus was fixed for the same derived allele across all sampled males. |

## Notes

- Mutation IDs start at 4 because founders use labels 0–3.
- SNP divergence currently counts loci where Y and X differ and at least one allele is derived; tweak `variant_mask` in `simulate_population()` if you want a different definition.
- `neoY_fixed_mutation_frequency.txt` checks both “all derived” and “all identical”, so only true fixation events contribute.
