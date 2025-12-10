import numpy as np
import plotly.express as px

from genomic_position_generator import load_map_data, GenomicPositionGenerator


class Copos:
    def __init__(self, map_file: str) -> None:
        map_data = load_map_data(map_file)
        self.generator = GenomicPositionGenerator(map_data)

    def return_copos(self, length: int) -> int:
        pos = int(length * self.generator.generate_position())
        return max(1, min(length - 1, pos))


def recombine_pair(src_a: np.ndarray, src_b: np.ndarray, pos: int) -> tuple[np.ndarray, np.ndarray]:
    child_a = np.empty_like(src_a)
    child_b = np.empty_like(src_b)
    if pos <= 0 or pos >= src_a.size:
        np.copyto(child_a, src_a)
        np.copyto(child_b, src_b)
        return child_a, child_b

    child_a[:pos] = src_a[:pos]
    child_a[pos:] = src_b[pos:]
    child_b[:pos] = src_b[:pos]
    child_b[pos:] = src_a[pos:]
    return child_a, child_b


def mutate_parents(
    arrays: tuple[np.ndarray, ...],
    mutation_counter: int,
    rng: np.random.Generator,
    chromosome_size: int,
    mutations_per_chromosome: int,
) -> int:
    for arr in arrays:
        idx = rng.integers(1, chromosome_size - 1, size=(arr.shape[0], mutations_per_chromosome))
        for row, positions in enumerate(idx):
            new_ids = np.arange(mutation_counter, mutation_counter + mutations_per_chromosome, dtype=np.int64)
            arr[row, positions] = new_ids
            mutation_counter += mutations_per_chromosome
    return mutation_counter


def inherit_chromatid(
    chrom_a: np.ndarray,
    chrom_b: np.ndarray,
    map_source: Copos,
    rng: np.random.Generator,
    chromosome_size: int,
    left_parent_only: bool = False,
) -> np.ndarray:
    """Copy one of the chromatids from a meiosis triggered on demand.

    If left_parent_only is True, limit choices to chromatids that retain the left flank from chrom_a
    (i.e., chrom_a itself and the recombinant receiving chrom_a's left arm)."""
    pos = map_source.return_copos(chromosome_size)
    rec_a, rec_b = recombine_pair(chrom_a, chrom_b, pos)
    if left_parent_only:
        choices = (chrom_a, rec_a)
    else:
        choices = (chrom_a, chrom_b, rec_a, rec_b)
    selected = choices[rng.integers(len(choices))]
    out = np.empty_like(selected)
    np.copyto(out, selected)
    return out


def simulate_generation(
    male_y: np.ndarray,
    male_x: np.ndarray,
    female_x2: np.ndarray,
    female_x3: np.ndarray,
    male_map: Copos,
    female_map: Copos,
    rng: np.random.Generator,
    chromosome_size: int,
    mutation_counter: int,
    mutations_per_chromosome: int,
) -> tuple[int, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    male_count = male_y.shape[0]
    female_count = female_x2.shape[0]

    mutation_counter = mutate_parents(
        (male_y, male_x, female_x2, female_x3),
        mutation_counter,
        rng,
        chromosome_size,
        mutations_per_chromosome,
    )

    next_male_y = np.empty_like(male_y)
    next_male_x = np.empty_like(male_x)
    next_female_x2 = np.empty_like(female_x2)
    next_female_x3 = np.empty_like(female_x3)

    male_parents = rng.integers(male_count, size=male_count)
    female_parents_for_males = rng.integers(female_count, size=male_count)
    male_parents_for_females = rng.integers(male_count, size=female_count)
    female_parents_for_females = rng.integers(female_count, size=female_count)

    for idx in range(male_count):
        mp = male_parents[idx]
        fp = female_parents_for_males[idx]

        np.copyto(
            next_male_y[idx],
            inherit_chromatid(
                male_y[mp],
                male_x[mp],
                male_map,
                rng,
                chromosome_size,
                left_parent_only=True,
            ),
        )

        female_chrom_a = female_x2[fp]
        female_chrom_b = female_x3[fp]
        np.copyto(
            next_male_x[idx],
            inherit_chromatid(female_chrom_a, female_chrom_b, female_map, rng, chromosome_size),
        )

    for idx in range(female_count):
        mp = male_parents_for_females[idx]
        fp = female_parents_for_females[idx]

        np.copyto(
            next_female_x2[idx],
            inherit_chromatid(
                male_x[mp],
                male_y[mp],
                male_map,
                rng,
                chromosome_size,
                left_parent_only=True,
            ),
        )

        female_chrom_a = female_x2[fp]
        female_chrom_b = female_x3[fp]
        np.copyto(
            next_female_x3[idx],
            inherit_chromatid(female_chrom_a, female_chrom_b, female_map, rng, chromosome_size),
        )

    return mutation_counter, next_male_y, next_male_x, next_female_x2, next_female_x3


def simulate_population(
    generations: int,
    chromosome_size: int,
    num_pairs: int,
    male_map: Copos,
    female_map: Copos,
    rng: np.random.Generator,
    mutations_per_chromosome: int,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    male_y = np.zeros((num_pairs, chromosome_size), dtype=np.int64)
    male_x = np.ones((num_pairs, chromosome_size), dtype=np.int64)
    female_x2 = 2 * np.ones((num_pairs, chromosome_size), dtype=np.int64)
    female_x3 = 3 * np.ones((num_pairs, chromosome_size), dtype=np.int64)

    yac = np.empty((generations, chromosome_size), dtype=np.int64)
    x1ac = np.empty((generations, chromosome_size), dtype=np.int64)
    x2ac = np.empty((generations, chromosome_size), dtype=np.int64)
    x3ac = np.empty((generations, chromosome_size), dtype=np.int64)
    q = np.zeros((generations, chromosome_size), dtype=np.float64)

    mutation_counter = 4

    for gen in range(generations):
        yac[gen] = male_y[0]
        x1ac[gen] = male_x[0]
        x2ac[gen] = female_x2[0]
        x3ac[gen] = female_x3[0]

        variant_mask = (male_y != male_x) & ((male_y > 3) | (male_x > 3))
        q[gen] += np.sum(variant_mask, axis=0)

        mutation_counter, male_y, male_x, female_x2, female_x3 = simulate_generation(
            male_y,
            male_x,
            female_x2,
            female_x3,
            male_map,
            female_map,
            rng,
            chromosome_size,
            mutation_counter,
            mutations_per_chromosome,
        )

    return yac, x1ac, x2ac, x3ac, q, male_y, male_x, female_x2, female_x3


def plot_results(yac: np.ndarray, x1ac: np.ndarray, x2ac: np.ndarray, x3ac: np.ndarray, q: np.ndarray) -> None:
    from PIL import Image

    spacer = np.full((10, yac.shape[1]), -2.0)
    outimg = np.concatenate((yac, spacer, x1ac, spacer, x2ac, spacer, x3ac), axis=0)
    colorscale = ["white"] + px.colors.sequential.Bluered
    px.imshow(outimg, aspect="auto", color_continuous_scale=colorscale).write_html("adj_snpimg.html")

    outplt = np.vstack((yac[-1], 1 + x1ac[-1], 2 + x2ac[-1], 3 + x3ac[-1]))
    px.scatter(outplt.T).write_html("adj_snpplot.html")

    px.scatter(q[-1]).write_html("adj_snpaccum.html")
    np.savetxt("adj_snp_accum.txt", q[-1], "%f")

    px.imshow(q, aspect="auto").write_html("adj_snpaccum_image.html")
    maskA = (yac > 10)
    maskB = (x1ac > 10)
    h, w = yac.shape
    rgb = np.zeros((h, w, 3), dtype=np.uint8)
    rgb[..., 0] = maskA.astype(np.uint8) * 255
    rgb[..., 1] = maskB.astype(np.uint8) * 255
    Image.fromarray(rgb).save("A_plus_B_rgb.png")


def write_population_vcf(
    path: str,
    chromosome_size: int,
    male_y: np.ndarray,
    male_x: np.ndarray,
    female_x2: np.ndarray,
    female_x3: np.ndarray,
) -> None:
    samples = (
        [f"maleY_{i}" for i in range(male_y.shape[0])]
        + [f"maleX_{i}" for i in range(male_x.shape[0])]
        + [f"femaleX2_{i}" for i in range(female_x2.shape[0])]
        + [f"femaleX3_{i}" for i in range(female_x3.shape[0])]
    )
    genotype_matrix = np.vstack((male_y, male_x, female_x2, female_x3))

    with open(path, "w", encoding="utf-8") as fh:
        fh.write("##fileformat=VCFv4.3\n")
        fh.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Haploid genotype based on allele IDs\">\n")
        header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *samples]
        fh.write("\t".join(header) + "\n")

        for pos in range(chromosome_size):
            alleles = genotype_matrix[:, pos]
            unique_ids = np.unique(alleles)
            ref_id = unique_ids[0]
            ref_allele = f"A{ref_id}"
            alt_ids = unique_ids[1:]
            alt_alleles = ",".join(f"A{aid}" for aid in alt_ids) if alt_ids.size else "."

            allele_lookup = {ref_id: "0"}
            for idx, aid in enumerate(alt_ids, start=1):
                allele_lookup[aid] = str(idx)

            genotypes = [allele_lookup[a] for a in alleles]
            line = (
                f"neoXY\t{pos + 1}\t.\t{ref_allele}\t{alt_alleles}\t.\tPASS\t.\tGT\t"
                + "\t".join(genotypes)
            )
            fh.write(line + "\n")


def main() -> None:
    iterations = 6
    generations = 27
    chromosome_size = 10_000
    num_pairs = 5
    mutations_per_chromosome = 1

    rng = np.random.default_rng()
    male_map = Copos("ivsqrt.txt")
    female_map = Copos("ivsqrt.txt")
    #female_map = Copos("flat.txt")
#    male_map = Copos("flat.txt")
#    female_map = Copos("ivs.txt")
#    male_map = Copos("ivs.txt")

    q_accum = np.zeros((generations, chromosome_size), dtype=np.float64)
    yac = np.empty((generations, chromosome_size), dtype=np.float64)
    x1ac = np.empty((generations, chromosome_size), dtype=np.float64)
    x2ac = np.empty((generations, chromosome_size), dtype=np.float64)
    x3ac = np.empty((generations, chromosome_size), dtype=np.float64)
    mutation_counts = np.zeros(chromosome_size, dtype=np.float64)

    final_male_y = final_male_x = final_female_x2 = final_female_x3 = None

    for _ in range(iterations):
        y_hist, x1_hist, x2_hist, x3_hist, q, male_y, male_x, female_x2, female_x3 = simulate_population(
            generations,
            chromosome_size,
            num_pairs,
            male_map,
            female_map,
            rng,
            mutations_per_chromosome,
        )
        yac = y_hist.astype(float)
        x1ac = x1_hist.astype(float)
        x2ac = x2_hist.astype(float)
        x3ac = x3_hist.astype(float)
        q_accum += q
        final_male_y = male_y
        final_male_x = male_x
        final_female_x2 = female_x2
        final_female_x3 = female_x3
        fixed_mask = (male_y[0] > 3) & np.all(male_y == male_y[0], axis=0)
        mutation_counts += fixed_mask.astype(np.float64)

    plot_results(yac, x1ac, x2ac, x3ac, q_accum)
    if final_male_y is not None:
        write_population_vcf(
            "adj_population.vcf",
            chromosome_size,
            final_male_y,
            final_male_x,
            final_female_x2,
            final_female_x3,
        )
        np.savetxt("neoY_fixed_mutation_frequency.txt", (mutation_counts / iterations)[None, :], fmt="%.6f")


if __name__ == "__main__":
    main()
