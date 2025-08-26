# Founder-effect-array
Analysis potential Founder effect through array sequencing data.

## Introduction
When analyzing genetic data, we often observe a significant difference in the allele frequency of certain pathogenic variants between disease and control populations. This observation can be an initial indicator of a founder effect, a genetic phenomenon where a particular allele is highly prevalent in a population due to its presence in a limited number of ancestral individuals. A strong **founder effect** can significantly accelerate the prevalence of a genetic disorder within a specific community.

To provide genetic evidence supporting a founder effect, researchers often employ methods such as homozygous mapping or linkage disequilibrium analysis. These approaches investigate the shared genetic segments (haplotypes) surrounding the variant of interest, which are inherited from a common ancestor. In our study, we focused on a selection of common, suspected pathogenic variants identified in the Inherited Retinal Degeneration (IRD) Project. Our primary objective was to determine whether these variants are associated with a founder effect within the Taiwanese population.

To achieve this, we developed a bioinformatics pipeline to **analyze the haplotype frequency of these variants**. We compared patients from the IRD database (case group) with a control group comprising individuals from the Taiwan Biobank and unaffected family members from the IRD cohort. This pipeline processes array genotype data (.bed, .bim, .fam) with standard quality control filtering, for example:

```
plink --geno 0.05 --mind 0.05 --maf 0.05 --hwe 1e-6 --snps-only just-acgt
```

The core of our analysis involves the following key steps:

**Haplotype Generation**: We first define the flanking genomic regions around each target variant. Using PLINK (**b_1.output_region.sh**), we extract the selected regions. Then, we use BEAGLE (**b_2.phased.sh**) to infer the haplotype phases, which are essential for subsequent comparisons.

**Phenotype Integratio**n: We integrate phenotype information into the family files (**c.phenotype_fam.py**) to clearly distinguish between case and control groups for downstream analysis.

**Haplotype Analysis**: Finally, we perform a comprehensive haplotype analysis (**d.final_haplotype.py**). We calculate the haplotype frequency differences between the case and control groups and determine their statistical significance using a permutation test (100,000 times). This rigorous statistical approach minimizes the risk of false positives and provides robust evidence to support or refute the presence of a founder effect.

The results of this analysis will provide crucial insights into the genetic architecture of IRD in Taiwan and can inform future genetic screening and diagnostic strategies.

## Usage

### Scripts

a.output_rsnumber.py: Generates flanking regions based on a given variant position. (Requires specifying the array prefix, chromosome, position, rs number, and desired flanking region SNPs.)
b_1.output_region.sh: Uses PLINK to extract the target region generated in the previous step.
b_2.phased.sh: Performs phasing using Beagle.
c.phenotype_fam.py: Adds phenotype information to the .fam file.
d.final_haplotype.py: Calculates haplotype frequency differences between case and control groups, including p-values and a permutation test (100,000 iterations).

```
Step1: python3 a.output_rsnumber.py
Step2: sh b_1.output_region.sh
Step3: sh b_2.phased.sh
Step4: python3 c.phenotype_fam.py
Step5: python3 d.final_haplotype.py
```

Author: Chien-Yu, Lin
