from alphagenome.data import genome
from alphagenome.models import dna_client
import numpy as np
import pandas as pd

# A genomic interval is specified using `genome.Interval`:
interval = genome.Interval(chromosome='chr1', start = 1_000, end = 1_010) #Default hg38

## Resizing
# Interval width and center
print(f"Interval is: {interval}")
print(f"Center of interval is: {interval.center()}")
print(f"Width of interval is: {interval.width}")

# Comparing Intervals

second_interval = genome.Interval(chromosome="chr1", start=1_005, end=1_015)
print(f"Second inteval is: {second_interval}")
print(f"Do they overlap?: {interval.overlaps(second_interval)}") #=True
print(f"Does interval1 contain interval2?: {interval.contains(second_interval)}") #=False
print(f"Intersection interval is: {interval.intersect(second_interval)}")

## Genomic Variants
print("------------------")
print("---- VARIANTS ----")
print("------------------")
# A genome.Variant specifies a genetic variant:
variant = genome.Variant(
        chromosome='chr3', position=10000, reference_bases='A', alternate_bases='C'
        )
print(f"A SNP: {variant}")
# This variant changes A -> at position 10,000 on chr3. 

# INSERTION:
variant = genome.Variant(
        chromosome='chr3',
        position=10000,
        reference_bases='T',
        alternate_bases="CGTCAAT",
        )
print(f"An Insertion: {variant}")
# DELETION:
variant = genome.Variant(
        chromosome='chr3',
        position=10_000,
        reference_bases='AGGGATC',
        alternate_bases='C',
        )
print(f"A Deletion: {variant}")


## Reference Interval
# Can we get the genome.Interval corresponding to the reference bases of the variant using genome.Variant.reference_interval
variant = genome.Variant(
        chromosome='chr3', position=10_000, reference_bases='A', alternate_bases='T'
        )

print(f"\nBy defining the variant: {variant}")
print(f"We can easily grab the reference interval: {variant.reference_interval}")


# A common use-case is to make predictions in a genome region around a variant, which involves resizing the `genome.Variant.reference_interval` to a sequence length compatible with AlphaGenome.
input_interval = variant.reference_interval.resize(
        dna_client.SEQUENCE_LENGTH_1MB
        )
print(f"\nAlphaGenome has predfined sequence lengths. Longest is: {input_interval.width}.\nWe can easily resize our intervals to this length.")


# we can also check if a variant's reference or alt alleles overlap an genome.Interval:

variant = genome.Variant(
        chromosome='chr3',
        position=10_000,
        reference_bases='T',
        alternate_bases='CGTCAAT',
        )
interval = genome.Interval(chromosome='chr3', start=10_005, end=10_010)

# print('Reference overlaps': variant.reference_overlaps(interval)) # False
# print('Alternative overlaps:', variant.alternate_overlaps(interval)) # True
