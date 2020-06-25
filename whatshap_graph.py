#!/usr/bin/env python
from whatshap.core import ReadSet, PedigreeDPTable, Pedigree, NumericSampleIds, PhredGenotypeLikelihoods
from whatshap.testhelpers import string_to_readset, canonic_index_list_to_biallelic_gt_list

import textwrap
from collections import defaultdict
from whatshap.core import Read, ReadSet, readselection, Genotype
from whatshap.merge import ReadMerger

import sys

def gfa_to_readset(gfa_filename, split_gap=100, w=None, sample_ids=None, source_id=0, scale_quality=None):
    rs = ReadSet()
    with open(gfa_filename) as gfa_file:
        for line in gfa_file:
            fields = line.strip().split("\t")
            if fields[0] != "P":
                continue
            path_name = fields[1]
            path_str = fields[2]
            # order it
            path = sorted(set([int(s[:-1]) for s in path_str.split(",")]))
            # break each path into pieces separated by > x nodes (todo: use actual distance in the graph)
            # for each, add it to the ReadSet
            path_length = len(path)
            segment_idx = 0
            i = 0
            # how do we find segments?
            while i < path_length:
                read = Read("{}.{}".format(path_name, segment_idx), 50, source_id)
                segment_idx += 1
                q = 1
                # while the distance to the next node is less than our split_gap threshold
                curr = path[i]
                read.add_variant(position=(curr * 10), allele=1, quality=q)
                last = curr
                i += 1
                while i < path_length:
                    curr = path[i]
                    if curr - last > split_gap:
                        break
                    else:
                        for node_id in range(last+1, curr):
                            read.add_variant(position=(node_id * 10), allele=0, quality=q)
                        read.add_variant(position=(curr * 10), allele=1, quality=q)
                        i += 1
                        last = curr
                #read.sort()  # not sure if needed
                rs.add(read)
    rs.sort()
    print(rs)
    return rs


print('INPUT READ SET')
readset = gfa_to_readset(sys.argv[1])
readset = readset.subset(
    [i for i, read in enumerate(readset) if len(read) >= 2]
)
#read_merging=False,
"""
read_merging_error_rate=0.15,
read_merging_max_error_rate=0.25,
read_merging_positive_threshold=1000000,
read_merging_negative_threshold=1000,
read_merger = ReadMerger(
    read_merging_error_rate,
    read_merging_max_error_rate,
    read_merging_positive_threshold,
    read_merging_negative_threshold,
)
merged_reads = read_merger.merge(readset)
"""
max_coverage = 20
selected_indices = readselection(readset, max_coverage)
selected_reads = readset.subset(selected_indices)

print(readset)

positions = readset.get_positions()

# create genotypes over your variants: all heterozygous (=1)
genotypes = canonic_index_list_to_biallelic_gt_list([1] * len(positions))
# genotype likelihoods are None
genotype_likelihoods = [None] * len(positions)
# create empty pedigree
pedigree = Pedigree(NumericSampleIds())
# add one individual to pedigree
pedigree.add_individual('individual0', genotypes, genotype_likelihoods)
# recombination cost vector, irrelevant if one using one individual
recombcost = [1] * len(positions) 

# run the core phasing algorithm, creating a DP table
dp_table = PedigreeDPTable(readset, recombcost, pedigree, distrust_genotypes=False)
phasing, transmission_vector = dp_table.get_super_reads()
print('PHASING')
print(phasing[0])
print('MEC Score:', dp_table.get_optimal_cost())

# In case the bi-partition of reads is of interest:
partition = dp_table.get_optimal_partitioning()
print(partition)

