#!/usr/bin/env python
from whatshap.core import ReadSet, PedigreeDPTable, Pedigree, NumericSampleIds, PhredGenotypeLikelihoods
from whatshap.testhelpers import string_to_readset, canonic_index_list_to_biallelic_gt_list

import textwrap
from collections import defaultdict
from whatshap.core import Read, ReadSet, readselection, Genotype
from whatshap.merge import ReadMerger

import sys
import math

def gfa_to_readset(gfa_filename, split_gap=100, w=None, sample_ids=None, source_id=0, scale_quality=None):
    rs = ReadSet()
    node_length = {}
    node_coverage = {}
    with open(gfa_filename) as gfa_file:
        for line in gfa_file:
            fields = line.strip().split("\t")
            if fields[0] != "S":
                continue
            node_length[int(fields[1])] = len(fields[2])
    with open(gfa_filename) as gfa_file:
        for line in gfa_file:
            fields = line.strip().split("\t")
            if fields[0] != "P":
                continue
            path_name = fields[1]
            path_str = fields[2]
            for i in [int(s[:-1]) for s in path_str.split(",")]:
                if i in node_coverage:
                    node_coverage[i] += 1
                else:
                    node_coverage[i] = 1
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
            longest_read = None
            while i < path_length:
                read = Read("{}\t{}".format(path_name, segment_idx), 50, source_id)
                #read = Read("{}".format(path_name), 50, source_id)
                segment_idx += 1
                q = 1
                # while the distance to the next node is less than our split_gap threshold
                curr = path[i]
                read.add_variant(position=curr, allele=1, quality=-10*math.log10(1-1.0/node_coverage[curr]+0.001))
                last = curr
                i += 1
                while i < path_length:
                    curr = path[i]
                    dist = 0
                    for node_id in range(last+1, curr):
                        dist += node_length[node_id]
                    #eprint("for", path_name, "dist is", dist)
                    if dist > split_gap:
                        break
                    else:
                        for node_id in range(last+1, curr):
                            #eprint(node_coverage[node_id])
                            read.add_variant(position=node_id, allele=0, quality=1) #-10*math.log10(1-1.0/node_coverage[node_id]+0.001))
                        read.add_variant(position=curr, allele=1, quality=-10*math.log10(1-1.0/node_coverage[curr]+0.001))
                        i += 1
                        last = curr
                #read.sort()  # not sure if needed
                #if len(read) > min_read_length:
                if longest_read is None or len(read) > len(longest_read):
                    longest_read = read
                #rs.add(read)
            rs.add(longest_read)
    rs.sort()
    #print(rs)
    return rs

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

#print('INPUT READ SET')
readset = gfa_to_readset(sys.argv[1], int(sys.argv[2]))
readset = readset.subset(
    [i for i, read in enumerate(readset) if len(read) >= 2]
)

max_coverage = 17
selected_indices = readselection(readset, max_coverage)
selected_reads = readset.subset(selected_indices)

readset_length = 0
for read in selected_reads: readset_length+=len(read)

#print(selected_reads)

def bipartition(reads):
    positions = reads.get_positions()
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
    dp_table = PedigreeDPTable(reads, recombcost, pedigree, distrust_genotypes=False)
    phasing, transmission_vector = dp_table.get_super_reads()
    print('PHASING')
    print(phasing[0])
    mec_score = dp_table.get_optimal_cost()
    eprint("MEC Score:", mec_score)
    eprint("MEC Score / readset length:", float(mec_score) / float(readset_length))

    # In case the bi-partition of reads is of interest:
    partition = dp_table.get_optimal_partitioning()
    #print(partition)
    eprint("partition fraction:", sum(partition)/float(len(partition)))

    return partition

partition = bipartition(selected_reads)

read_partitions = {}

readsets = [ReadSet(), ReadSet()]

for read_index, read in enumerate(selected_reads):
    read_partitions[read.name] = [partition[read_index]]
    readsets[partition[read_index]].add(read)
    #print(partition[read_index], read.name)

partitions = [bipartition(readsets[0]), bipartition(readsets[1])]

for i in range(0,2):
    for read_index, read in enumerate(readsets[i]):
        read_partitions[read.name].append(partitions[i][read_index])

for read_index, read in enumerate(selected_reads):
    print(" ".join(map(str, read_partitions[read.name])), read.name)
#print(len(partition), len(selected_reads))
#print(selected_reads.subset(partition))
#print([i.name for i in selected_reads.subset(partition == 0)])

