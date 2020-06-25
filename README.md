# graphappy

## A variation graph-oriented fork of WhatsHap.

This fork for WhatsHap explores the application of the minimum error correction approach to sequences embedded in graphs.
Its goal is to enable the phasing of complex loci that arise during genome assembly, which cannot be assembled by alignment, and which instead require phasing and haplotype inference to resolve.

## Usage

To partition a path set in a graph, simply provide a sorted GFA file to `whatshap_graph.py`:

```
python whatshap_graph.py graph.gfa >graph.partition.tsv
```

The results can be used to subset the sequences embedded in the graph for reassembly.

## About WhatsHap

WhatsHap is a software for phasing genomic variants using DNA sequencing
reads, also called *read-based phasing* or *haplotype assembly*. It is
especially suitable for long reads, but works also well with short reads.

For documentation and information on how to cite WhatsHap, please visit the [WhatsHap Homepage](https://whatshap.readthedocs.io/)
