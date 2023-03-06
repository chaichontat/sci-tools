# sci-tools

[Optimized single-nucleus transcriptional profiling by combinatorial indexing](https://www.nature.com/articles/s41596-022-00752-0)

# Changes to Primers

### Common

- Indices are chosen based on the Levenshtein distance, which includes deletions and insertions, unlike the Hamming distance. This is very helpful for demultiplexing the RT indices as the ligation indices are of variable length.
  - 4 mismatches minimum. 30-70% GC. Indices do not start with `GG`.


### RT

- Switch to anchored dT23. [Reduces quantification error from internal priming](https://academic.oup.com/nargab/article/4/2/lqac035/6592171). dT23 to reduce oligo length and to increase proportion of full-length oligos.
- Swap the position of the UMI and index. Makes demultiplexing easier as all indices are contiguous.
- Lengthen the UMI to 12 bp.
- Change the sticky handle to `CTCACTG`.
- Indices are selected such that >94% of oligos do not have internal hairpins. Since the UMI is on this set, the oligo sequence is not fully deterministic (compared with 0-40% from the previous version).
  - This necessitates the removal of all Gs from the UMI and indices.
  - 5% PhiX spike-in would be necessary on the NovaSeq. Seems like a good trade-off for a more complex library.

### Ligation

- Change the sticky handle to `CAGTGAG`.
- The indices are now 10-11 bp (9-10 bp before) to account for potentially poor base calling toward the start of the read.
- Full TruSeq Read 1 sequence is included (previous version has a partial sequence, which is extended using PCR).
  - Phosphoramidite synthesis proceeds from 3' to 5'. This means that truncated oligos will not have the sticky handle, effectively eliminating those oligos from the reaction. On the other hand, truncated oligos can still proceed in a PCR reaction.

### PCR

- Indices are chosen such that they do not contribute to the weak hairpin between P7 and the Nextera handle.
- Dual index pairs are matched using minimum cost bipartite matching with the cost being heterodimer concentration from NUPACK.
- P7 small RNA primers are included for TotalSeq-A-derived ADTs.
