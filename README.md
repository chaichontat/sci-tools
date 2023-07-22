# sci-tools

[Optimized single-nucleus transcriptional profiling by combinatorial indexing](https://www.nature.com/articles/s41596-022-00752-0)

# Changes to Primers

- Read 1: 37 cycles (dark cycles \[11-17\] bp inclusive, first base is 1)
- Index 1: 10 cycles
- Index 2: 10 cycles
- Read 2: remaining

### Common

- Indices are chosen based on the Levenshtein distance, which includes deletions and insertions, unlike the Hamming distance (even though deletions are much rarer than substitutions in NovaSeq runs).
  - 4 mismatches minimum. 30-70% GC. Indices do not start with `GG`.


### RT

- Switch to anchored dT28VN. [Reduces quantification error from internal priming](https://academic.oup.com/nargab/article/4/2/lqac035/6592171). dT28 to reduce oligo length and to increase proportion of full-length oligos.
- Swap the position of the UMI and index. Makes demultiplexing easier as all indices are contiguous.
- Lengthen the UMI to 10 bp.
- Change the sticky handle to `CTCACTG`.

### Ligation

- Change the sticky handle to `CAGTGAG`.
- The indices are now 10 bp (no more staggering 9-10 bp). Use dark cycles to skip through the ligation handle.
- Full TruSeq Read 1 sequence is included (previous version has a partial sequence, which is extended using PCR).
  - Phosphoramidite synthesis proceeds from 3' to 5'. This means that truncated oligos will not have the sticky handle, effectively eliminating those oligos from the reaction. On the other hand, truncated oligos can still proceed in a PCR reaction.

### PCR
- Full dual-indexing, same barcode selection criteria.
- Barcoded P7 small RNA primers are included for TotalSeq-A-derived ADTs.
