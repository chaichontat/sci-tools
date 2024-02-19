# sci-tools

[Optimized single-nucleus transcriptional profiling by combinatorial indexing](https://www.nature.com/articles/s41596-022-00752-0)

# Changes to Primers

- Read 1: 34 cycles
- Index 1: 10 cycles
- Index 2: 10 cycles
- Read 2: remaining


### RT
- Switch to anchored dT28VN. [Reduces quantification error from internal priming](https://academic.oup.com/nargab/article/4/2/lqac035/6592171). dT28 to reduce oligo length and to increase proportion of full-length oligos.

### PCR
- Full dual-indexing, same barcode selection criteria.
- Barcoded P7 small RNA primers are included for TotalSeq-A-derived ADTs.
- Indices are chosen based on the Levenshtein distance, which includes deletions and insertions, unlike the Hamming distance (even though deletions are much rarer than substitutions in NovaSeq runs).

## Analysis

1. 
