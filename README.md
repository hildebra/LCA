# LCA
least common ancestor algorithm that is used in lOTUs pipelin (http://psbweb05.psb.ugent.be/lotus/).
If used, please cite 1 Hildebrand F, Moitinho-Silva L, Blasche S, et al. Antibiotics-induced monodominance of a novel gut bacterial order. Gut 2019;:gutjnl-2018-317715. doi:10.1136/gutjnl-2018-317715

The preferred mapping input is a custom, query-grouped m8-style whitespace-delimited file with exactly these columns (tabs are conventional):

`qseqid`, `sseqid`, `pident`, `length`, `mismatch`, `gapopen`, `qstart`, `qend`, `sstart`, `send`, `qlen`.

A matching header row is optional. Comment lines beginning with `#` and blank lines are ignored; malformed rows fail with their line number.

For legacy compatibility, the standard 12-column BLAST layout is also accepted:

`qseqid`, `sseqid`, `pident`, `length`, `mismatch`, `gapopen`, `qstart`, `qend`, `sstart`, `send`, `evalue`, `bitscore`.

Legacy rows do not contain `qlen`, so query-coverage filtering is skipped for them. A file must consistently use one of the two layouts.

Build with a C++20 compiler and zlib using `make`. The current implementation is single-threaded.
Run the focused correctness checks with `make check`.

Use `-tdep N` to read and report a non-default number of taxonomy levels (default: 7). Recognized `d__`/`k__`, `p__`, `c__`, `o__`, `f__`, `g__`, `s__`, and `t__` prefixes are assigned to their explicit ranks; unprefixed entries retain sequential behavior.
