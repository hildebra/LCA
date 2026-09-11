# LCA
least common ancestor algorithm that is used in lOTUs pipelin (http://psbweb05.psb.ugent.be/lotus/).
If used, please cite 1 Hildebrand F, Moitinho-Silva L, Blasche S, et al. Antibiotics-induced monodominance of a novel gut bacterial order. Gut 2019;:gutjnl-2018-317715. doi:10.1136/gutjnl-2018-317715

The preferred mapping input is a custom, query-grouped m8-style tab-separated file with exactly these columns:

`qseqid`, `sseqid`, `pident`, `length`, `mismatch`, `gapopen`, `qstart`, `qend`, `sstart`, `send`, `qlen`.

Fields must be separated by tabs. Spaces in query and subject identifiers are preserved, including Illumina header metadata such as `LH00409:413:22WTJGLT4:6:2481:9713:19847 1:N:0:CGTATCTC+CTCGAACA`. Leading and trailing whitespace around numeric values is accepted; it does not create extra columns. Empty required fields are rejected.

A matching tab-separated header row is optional. Comment lines beginning with `#` and blank lines are ignored; malformed rows fail with their line number.

For legacy compatibility, the standard 12-column BLAST layout is also accepted:

`qseqid`, `sseqid`, `pident`, `length`, `mismatch`, `gapopen`, `qstart`, `qend`, `sstart`, `send`, `evalue`, `bitscore`.

Legacy rows do not contain `qlen`, so query-coverage filtering is skipped for them. A file must consistently use one of the two layouts.

Build with a C++20 compiler and zlib using `make`. The current implementation is single-threaded.
Run the C++ and end-to-end correctness checks with `make check` (Python 3 is required for the latter).

Use `-tdep N` to read and report a non-default number of taxonomy levels (default: 7). Recognized `d__`/`k__`, `p__`, `c__`, `o__`, `f__`, `g__`, `s__`, and `t__` prefixes are assigned to their explicit ranks; unprefixed entries retain sequential behavior.

`-reportBestHit` selects one eligible alignment using identity, alignment length, query coverage, and then subject identifier as deterministic tie-breakers. In normal LCA mode, one retained alignment per subject contributes a vote; duplicate alignments are collapsed after filtering so that an ineligible alignment cannot hide a valid one.

Identity cutoffs are applied consistently to the main taxonomy output, hit-pattern depth, abundance matrices, and selection across databases. Output paths must be distinct from all input paths and from each other, including generated matrix files. Input read failures and corrupt gzip streams return a nonzero exit status.

See [AUDIT.md](AUDIT.md) for the September 2026 audit findings and validation scope.
