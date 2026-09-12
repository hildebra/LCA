# LCA

[![Anaconda-Server Badge](https://anaconda.org/bioconda/LCA/badges/downloads.svg)](https://anaconda.org/bioconda/LCA)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/LCA/badges/latest_release_relative_date.svg)](https://anaconda.org/bioconda/LCA)


LCA is a command-line program for assigning taxonomy to sequence reads, operational taxonomic units (OTUs), and amplicon sequence variants (ASVs) from their alignments to a reference database. It combines the reference hits using a least common ancestor approach with a configurable consensus fraction, reporting the deepest supported taxonomy for each query. LCA is used in the LotuS1 - LotuS3 pipeline.

The program takes two inputs: a table of sequence alignments and a table linking reference sequence identifiers to taxonomic lineages. Run your sequence search or alignment tool before running LCA. Outputs include per-query assignments, optional identity and reference-hit columns, and counts at each taxonomic rank.

## Installation

Build requirements:

- A C++20 compiler, such as a compatible GCC or Clang installation.
- GNU Make.
- zlib headers and library for the default build.
- Python 3 to run the end-to-end tests; it is not needed to run LCA itself.

From the source directory, build the executable and check its version:

```sh
make
./LCA -v
```

## Quick example

The following shell commands create a small, synthetic reference taxonomy and four alignments for three queries. The taxonomy labels are illustrative. Run the commands from the source directory after building LCA; `printf` writes the required literal tabs.

```sh
mkdir -p example

printf '%b\n' \
  'refA\tk__Bacteria;p__P;c__C;o__O;f__F;g__Genus;s__Genus alpha' \
  'refB\tk__Bacteria;p__P;c__C;o__O;f__F;g__Genus;s__Genus beta' \
  > example/reference.tax

printf '%b\n' \
  'qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tqlen' \
  'query1\trefA\t99\t100\t1\t0\t1\t100\t1\t100\t100' \
  'query2\trefA\t99\t100\t1\t0\t1\t100\t1\t100\t100' \
  'query2\trefB\t99\t100\t1\t0\t1\t100\t1\t100\t100' \
  'query3\trefA\t96\t100\t4\t0\t1\t100\t1\t100\t100' \
  > example/hits.m8

./LCA -i example/hits.m8 -r example/reference.tax \
  -o example/assignments.tsv -showHitRead -reportID

cat example/assignments.tsv
```

Expected contents of `example/assignments.tsv`:

```text
OTU	Domain	Phylum	Class	Order	Family	Genus	Species	Hit2DB	%ID
query1	Bacteria	P	C	O	F	Genus	Genus alpha	refA	99.000000
query2	Bacteria	P	C	O	F	Genus	?	?	99.000000
query3	Bacteria	P	C	O	F	Genus	?	refA	96.000000
```

- `query1` has one eligible reference hit and receives its species assignment.
- `query2` has equally strong hits to two species in the same genus. Neither species meets the default 90% consensus requirement, so the assignment stops at genus. `Hit2DB` is `?` because multiple reference subjects contributed.
- `query3` has one hit at 96% identity. It meets the default genus cutoff of 95%, but not the species cutoff of 97%.

## Input files

### Alignment table (`-i`)

The preferred input is a custom m8-style table with **exactly 11 tab-separated columns**, in this order:

```text
qseqid  sseqid  pident  length  mismatch  gapopen  qstart  qend  sstart  send  qlen
```

The column list above is spaced for readability; actual files must use tabs.

| Field | Meaning |
| --- | --- |
| `qseqid` | Query sequence identifier. |
| `sseqid` | Reference sequence identifier, matching the taxonomy file exactly. |
| `pident` | Percentage identity, from 0 to 100. |
| `length` | Positive integer alignment length. |
| `mismatch`, `gapopen` | Nonnegative integer mismatch and gap-opening counts. |
| `qstart`, `qend` | Positive, one-based inclusive query coordinates, within `qlen`. Reverse alignments are accepted. |
| `sstart`, `send` | Positive, one-based inclusive reference coordinates. |
| `qlen` | Positive integer query sequence length. |

Query coverage is calculated as `min(abs(qend - qstart) + 1, length) / qlen`.

For compatibility, LCA also accepts the standard **12-column BLAST layout**:

```text
qseqid  sseqid  pident  length  mismatch  gapopen  qstart  qend  sstart  send  evalue  bitscore
```

Legacy rows have no `qlen`, so query-coverage filtering is skipped, even when `-cover` is supplied. E-values and bit scores are parsed and validated, but are not used to rank hits. Each file must use one layout consistently; the layout is detected automatically.

Additional requirements:

- All alignments for a query must occur in one contiguous block. Query blocks do not need to be alphabetically sorted, and hits within a block do not need to be ranked. A query that reappears in a later block causes an error.
- A matching tab-separated column header is optional and must precede the data. Blank lines and comment lines beginning with `#` are ignored.
- Spaces in query and subject identifiers are preserved, including Illumina metadata such as `LH00409:413:22WTJGLT4:6:2481:9713:19847 1:N:0:CGTATCTC+CTCGAACA`. Spaces are not field separators.
- Numeric fields and column-header labels may have surrounding whitespace. CRLF line endings are accepted. Empty required fields, invalid numbers, and extra columns are rejected with a line number.
- Mapping files may be plain text or gzip-compressed with a `.gz` extension when gzip support is enabled.

### Reference taxonomy (`-r`)

Supply an **uncompressed text file** with one reference identifier and one semicolon-separated lineage per line, separated by a tab:

```text
reference_identifier<TAB>k__Bacteria;p__P;c__C;o__O;f__F;g__Genus;s__Genus alpha
```

Replace `<TAB>` with a literal tab. Use no column header. Empty lines and lines starting with `#` are ignored. Reference identifiers, including spaces, must match `sseqid` exactly. Every subject retained after alignment filtering must have a taxonomy entry; a missing entry causes an error.

By default, LCA reads and reports seven ranks:

| Rank | Recognized prefix | Default minimum identity |
| --- | --- | --- |
| Domain | `d__` or `k__` | 0% |
| Phylum | `p__` | 78% |
| Class | `c__` | 88% |
| Order | `o__` | 91% |
| Family | `f__` | 93% |
| Genus | `g__` | 95% |
| Species | `s__` | 97% |

Recognized prefixes assign entries to their explicit ranks, even when ranks are omitted or listed out of order. Unprefixed lineages are also accepted in rank order, starting at Domain; preserve empty positions for missing ranks. Empty names, `?`, `unclassified`, `uncultured`, and `uncultured bacterium` are treated as unknown by default.

Use `-tdep 8` to include Strain (`t__`). Additional ranks default to a 97% identity cutoff; custom cutoffs can be supplied with `-id`. Taxonomy beyond the configured depth is ignored.

Species names undergo normalization: ordinary names are shortened to genus and species, with an extra word retained for `Candidatus` names. Names containing a standalone `sp.` or starting with `uncultured` or `unclassified` are treated as uncertain; their species votes can be suppressed when a retained hit supplies a known, certain species.

## Usage and options

```sh
./LCA -i hits.m8 -r reference.tax -o assignments.tsv [options]
./LCA -h
./LCA -v
```

`-i`, `-r`, and `-o` are required for an analysis. `-h` prints help and `-v` prints the version. Boolean flags take no value.

### Filtering and classification

In normal LCA mode, the program filters alignments by minimum length and query coverage, then applies relative length and identity filters around a selected anchor hit. The relative length cutoff is 85% of the anchor's alignment length; the identity window ranges from 0.05 to 1.5 percentage points depending on anchor identity. The anchor normally favors the highest-identity hit, but may favor an alignment at least 20% longer whose identity is at least 90% of that highest identity.

Duplicate alignments to the same subject are collapsed after filtering, giving each retained subject one vote. Each hit's identity limits its eligible taxonomy depth. Consensus proceeds from Domain toward deeper ranks, requiring the configured fraction of known votes at each rank and retaining only descendants of the accepted parent. The mean identity of retained hits also limits the final reported depth.

| Option | Default | Description |
| --- | --- | --- |
| `-minAlignLen N` | `75` | Minimum alignment length; a nonnegative integer. |
| `-cover F` | `0.5` | Minimum query coverage, from 0 to 1. Applies only to input with `qlen`. |
| `-LCAfrac F` | `0.9` | Required fraction of matching known taxonomies at each rank; greater than 0 and at most 1. Set to `1` to require all known votes to agree. |
| `-id LIST` | `97,95,93,91,88,78,0` | Comma-separated identity cutoffs, from the deepest configured rank back to Domain. Supply exactly `-tdep` values, from 0 to 100, in nonincreasing order. |
| `-tdep N` | `7` | Number of taxonomy ranks to read and report, from 1 to 64. Rank 8 is Strain; later ranks are named `Level9`, `Level10`, etc. |
| `-no_bl_filter` | Off | Skip minimum length, coverage, and relative alignment filters when input has already been filtered. Subject deduplication and taxonomy identity cutoffs still apply. |
| `-reportBestHit` | Off | Select one eligible hit instead of computing consensus across subjects. See the cutoff behavior below. |
| `-no_taxDB_filter` | Off | Preserve taxonomy strings otherwise recognized as unknown. Rank parsing and species-name normalization still apply. |

**Best-hit mode changes identity cutoffs:** `-reportBestHit` replaces all rank cutoffs with **1%**, overriding any supplied `-id` values. It still applies minimum length and coverage filters unless `-no_bl_filter` is also set. Hits are ranked by higher identity, longer alignment, higher query coverage when available, and finally lexicographically smaller subject identifier. This mode reports the selected hit's eligible lineage and does not use the normal relative alignment filters.

For example, to require 99% identity at species rank while retaining the other default cutoffs in normal LCA mode:

```sh
./LCA -i hits.m8 -r reference.tax -o assignments.tsv \
  -id 99,95,93,91,88,78,0
```

### Reporting and compatibility

| Option | Description |
| --- | --- |
| `-reportID` | Append `%ID`, the arithmetic mean identity of retained subject hits. In best-hit mode, this is the selected hit's identity. |
| `-showHitRead` | Append `Hit2DB`, the reference identifier when one subject contributes, or `?` when multiple subjects contribute. |
| `-reportHitPattern FILE` | Write query identifier, mean identity, and eligible taxonomy depth to a separate tab-separated file. |
| `-matHigh` | Write one lineage-count file per configured taxonomy rank. |
| `-readInput` | Label the first assignment column `Reads` instead of the default `OTU`; classification and counting are unchanged. |
| `-t 1` | Compatibility option; only one thread is supported. Other values are rejected. |
| `-SLVfmt` | Compatibility flag; recognized taxonomy prefixes are automatically detected without it. |
| `-f bl8` | Compatibility option; `bl8` is the only accepted value and covers both supported mapping layouts. |

### Multiple reference databases

Pass comma-separated mapping files and taxonomy files in corresponding order:

```sh
./LCA -i hits_db1.m8,hits_db2.m8 \
  -r taxonomy_db1.tax,taxonomy_db2.tax \
  -o combined.tsv -reportID
```

The -i and -r lists must have the same length. LCA computes assignments separately for each database and selects one result per query. Selection prioritizes more known eligible ranks, then greater eligible depth, then higher mean identity, with deterministic tie-breakers. Hits from different databases are not pooled into a single consensus. Combined output is sorted by query identifier. File paths containing commas cannot be represented in these lists.

## Output files

The main output (`-o`) is a tab-separated table with a header: `OTU` (or `Reads`), the configured taxonomy ranks, and any requested `Hit2DB` and `%ID` columns. Unknown or unsupported ranks are written as `?`. Queries absent from the mapping input or with no alignments surviving filtering have no output row. A surviving query can have all ranks reported as unknown.

With `-reportHitPattern FILE`, the header is `OTU/ASV`, `ID`, `TaxDepth`. Here `ID` means percentage identity, and `TaxDepth` is the eligible depth, with Domain at depth 1 and Species at depth 7; it is not necessarily a count of nonempty ranks.

With `-matHigh`, files are named by appending the rank to the full output path, for example `assignments.tsv_Domain` and `assignments.tsv_Species`. Each headerless file contains a semicolon-separated lineage followed by a tab and its query count, sorted by lineage. Each reported query contributes one count; OTU abundance encoded in an identifier is not interpreted. Combining `-matHigh` with `-showHitRead` also creates `assignments.tsv_Hit2DB`.

To generate these additional reports for the quick example:

```sh
./LCA -i example/hits.m8 -r example/reference.tax \
  -o example/reported.tsv -reportID -showHitRead -matHigh \
  -reportHitPattern example/hit-pattern.tsv
```

Output paths must be distinct from all inputs and from one another, including generated count files and existing symbolic or hard-link aliases. Existing output files are overwritten. Input read errors, corrupt gzip streams, and output write failures return a nonzero exit status; a failed run may leave partial output files.

## Tests

Run the C++ regression suite and the Python end-to-end checks with:

```sh
make check
```
## Citation

If you use LCA, please cite:

Hildebrand F, Moitinho-Silva L, Blasche S, et al. *Antibiotics-induced monodominance of a novel gut bacterial order.* Gut (2019). DOI: [10.1136/gutjnl-2018-317715](https://doi.org/10.1136/gutjnl-2018-317715).

## License

LCA is distributed under the GNU General Public License, version 3. See [LICENSE](LICENSE).
