# LCA

[![Anaconda-Server Badge](https://anaconda.org/bioconda/LCA/badges/downloads.svg)](https://anaconda.org/bioconda/LCA)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/LCA/badges/latest_release_relative_date.svg)](https://anaconda.org/bioconda/LCA)


LCA is a command-line program for assigning taxonomy to sequence reads, operational taxonomic units (OTUs), and amplicon sequence variants (ASVs) from their alignments to a reference database. It combines the reference hits using a least common ancestor approach with a configurable consensus fraction, reporting the deepest supported taxonomy for each query. LCA is used in the LotuS1 - LotuS3 pipeline.

The program takes two inputs: a table of sequence alignments and a table linking reference sequence identifiers to taxonomic lineages. Run your sequence search or alignment tool before running LCA. Outputs include per-query assignments, optional identity and reference-hit columns, and counts at each taxonomic rank.

## Installation

Released versions are available from Bioconda, which can lag behind this repository:

```sh
conda install -c bioconda lca
```

To build from source, you need:

- A C++20 compiler; tested with GCC 13.3 and Clang 18.1 on Linux.
- GNU Make.
- zlib headers and library for the default build.
- Perl 5 with `prove` to run the end-to-end tests; it is not needed to run LCA itself.

From the source directory, build the executable and check its version:

```sh
make
./LCA -v
```

To build without zlib, and therefore without gzip input support:

```sh
make CPPFLAGS=-DLCA_NO_GZIP LDLIBS=
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

- All alignments for a query must occur in one contiguous block. Query blocks do not need to be alphabetically sorted, and hits within a block do not need to be ranked. A query that reappears within the next 100,000 queries causes an error. To bound memory, only the 100,000 most recently completed queries are tracked; a query that reappears later is not detected and receives an additional output row.
- Matching tab-separated column headers are optional and may appear anywhere, for example in concatenated files. Blank lines and comment lines beginning with `#` are ignored.
- Spaces in query and subject identifiers are preserved, including Illumina metadata such as `LH00409:413:22WTJGLT4:6:2481:9713:19847 1:N:0:CGTATCTC+CTCGAACA`. Spaces are not field separators.
- Numeric fields and column-header labels may have surrounding whitespace. CRLF line endings are accepted.
- Malformed records are skipped with a warning: wrong column counts, a row whose layout differs from earlier records, empty required fields, invalid numbers, or out-of-range values such as the zero coordinates of no-hit placeholder rows. The first five skipped records are listed by line number, followed by a count per file. A file with no valid records causes an error.
- Mapping files may be plain text or gzip-compressed with a `.gz` extension when gzip support is enabled.

### Producing the alignment table

These commands, as used by LotuS, write the 11-column table and keep each query's hits together. Replace the file names with your own; `reference` is a database made with `makeblastdb`, and `reference.lba.gz` an index made with `lambda3 mkindexn`:

```sh
blastn -query queries.fasta -db reference -out hits.m8 \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen'
```

```sh
vsearch --usearch_global queries.fasta --db reference.fasta --id 0.75 --maxaccepts 100 \
  --userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+ql --userout hits.m8
```

```sh
lambda3 searchn -q queries.fasta -i reference.lba.gz -o hits.m8 \
  --output-columns 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen'
```

Report enough hits per query for the consensus to see competing references; LotuS uses 100 to 200. BLAST and VSEARCH report only the first word of each FASTA header by default, so the taxonomy file must use those words as identifiers. Plain `blastn -outfmt 6` output is the 12-column layout and works without query-coverage filtering.

### Reference taxonomy (`-r`)

Supply an **uncompressed text file** with one reference identifier and one semicolon-separated lineage per line, separated by a tab:

```text
reference_identifier<TAB>k__Bacteria;p__P;c__C;o__O;f__F;g__Genus;s__Genus alpha
```

Replace `<TAB>` with a literal tab. Use no column header. Empty lines and lines starting with `#` are ignored. Columns after the lineage are ignored with a warning. Lines without an identifier, a tab, and a lineage are skipped with a warning; a file with no valid lines causes an error. Reference identifiers, including spaces, must match `sseqid` exactly. Every subject retained after alignment filtering must have a taxonomy entry; a missing entry causes an error.

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

A hit is eligible for a rank when its identity is at least the rank's cutoff. Before version 0.28, the identity had to exceed the cutoff, and `-id` values were truncated to integers; hits at exactly a cutoff, such as 97.0%, are therefore assigned one rank deeper than in earlier versions.

Recognized prefixes assign entries to their explicit ranks, even when ranks are omitted or listed out of order. Unprefixed lineages are also accepted in rank order, starting at Domain; preserve empty positions for missing ranks. Empty names, `?`, `unclassified`, `uncultured`, `uncultured bacterium`, `unknown`, `unidentified`, and `unassigned` are treated as unknown by default, regardless of case.

Use `-tdep 8` to include Strain (`t__`). Additional ranks default to a 97% identity cutoff; custom cutoffs can be supplied with `-id`. Taxonomy beyond the configured depth is ignored.

Species names undergo normalization: ordinary names are shortened to genus and species, with an extra word retained for `Candidatus` names. Placeholder names are treated as uncertain and are not shortened: names starting with `uncultured` or `unclassified`, or containing any of the words `sp`, `sp.`, `spp`, `spp.`, `cf`, `cf.`, `aff`, `aff.`, `bacterium`, `archaeon`, `metagenome`, `unidentified`, `unknown`, `environmental`, `uncultured`, `unclassified`, or `incertae`, regardless of case. Words are separated by spaces or underscores, so `Lachnospiraceae bacterium 3_1_46FAA` and `Bacteroides_sp._XY1` are placeholders, while GTDB names such as `Bacteroides_A sp000123456` are not. Uncertain species votes are suppressed when a retained hit supplies a known, certain species.

## Usage and options

```sh
./LCA -i hits.m8 -r reference.tax -o assignments.tsv [options]
./LCA -h
./LCA -v
```

`-i`, `-r`, and `-o` are required for an analysis. `-h` prints help and `-v` prints the version. Boolean flags take no value.

### Filtering and classification

In normal LCA mode, the program filters alignments by minimum length and query coverage, then applies relative coverage and identity filters around an anchor hit. The anchor is the highest-identity hit; query coverage, then alignment length, only break identity ties. Retained hits must cover at least 85% as much of the query as the anchor, and their identity must fall within a window of 0.05 to 1.5 percentage points below the anchor, depending on anchor identity. Legacy rows without `qlen` use 85% of the anchor's alignment length instead of coverage.

Duplicate alignments to the same subject are collapsed after filtering, giving each retained subject one vote. Each hit's identity limits its eligible taxonomy depth. Consensus proceeds from Domain toward deeper ranks, requiring the configured fraction of known votes at each rank. Hits naming a different taxon at an accepted rank stop voting. Hits without a name at a rank abstain there but still vote at deeper ranks, so lineages with a missing rank, such as an absent order, are resolved the same way for one subject or several, and their matching deeper names are used. The mean identity of retained hits also limits the final reported depth.

| Option | Default | Description |
| --- | --- | --- |
| `-minAlignLen N` | `75` | Minimum alignment length; a nonnegative integer. |
| `-cover F` | `0.5` | Minimum query coverage, from 0 to 1. Applies only to input with `qlen`. |
| `-LCAfrac F` | `0.9` | Required fraction of matching known taxonomies at each rank; greater than 0 and at most 1. Set to `1` to require all known votes to agree. |
| `-id LIST` | `97,95,93,91,88,78,0` | Comma-separated minimum identities (inclusive), from the deepest configured rank back to Domain. Supply exactly `-tdep` values, from 0 to 100, in nonincreasing order; decimals are allowed. |
| `-tdep N` | `7` | Number of taxonomy ranks to read and report, from 1 to 64. Rank 8 is Strain; later ranks are named `Level9`, `Level10`, etc. |
| `-no_bl_filter` | Off | Skip minimum length, coverage, and relative coverage and identity filters when input has already been filtered. Subject deduplication and taxonomy identity cutoffs still apply. |
| `-reportBestHit` | Off | Use the best eligible hit instead of computing consensus across all retained subjects. See the tie and cutoff behavior below. |
| `-no_taxDB_filter` | Off | Preserve taxonomy strings otherwise recognized as unknown. Rank parsing and species-name normalization still apply. |

**Best-hit mode changes identity cutoffs:** `-reportBestHit` replaces all rank cutoffs with **1%**, overriding any supplied `-id` values. It still applies minimum length and coverage filters unless `-no_bl_filter` is also set. Hits are ranked by higher identity, then higher query coverage when available, then longer alignment. This mode reports the best hit's eligible lineage and does not use the normal relative filters. When several subjects tie for best, their lineages are combined by the normal consensus with `-LCAfrac`, so tied subjects of different species are reported at genus.

For example, to require 99% identity at species rank while retaining the other default cutoffs in normal LCA mode:

```sh
./LCA -i hits.m8 -r reference.tax -o assignments.tsv \
  -id 99,95,93,91,88,78,0
```

### Reporting and compatibility

| Option | Description |
| --- | --- |
| `-reportID` | Append `%ID`, the arithmetic mean identity of retained subject hits. In best-hit mode, this is the selected hit's identity. |
| `-showHitRead` | Append `Hit2DB`, the reference identifier when one subject contributes, or `?` when multiple subjects contribute. In best-hit mode with tied subjects, it names the top-ranked tied subject whose lineage contains the reported consensus, lexicographically first among equals. |
| `-reportHitPattern FILE` | Write query identifier, mean identity (formatted like `%ID`), and eligible taxonomy depth to a separate tab-separated file. |
| `-matHigh` | Write one lineage-count file per configured taxonomy rank. |
| `-readInput` | Label the first assignment column `Reads` instead of the default `OTU`; classification and counting are unchanged. |
| `-t N` | Threads for loading the reference taxonomy and computing assignments; default 1. Output is identical for any value. Reading the mapping file stays on one thread, so gains level off at about 2-4 threads. Values that are not a positive integer mean one thread. |
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

Output paths must be distinct from all inputs and from one another, including generated count files and existing symbolic or hard-link aliases. Existing output files are overwritten. Input read errors, corrupt gzip streams, and output write failures return a nonzero exit status; a failed run may leave partial output files. Skipped malformed input lines are reported on standard error and do not change the exit status.

## Exit status

LCA exits with 0 on success, including for `-h` and `-v`. Other values identify the problem; the message on standard error gives the details.

| Status | Meaning |
| --- | --- |
| 5, 6, 24, 39 | Invalid command line: a required argument is missing, an option is unknown or has an invalid value, the `-i` and `-r` lists differ in length, or the `-id` cutoffs are invalid. |
| 11, 12, 13, 14 | Taxonomy file: it cannot be opened, has no valid records, cannot be read, or has too many entries. |
| 23, 25, 27, 28, 50 | Mapping file: it cannot be opened, has no valid records, lists a query in non-contiguous blocks, cannot be read or has corrupt gzip data, or is gzip-compressed in a build without gzip support. |
| 29-34 | Output: a path overlaps an input or another output, a file cannot be created, or writing failed. |
| 74 | The subject of a retained hit is missing from the taxonomy file. |

## Performance and memory

- The reference taxonomy needs about 100 bytes of memory per reference, about 100 MB for one million references. With multiple databases, one taxonomy is loaded at a time.
- `-t N` loads the taxonomy and computes assignments in parallel, and the output does not depend on the thread count. Reading the mapping file stays on one thread, so 2 to 4 threads give most of the gain. In a benchmark with one million references, an ASV run with 20,000 queries took 1.3 s with `-t 2` and 2.0 s with one thread.
- Pipelines that call LCA, such as LotuS, need to pass `-t` to use threads.
- Gzip-compressed mapping files are read almost as fast as plain text. Linking against zlib-ng in zlib-compatible mode can speed up decompression further.

## Using LCA as a library

`src/LCAlib.h` gives a program that holds a query's hits in memory the same assignment that LCA computes from files. It is a single C++20 header that needs only the standard library, so it can be copied into another source tree. protal uses it to place rRNA reads.

```cpp
#include "LCAlib.h"

lca::Taxonomy tax;  // seven ranks, Domain to Species
tax.AddLineage("k__Bacteria;p__P;c__C;o__O;f__F;g__Genus;s__Genus alpha");  // subject 0
tax.AddLineage("k__Bacteria;p__P;c__C;o__O;f__F;g__Genus;s__Genus beta");   // subject 1

lca::Options opt;      // LCA's defaults
lca::Scratch scratch;  // one per thread, reused for every query
std::vector<lca::Hit> hits = {{0, 99.0, 100, 100, 100}, {1, 99.0, 100, 100, 100}};
lca::Assignment a = lca::Assign(hits, tax, opt, scratch);
// a.depth == 6; tax.Name(a.name_ids[rank]) gives Bacteria ... Genus for rank 0 to 5
```

- `Taxonomy` stores lineages by subject index, in the order they are added. `AddLineage` reads a lineage as the taxonomy file has it, and `Add` takes names per rank from Domain down, without prefixes. Both apply the naming rules described under [Reference taxonomy](#reference-taxonomy--r), so unknown names become `0` and species names are shortened or marked as placeholders. `NameId(subject, rank)` returns an interned name ID, where `0` means unknown. `Name(id)` returns its text, and `Name(0)` is `?`.
- `Hit` holds `subject` (the index into the taxonomy), `identity` in percent, `align_len`, `query_len` and `query_covered`. A `query_len` of 0 means the query length is unknown, as in 12-column BLAST input; the coverage filter is then skipped. Several hits may share a subject. Hits of subjects outside the taxonomy are ignored.
- `Options` mirrors the command-line options: `min_identity` (`-id`, deepest rank first), `lca_frac` (`-LCAfrac`), `min_align_len` (`-minAlignLen`), `min_cover` (`-cover`), `best_hit` (`-reportBestHit`), `no_filters` (`-no_bl_filter`) and `depth` (`-tdep`). Unlike `-reportBestHit`, `best_hit` keeps the supplied cutoffs; set every `min_identity` to 1 to match the program.
- `Assignment` reports `depth`, the number of ranks assigned from Domain down, ending at the deepest named rank. It also reports `name_ids` per rank, where `0` marks a rank that no retained hit names, and `mean_identity` (as `%ID`). `hits_used` is the number of retained subjects; 0 means no hit passed the filters, where the program writes no row. `species_uncertain` is set when the assigned species is a placeholder name.
- `Parts()` copies out a `Taxonomy`'s names and name IDs so they can be saved, and `Taxonomy::FromParts` rebuilds the `Taxonomy` from them without looking up every subject's names again. Saved parts are valid only under the `kNamingRulesVersion` they were made with; under another version, build the `Taxonomy` from the lineages again.
- A `Taxonomy` is read-only once built, so any number of threads can call `Assign` with it at the same time, each with its own `Scratch`. `Assign` allocates no memory once its `Scratch` has grown to the largest query.

`make check` compares the library with the program on random queries in every filtering mode, so a change to the rules in `src/LCAimpl.cpp` must be made in `src/LCAlib.h` as well.

## Changes in recent versions

Results can differ from earlier versions for these reasons, described in the sections above:

- 0.28: identity cutoffs are inclusive, and `-id` accepts decimals.
- 0.31: the highest-identity hit anchors the filters, and the relative filter uses query coverage instead of alignment length; names below a missing rank are used; tied best hits in `-reportBestHit` mode are combined by consensus; more placeholder species names are recognized; malformed input lines are skipped with warnings instead of stopping the run.
- 0.32: faster input, about 4x less memory, and `-t` threads; assignments are identical to 0.31.
- 0.33: `src/LCAlib.h` provides the assignment as a library, and the sources moved to `src/`; the program's output is identical to 0.32.
- 0.34: `LCAlib.h` saves and rebuilds a Taxonomy (`Taxonomy::Parts`, `Taxonomy::FromParts`, `kNamingRulesVersion`); the program's output is identical to 0.33.

## Tests

`make check` builds LCA and runs the C++ unit tests in `tests/regression.cpp` and the end-to-end tests in `tests/cli_regression.t`, which use `prove`:

```sh
make check
```

## Citation

If you use LCA, please cite:

Hildebrand F, Moitinho-Silva L, Blasche S, et al. *Antibiotics-induced monodominance of a novel gut bacterial order.* Gut (2019). DOI: [10.1136/gutjnl-2018-317715](https://doi.org/10.1136/gutjnl-2018-317715).

## License

LCA is distributed under the GNU General Public License, version 2 (GPL-2.0-only). See [LICENSE](LICENSE).

Report bugs and contribute through the [LCA GitHub repository](https://github.com/hildebra/LCA). For questions, contact Falk Hildebrand at `falk.hildebrand@gmail.com`.

