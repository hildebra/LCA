#include "../src/LCAimpl.h"
#include "../src/LCAlib.h"
#include "../src/Matrix.h"
#include "../src/Pipeline.h"

#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <memory>
#include <new>
#include <numeric>
#include <random>

// Counts the program's heap allocations, so that a check can tell that a call allocated nothing. Not inlined, so
// that GCC does not take the malloc and free inside for a mismatch with new and delete (-Wmismatched-new-delete).
static std::atomic<size_t> heapAllocations{0};
[[gnu::noinline]] void* operator new(std::size_t size) {
	++heapAllocations;
	if (void* memory = std::malloc(size == 0 ? 1 : size)) { return memory; }
	throw std::bad_alloc();
}
[[gnu::noinline]] void operator delete(void* memory) noexcept { std::free(memory); }
[[gnu::noinline]] void operator delete(void* memory, std::size_t) noexcept { std::free(memory); }

namespace {
int failures = 0;
const char* currentGroup = "";

// Reports through stdio so failures stay visible while iostreams are silenced.
void check(bool condition, const char* message) {
	if (!condition) {
		std::fprintf(stderr, "FAIL [%s]: %s\n", currentGroup, message);
		++failures;
	}
}

BlastRes hit(const char* query, const char* subject, double identity,
	int alignmentLength, int queryStop, int queryLength) {
	std::ostringstream line;
	line << query << '\t' << subject << '\t' << identity << '\t' << alignmentLength
		<< "\t0\t0\t1\t" << queryStop << "\t1\t" << alignmentLength << '\t' << queryLength;
	return BlastRes(line.str(), 0);
}

options makeOptions(std::vector<std::string> extra = {}) {
	std::vector<std::string> args = {"regression", "-i", "input.m8", "-r", "reference.tax", "-o", "output.tsv"};
	args.insert(args.end(), extra.begin(), extra.end());
	std::vector<char*> argv;
	for (auto& arg : args) { argv.push_back(arg.data()); }
	return options(static_cast<int>(argv.size()), argv.data(), 7);
}

// Rank names of the TaxObjs and RefTaxes built here; TaxObjs point into it, so it lives as long as the tests.
NamePool testNames;
const std::string* nameOf(std::string_view name) { return &testNames.name(testNames.intern(name)); }

TaxObj* lineageOf(std::initializer_list<const char*> ranks) {
	TaxObj* tax = new TaxObj(7);
	int rank = 0;
	for (const char* name : ranks) {
		if (*name != '\0') { tax->set(rank, nameOf(name)); }
		++rank;
	}
	return tax;
}

void deleteAll(std::vector<TaxObj*>& taxa) {
	for (TaxObj* tax : taxa) { delete tax; }
	taxa.clear();
}

// Captures cout and cerr, e.g. the taxonomy loading report and input warnings.
struct Silence {
	std::ostringstream sink;
	std::streambuf* out;
	std::streambuf* err;
	Silence() : sink(), out(std::cout.rdbuf(sink.rdbuf())), err(std::cerr.rdbuf(sink.rdbuf())) {}
	~Silence() { std::cout.rdbuf(out); std::cerr.rdbuf(err); }
};

struct TempFile {
	std::string path;
	explicit TempFile(const std::string& content, const std::string& suffix = ".txt") {
		static int counter = 0;
		const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
		path = (std::filesystem::temp_directory_path() /
			("lca-regression-" + std::to_string(stamp) + "-" + std::to_string(counter++) + suffix)).string();
		std::ofstream out(path, std::ios::binary);
		out << content;
	}
	~TempFile() { std::error_code ignored; std::filesystem::remove(path, ignored); }
};

void taxonomyParsing() {
	currentGroup = "taxonomy parsing";
	TaxObj prefixed("s__Genus sp.;k__Bacteria;g__Genus", 7, false, false, testNames);
	check(prefixed.get(0) == "Bacteria", "k__ must map to Domain even after s__");
	check(prefixed.get(5) == "Genus", "g__ must map to Genus even after s__");
	check(prefixed.get(6) == "Genus sp.", "s__ must map to Species");
	check(prefixed.speciesUncertain, "a species containing sp. must be marked uncertain");

	TaxObj unprefixed("Bacteria;P;;O;F;G;Genus species", 7, false, false, testNames);
	check(unprefixed.get(2) == __unkwnTax && unprefixed.get(3) == "O" && unprefixed.depth == 7,
		"empty positions in unprefixed lineages must keep later ranks in place");
	TaxObj padded(" k__ Bacteria ; p__ P ;", 7, false, false, testNames);
	check(padded.get(0) == "Bacteria" && padded.get(1) == "P" && padded.depth == 2,
		"whitespace around names and a trailing separator must be ignored");
	TaxObj upper("D__Bacteria;P__Firmicutes", 7, false, false, testNames);
	check(upper.get(0) == "Bacteria" && upper.get(1) == "Firmicutes", "rank prefixes must be case-insensitive");
	TaxObj kept("k__Bacteria;g__uncultured", 7, false, true, testNames);
	TaxObj filtered("k__Bacteria;g__uncultured", 7, false, false, testNames);
	check(kept.get(5) == "uncultured" && filtered.get(5) == __unkwnTax && filtered.depth == 1,
		"-no_taxDB_filter must keep names that are otherwise treated as unknown");
	TaxObj unknownGenus("k__Bacteria;p__P;g__Unknown;s__Unidentified", 7, false, false, testNames);
	check(unknownGenus.get(5) == __unkwnTax && unknownGenus.get(6) == __unkwnTax && unknownGenus.depth == 2,
		"unknown and unidentified rank names must be treated as missing");
	TaxObj strain8("k__Bacteria;s__Genus species;t__strain X", 8, false, false, testNames);
	TaxObj strain7("k__Bacteria;s__Genus species;t__strain X", 7, false, false, testNames);
	check(strain8.get(7) == "strain X" && strain7.depth == 7, "t__ must be read only when Strain is configured");
}

void speciesNames() {
	currentGroup = "species names";
	TaxObj candidatus("k__Bacteria;p__P;c__C;o__O;f__F;g__Liberibacter;s__Candidatus Liberibacter asiaticus strain X",
		7, false, false, testNames);
	check(candidatus.get(6) == "Candidatus Liberibacter asiaticus",
		"Candidatus species names must retain their three-part species name");
	TaxObj serotype("k__Bacteria;s__Escherichia coli O157:H7", 7, false, false, testNames);
	TaxObj subspecies("k__Bacteria;s__Lactobacillus casei subsp. casei", 7, false, false, testNames);
	check(serotype.get(6) == "Escherichia coli" && subspecies.get(6) == "Lactobacillus casei",
		"certain species names must be shortened to genus and species");
	const char* placeholders[] = {"Lachnospiraceae bacterium 3_1_46FAA", "Methanomassiliicoccales archaeon RumEn M1",
		"human gut metagenome", "Blautia sp", "Bacteroides_sp._XY1", "Lactobacillus cf. casei",
		"Streptococcus aff. mitis", "unidentified Prevotella", "Ruminococcus spp."};
	for (const char* name : placeholders) {
		TaxObj tax(std::string("k__Bacteria;g__G;s__") + name, 7, false, false, testNames);
		check(tax.speciesUncertain && tax.get(6) == name, "placeholder species names must be marked uncertain, unshortened");
	}
	const char* species[] = {"Mycobacterium tuberculosis", "Bacteroides_A sp000123456", "Blautia obeum",
		"Spirochaeta spectabilis"};
	for (const char* name : species) {
		TaxObj tax(std::string("k__Bacteria;g__G;s__") + name, 7, false, false, testNames);
		check(!tax.speciesUncertain, "real and GTDB placeholder-accession species must stay certain");
	}
}

void mappingRecords() {
	currentGroup = "mapping records";
	BlastRes crlf("q\ts\t99\t100\t0\t0\t100\t1\t200\t101\t100\r", 0);
	check(!crlf.fail && crlf.Query == "q" && crlf.Sbj == "s", "CRLF m8 records must parse");
	check(std::fabs(crlf.Qcoverage - 1.0) < 1e-12, "reverse query coordinates must produce valid coverage");
	check(BlastRes::isColumnHeader(
		"qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tqlen\r\n"),
		"the documented custom m8 header must be recognized");
	check(BlastRes::isColumnHeader(
		"qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"),
		"the legacy 12-column BLAST header must be recognized");
	check(!BlastRes::isColumnHeader("qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen"),
		"headers must use tabs as column separators");
	check(BlastRes("q s 99 100 0 0 1 100 1 100 100", 0).fail,
		"ordinary spaces must not be interpreted as column separators");

	BlastRes legacy(
		"HSQ-700358:311:H23FFBCXX:1:2210:10499:32978\tAACY020299484\t96.70\t91\t3\t0\t203\t293\t1912\t2002\t1e-35\t 152",
		0);
	check(!legacy.fail && !legacy.queryCoverageKnown && legacy.eval == 1e-35 && legacy.score == 152.0,
		"legacy qseqid-to-bitscore 12-column BLAST records must parse");
	check(BlastRes("q\ts\t99\t100\t0\t0\t1\t100\t1\t100\t-1\t50", 0).fail, "negative e-values must be rejected");

	BlastRes gapped("q\ts\t99\t100\t0\t1\t1\t90\t1\t100\t100", 0);
	check(!gapped.fail && std::fabs(gapped.Qcoverage - 0.9) < 1e-12 && gapped.queryCovered == 90,
		"coverage must use the query-coordinate span when it is shorter than the alignment");
	BlastRes overhang = hit("q", "s", 99, 100, 120, 200);
	check(!overhang.fail && overhang.queryCovered == 100 && overhang.Qcoverage == 0.5,
		"coverage must use the alignment length when it is shorter than the query span");

	BlastRes padded("  q description \t s description \t 99 \t 100 \t0\t0\t1\t100\t1\t100\t 100 \r\n", 0);
	check(!padded.fail && padded.Query == "  q description " && padded.Sbj == " s description ",
		"numeric padding and CRLF must be accepted without trimming identifiers");
	const std::string illumina = "LH00409:413:22WTJGLT4:6:2481:9713:19847 1:N:0:CGTATCTC+CTCGAACA";
	const std::string reportedLine = illumina + "\tOBEP010038528\t 95.3 \t150\t7\t0\t1\t150\t795\t944\t150";
	BlastRes reported(reportedLine, 0);
	check(!reported.fail && reported.Query == illumina && reported.Sbj == "OBEP010038528" &&
		reported.perID == 95.3 && reported.queryCoverageKnown && reported.Qcoverage == 1.0,
		"the reported Usearch row must preserve its complete Illumina query header");
	check(BlastRes::supportedColumnCount(reportedLine) == 11,
		"an embedded identifier space must not change the detected layout");
	std::string queryToken = "stale";
	check(!BlastRes::extractQueryToken(" \t\r\n", queryToken), "blank lines must not contain a query token");
	check(BlastRes::extractQueryToken("  query with barcode \tsubject", queryToken) &&
		queryToken == "  query with barcode ", "query extraction must preserve the complete first tab-separated field");
	check(BlastRes::extractQueryToken(reportedLine, queryToken) && queryToken == illumina,
		"query extraction must agree with parsing for Illumina metadata");

	std::vector<std::string> requiredFields = {"q", "s", "99", "100", "0", "0", "1", "100", "1", "100", "100"};
	for (size_t missing = 0; missing < requiredFields.size(); ++missing) {
		std::string line;
		for (size_t column = 0; column < requiredFields.size(); ++column) {
			if (column != 0) { line += '\t'; }
			if (column != missing) { line += requiredFields[column]; }
		}
		check(BlastRes::supportedColumnCount(line) == 11 && BlastRes(line, 0).fail,
			"empty tab-separated fields must be preserved and rejected, never collapsed");
	}
	const char* invalid[] = {
		"q\ts\tbad\t100\t0\t0\t1\t90\t1\t100\t100",            // non-numeric identity
		"q\ts\t9 9\t100\t0\t0\t1\t100\t1\t100\t100",           // space inside a number
		"q\ts\tnan\t100\t0\t0\t1\t100\t1\t100\t100",           // not finite
		"q\ts\tinf\t100\t0\t0\t1\t100\t1\t100\t100",
		"q\ts\t100.5\t100\t0\t0\t1\t100\t1\t100\t100",         // identity above 100
		"q\ts\t-1\t100\t0\t0\t1\t100\t1\t100\t100",
		"q\ts\t99\t99999999999\t0\t0\t1\t100\t1\t100\t100",    // integer overflow
		"q\ts\t99\t100.5\t0\t0\t1\t100\t1\t100\t100",          // fractional integer
		"q\ts\t99\t100\t0\t0\t0\t99\t1\t100\t100",             // zero coordinate
		"q\ts\t99\t100\t0\t0\t1\t101\t1\t100\t100",            // qend beyond qlen
		"q\ts\t99\t100\t0\t0\t1\t100\t1\t100\t100\t1e-20\t50", // 13 columns
		"q\ts\t99\t100\t0\t0\t1\t100\t1\t100\t100\t1\t2\t3",
	};
	for (const char* line : invalid) { check(BlastRes(line, 0).fail, line); }
	BlastRes plus("q\ts\t+99\t100\t0\t0\t1\t100\t1\t100\t100", 0);
	BlastRes exponent("q\ts\t9.9e1\t100\t0\t0\t1\t100\t1\t100\t100", 0);
	check(!plus.fail && plus.perID == 99.0 && !exponent.fail && exponent.perID == 99.0,
		"a leading plus sign and exponent notation must be accepted as before");
}

void hitRanking() {
	currentGroup = "hit ranking";
	BlastRes a = hit("q", "a", 99, 100, 100, 100);
	BlastRes b = hit("q", "b", 99, 100, 100, 100);
	check(a.compareRank(b) == 0 && a.betterThan(b) && !b.betterThan(a),
		"equally ranked hits must be ordered only by subject identifier");
	BlastRes higherIdentity = hit("q", "z", 99.1, 100, 60, 100);
	BlastRes fullCoverage = hit("q", "a", 99.0, 100, 100, 100);
	check(higherIdentity.betterThan(fullCoverage), "identity must outrank query coverage");
	std::vector<BlastRes> sameIdentity = {hit("q", "longer_gapped", 99.0, 140, 100, 120),
		hit("q", "covering", 99.0, 115, 115, 120)};
	check(sameIdentity[1].betterThan(sameIdentity[0]),
		"at equal identity, query coverage must outrank alignment length");
}

void primaryFilter() {
	currentGroup = "primary filter";
	options opt = makeOptions();
	double bestID = 0.0;
	opt.minCover = 1.0;
	std::vector<BlastRes> legacyHits = {BlastRes(
		"q\ts\t96.70\t91\t3\t0\t203\t293\t1912\t2002\t1e-35\t152", 0)};
	filterBlastPrimary(legacyHits, &opt, bestID);
	check(legacyHits.size() == 1, "coverage filtering must not reject legacy records whose query length is unavailable");
	opt.minCover = 0.5;

	std::vector<BlastRes> forward = {hit("q", "high_identity", 100.0, 100, 100, 200), hit("q", "longer", 95.0, 120, 120, 200)};
	std::vector<BlastRes> reverse = {forward[1], forward[0]};
	double forwardBest = 0.0, reverseBest = 0.0;
	filterBlastPrimary(forward, &opt, forwardBest);
	filterBlastPrimary(reverse, &opt, reverseBest);
	check(forward.size() == 1 && reverse.size() == 1 && forward.front().Sbj == "high_identity" &&
		reverse.front().Sbj == "high_identity",
		"a longer alignment must not displace the highest-identity anchor, in either input order");
	check(forwardBest == reverseBest && forwardBest == 100.0, "order must not change the anchor identity");

	std::vector<BlastRes> partialReference = {hit("q", "close", 99.0, 200, 200, 250), hit("q", "far", 90.0, 245, 245, 250)};
	filterBlastPrimary(partialReference, &opt, bestID);
	check(partialReference.size() == 1 && partialReference.front().Sbj == "close",
		"a higher-coverage, lower-identity hit must not remove the closest hit");
	std::vector<BlastRes> relativeCoverage = {hit("q", "anchor", 99.0, 200, 200, 250),
		hit("q", "enough", 98.9, 170, 170, 250), hit("q", "short", 98.9, 169, 169, 250)};
	filterBlastPrimary(relativeCoverage, &opt, bestID);
	check(relativeCoverage.size() == 2 && relativeCoverage[1].Sbj == "enough",
		"the relative filter must keep exactly 85% of the anchor's query coverage");
	std::vector<BlastRes> boundary = {hit("q", "a", 97, 100, 100, 100), hit("q", "b", 97, 85, 85, 100)};
	filterBlastPrimary(boundary, &opt, bestID);
	check(boundary.size() == 2, "an alignment at the exact 85% coverage cutoff must be retained");
	std::vector<BlastRes> legacyRelative = {
		BlastRes("q\tlong\t99\t200\t0\t0\t1\t200\t1\t200\t1e-50\t300", 0),
		BlastRes("q\tenough\t99\t170\t0\t0\t1\t170\t1\t170\t1e-40\t250", 0),
		BlastRes("q\tshort\t99\t169\t0\t0\t1\t169\t1\t169\t1e-40\t250", 0)};
	filterBlastPrimary(legacyRelative, &opt, bestID);
	check(legacyRelative.size() == 2 && legacyRelative[1].Sbj == "enough",
		"legacy rows must fall back to 85% of the anchor's alignment length");

	std::vector<BlastRes> absolute = {hit("q", "cov50", 99, 100, 100, 200), hit("q", "cov49", 99, 98, 98, 200),
		hit("q", "len75", 99, 75, 75, 100), hit("q", "len74", 99, 74, 74, 100)};
	opt.BLfilter = true;
	filterBlastPrimary(absolute, &opt, bestID);
	std::vector<std::string> kept;
	for (const auto& h : absolute) { kept.push_back(h.Sbj); }
	check(std::find(kept.begin(), kept.end(), "cov49") == kept.end() &&
		std::find(kept.begin(), kept.end(), "len74") == kept.end() &&
		std::find(kept.begin(), kept.end(), "len75") != kept.end(),
		"minimum coverage and alignment length must be inclusive cutoffs");

	std::vector<BlastRes> window = {hit("q", "best", 99, 100, 100, 100), hit("q", "inside", 98.75, 100, 100, 100),
		hit("q", "outside", 98.74, 100, 100, 100)};
	filterBlastPrimary(window, &opt, bestID);
	check(window.size() == 2 && window[1].Sbj == "inside", "the identity window below 99% must be 0.25 points, inclusive");

	std::vector<BlastRes> duplicates = {hit("q", "a", 100, 50, 50, 100), hit("q", "a", 99, 100, 100, 100)};
	filterBlastPrimary(duplicates, &opt, bestID);
	check(duplicates.size() == 1 && duplicates[0].perID == 99,
		"a short duplicate must not hide an eligible alignment for the same subject");
	duplicates = {hit("q", "a", 100, 100, 100, 300), hit("q", "a", 99, 100, 100, 100)};
	filterBlastPrimary(duplicates, &opt, bestID);
	check(duplicates.size() == 1 && duplicates[0].perID == 99,
		"a low-coverage duplicate must not hide an eligible alignment for the same subject");
	duplicates = {hit("q", "a", 97.9, 100, 100, 100), hit("q", "b", 98, 100, 100, 100), hit("q", "a", 98, 100, 100, 100)};
	filterBlastPrimary(duplicates, &opt, bestID);
	check(duplicates.size() == 2 && duplicates[0].Sbj == "a" && duplicates[0].perID == 98 && duplicates[1].Sbj == "b",
		"deduplication must keep first-occurrence order and the best alignment per subject");

	opt.BLfilter = false;
	std::vector<BlastRes> prefiltered = {hit("q", "a", 99, 100, 100, 100), hit("q", "b", 80, 30, 30, 100),
		hit("q", "a", 98, 100, 100, 100)};
	filterBlastPrimary(prefiltered, &opt, bestID);
	check(prefiltered.size() == 2 && prefiltered[1].Sbj == "b",
		"-no_bl_filter must skip length, coverage and identity filters but still deduplicate");
}

void bestHitMode() {
	currentGroup = "best-hit mode";
	options opt = makeOptions({"-reportBestHit"});
	double bestID = 0.0;
	for (bool filtered : {true, false}) {
		opt.BLfilter = filtered;
		std::vector<BlastRes> tied = {hit("q", "b", 99, 100, 100, 100), hit("q", "a", 99, 100, 100, 100),
			hit("q", "c", 98.9, 100, 100, 100), hit("q", "a", 99, 100, 100, 100)};
		filterBlastPrimary(tied, &opt, bestID);
		check(tied.size() == 2 && tied[0].Sbj == "b" && tied[1].Sbj == "a",
			"best-hit mode must keep every equally ranked subject once, with or without primary filtering");
		std::vector<BlastRes> ranked = {hit("q", "partial", 99, 90, 90, 100), hit("q", "full", 99, 100, 100, 100)};
		filterBlastPrimary(ranked, &opt, bestID);
		check(ranked.size() == 1 && ranked[0].Sbj == "full", "best-hit mode must prefer higher query coverage at equal identity");
	}
	check(std::all_of(opt.idThr.begin(), opt.idThr.end(), [](double t) { return t == 1.0; }),
		"best-hit mode must replace all identity cutoffs with 1%");
}

void consensus() {
	currentGroup = "consensus";
	bool singular = false;
	std::vector<TaxObj*> oneGap = {lineageOf({"B", "P", "C", "", "F", "G", "S"})};
	std::vector<TaxObj*> twoGaps = {lineageOf({"B", "P", "C", "", "F", "G", "S"}), lineageOf({"B", "P", "C", "", "F", "G", "S"})};
	TaxObj* single = LCAcore(oneGap, singular, 0.9, 7);
	check(singular, "one subject must be reported as a singular hit");
	TaxObj* pair = LCAcore(twoGaps, singular, 0.9, 7);
	check(!singular, "two subjects are not a singular hit");
	check(single->SavedTaxs == pair->SavedTaxs && pair->get(3) == __unkwnTax && pair->get(6) == "S",
		"identical lineages with a missing rank must give the same result for one or several subjects");
	delete single; delete pair;
	deleteAll(oneGap); deleteAll(twoGaps);

	std::vector<TaxObj*> conflict = {lineageOf({"B", "P", "C", "", "F1", "G1", "S1"}), lineageOf({"B", "P", "C", "O", "F2", "G2", "S2"})};
	TaxObj* conflicted = LCAcore(conflict, singular, 0.9, 7);
	check(conflicted->depth == 4 && conflicted->get(3) == "O" && conflicted->get(4) == __unkwnTax,
		"a hit with a missing rank must still vote against conflicting deeper ranks");
	delete conflicted;
	deleteAll(conflict);

	std::vector<TaxObj*> parentTest;
	for (int i = 0; i < 8; ++i) { parentTest.push_back(lineageOf({"A"})); }
	parentTest.push_back(lineageOf({"", "child-a"}));
	parentTest.push_back(lineageOf({"", "child-b"}));
	TaxObj* lineage = LCAcore(parentTest, singular, 0.9, 2);
	check(lineage->depth == 1 && lineage->get(0) == "A" && lineage->get(1) == __unkwnTax,
		"conflicting children below a missing parent must stop the LCA");
	delete lineage;
	parentTest[9]->set(1, nameOf("child-a"));
	lineage = LCAcore(parentTest, singular, 0.9, 2);
	check(lineage->depth == 2 && lineage->get(1) == "child-a", "matching children below a missing parent must be used");
	delete lineage;
	deleteAll(parentTest);

	std::vector<TaxObj*> votes;
	for (int i = 0; i < 9; ++i) { votes.push_back(lineageOf({"A"})); }
	votes.push_back(lineageOf({"B"}));
	lineage = LCAcore(votes, singular, 0.9, 7);
	check(lineage->get(0) == "A", "exactly the -LCAfrac share of known votes must be accepted");
	delete lineage;
	votes[8]->set(0, nameOf("B"));
	lineage = LCAcore(votes, singular, 0.9, 7);
	check(lineage->depth == 0, "less than the -LCAfrac share of known votes must stop the LCA");
	delete lineage;
	deleteAll(votes);
	// In floating point, 25 * 0.56 is slightly more than 14.
	for (int i = 0; i < 25; ++i) { votes.push_back(lineageOf({i < 14 ? "A" : "B"})); }
	lineage = LCAcore(votes, singular, 0.56, 7);
	check(lineage->get(0) == "A", "an exact -LCAfrac share must be accepted despite rounding");
	delete lineage;
	deleteAll(votes);

	std::vector<TaxObj*> tie = {lineageOf({"B"}), lineageOf({"A"}), lineageOf({"B"}), lineageOf({"A"})};
	lineage = LCAcore(tie, singular, 0.5, 7);
	check(lineage->get(0) == "A", "equally supported taxa must be resolved lexicographically");
	delete lineage;
	deleteAll(tie);

	std::vector<TaxObj*> unknown = {lineageOf({}), lineageOf({})};
	lineage = LCAcore(unknown, singular, 0.9, 7);
	check(lineage->depth == 0, "lineages without names must give an empty result");
	delete lineage;
	deleteAll(unknown);

	TaxObj uncertain("k__Bacteria;s__Genus sp. X", 7, false, false, testNames);
	std::vector<TaxObj*> alone = {&uncertain};
	lineage = LCAcore(alone, singular, 0.9, 7);
	check(lineage->speciesUncertain && lineage->get(6) == "Genus sp. X", "a single hit must keep its species and uncertainty");
	delete lineage;
}

void identityLimits() {
	currentGroup = "identity limits";
	const std::vector<double> thresholds = makeOptions().idThr;
	const std::pair<double, int> cases[] = {{0.0, 1}, {77.99, 1}, {78.0, 2}, {92.99, 4}, {93.0, 5}, {96.99, 6}, {97.0, 7}, {100.0, 7}};
	for (const auto& [identity, depth] : cases) {
		TaxObj tax(7);
		for (int i = 0; i < 7; ++i) { tax.set(i, nameOf("rank" + std::to_string(i))); }
		tax.limitByIdentity(identity, thresholds);
		check(tax.depth == depth, "each rank must require at least its identity cutoff");
	}
	TaxObj shallow(7);
	shallow.set(0, nameOf("B")); shallow.set(1, nameOf("P"));
	shallow.limitByIdentity(100.0, thresholds);
	check(shallow.depth == 2, "identity must not deepen a lineage beyond its known ranks");
}

void outputFormatting() {
	currentGroup = "output formatting";
	TaxObj output(7);
	for (int i = 0; i < 7; ++i) { output.set(i, nameOf("rank" + std::to_string(i))); }
	check(output.getWriteString(std::vector<double>(7, 0.0)) == "rank0\trank1\trank2\trank3\trank4\trank5\trank6",
		"no Hit2DB or %ID column may be written unless requested");
	output.setHitDB("subject-read");
	output.setRepID(true);
	output.perID = 98.5;
	check(output.getWriteString(std::vector<double>(7, 0.0)) ==
		"rank0\trank1\trank2\trank3\trank4\trank5\trank6\tsubject-read\t98.500000",
		"hit read and identity must follow all normal taxonomy columns");
	output.perID = 94.0;
	check(output.getWriteString(makeOptions().idThr) == "rank0\trank1\trank2\trank3\trank4\t?\t?\tsubject-read\t94.000000",
		"ranks above the reported identity must be written as ?");
	TaxObj partial(7);
	partial.set(0, nameOf("B"));
	partial.set(2, nameOf("C"));
	partial.setHitDB(__unkwnTax);
	check(partial.getWriteString(std::vector<double>(7, 0.0)) == "B\t?\tC\t?\t?\t?\t?\t?",
		"missing ranks, ranks beyond the depth and an unknown Hit2DB must be written as ?");
}

TaxObj assignment(std::initializer_list<const char*> ranks, double identity, const char* subject) {
	TaxObj* built = lineageOf(ranks);
	TaxObj result(built);
	delete built;
	result.perID = identity;
	result.setHitDB(subject);
	result.setRepID(true);
	return result;
}

void databaseSelection() {
	currentGroup = "database selection";
	TaxObj tieA = assignment({"r0", "r1", "r2", "r3", "r4", "r5", "r6"}, 98.5, "z-subject");
	TaxObj tieB = assignment({"r0", "r1", "r2", "r3", "r4", "r5", "r6"}, 98.5, "a-subject");
	check(tieA.evalAcpyTax(&tieB) && tieA.getHitDB() == "a-subject",
		"equal multi-database assignments must select hit reads deterministically");

	TaxObj truncated(7);
	for (int i = 0; i < 7; ++i) { truncated.set(i, nameOf("rank" + std::to_string(i))); }
	truncated.limitByIdentity(94, makeOptions().idThr);
	TaxObj moreKnown(7);
	for (int i = 0; i < 6; ++i) { moreKnown.set(i, nameOf("rank" + std::to_string(i))); }
	check(truncated.depth == 5 && truncated.evalAcpyTax(&moreKnown) && truncated.depth == 6,
		"taxonomy beyond the identity cutoff must not count toward database selection");

	std::vector<TaxObj> candidates = {
		assignment({"B", "P", "C", "O", "F", "", "S"}, 97, "deeper"),        // 6 known ranks, depth 7
		assignment({"B", "P", "C", "O", "F", "G"}, 99, "shallower"),        // 6 known ranks, depth 6
		assignment({"B", "P", "C", "O", "F", "G"}, 98, "lower-identity"),
		assignment({"B", "P", "C", "O", "F"}, 100, "fewer-known"),
		assignment({"B", "P", "C", "O", "F", "", "S"}, 97, "another-subject"),
		assignment({"B", "P", "C", "O", "F", "", "R"}, 97, "other-lineage"),
	};
	candidates[4].speciesUncertain = true;
	for (size_t i = 0; i < candidates.size(); ++i) {
		for (size_t j = 0; j < candidates.size(); ++j) {
			TaxObj first(&candidates[i]);
			first.evalAcpyTax(&candidates[j]);
			TaxObj second(&candidates[j]);
			second.evalAcpyTax(&candidates[i]);
			check(first.getWriteString(std::vector<double>(7, 0.0)) == second.getWriteString(std::vector<double>(7, 0.0)),
				"database selection must not depend on database order");
		}
	}
	TaxObj best(&candidates[0]);
	for (const auto& candidate : candidates) { best.evalAcpyTax(&candidate); }
	check(best.getHitDB() == "other-lineage",
		"more known ranks, then depth, then identity, then a certain species, then the smaller lineage must win");

	TaxObj hiddenZ = assignment({"r0", "r1", "r2", "r3", "r4", "zzz", "r6"}, 94, "a-subject");
	TaxObj hiddenA = assignment({"r0", "r1", "r2", "r3", "r4", "aaa", "r6"}, 94, "b-subject");
	hiddenZ.limitByIdentity(94, makeOptions().idThr);
	hiddenA.limitByIdentity(94, makeOptions().idThr);
	check(hiddenA.evalAcpyTax(&hiddenZ) && hiddenA.getHitDB() == "a-subject",
		"ranks hidden by the identity cutoff must not decide between databases");
}

void optionParsing() {
	currentGroup = "options";
	options depth8 = makeOptions({"-tdep", "8"});
	check(depth8.taxDepth == 8 && depth8.idThr.size() == 8 && depth8.Taxlvls.size() == 8 && depth8.Taxlvls[7] == "Strain",
		"-tdep must resize taxonomy thresholds and labels consistently");
	options depth10 = makeOptions({"-tdep", "10"});
	check(depth10.Taxlvls[8] == "Level9" && depth10.Taxlvls[9] == "Level10" && depth10.idThr[9] == 97,
		"ranks beyond Strain must be labelled LevelN with a 97% default cutoff");
	options custom = makeOptions({"-id", "99.5,95,93,91,88,78,0", "-cover", "0.75", "-LCAfrac", "0.8", "-minAlignLen", "60"});
	check(custom.idThr[6] == 99.5 && custom.idThr[5] == 95 && custom.idThr[0] == 0,
		"-id must list cutoffs from the deepest rank back to Domain, with decimals");
	check(custom.minCover == 0.75 && custom.LCAfract == 0.8 && custom.minAliLen == 60,
		"numeric options must be stored as given");
	options threads = makeOptions({"-t", "16", "-SLVfmt", "-f", "bl8"});
	check(threads.taxDepth == 7 && threads.input_format == "bl8" && threads.numThr == 16,
		"-t must set the thread count alongside compatibility options");
	check(makeOptions().numThr == 1 && makeOptions({"-t", "0"}).numThr == 1 && makeOptions({"-t", "four"}).numThr == 1 &&
		makeOptions({"-t", "1000"}).numThr == 256, "invalid -t values must mean one thread, and large ones are capped");
	check(makeOptions().TaxLvl2string() == "Domain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies",
		"the output header must list the configured ranks");
}

void readers() {
	currentGroup = "readers";
	{
		TempFile m8("# comment\r\n"
			"qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tqlen\r\n"
			"q1\ta\t99\t100\t0\t0\t1\t100\t1\t100\t100\r\n"
			"q1\tb\t98\t100\t0\t0\t1\t100\t1\t100\t100\n"
			"q1 malformed\n"
			"q2\ta\t97\t100\t0\t0\t1\t100\t1\t100\t100");
		Silence quiet;
		BlastReader reader(m8.path, "bl8");
		const auto first = reader.getResBatch();
		const auto second = reader.getResBatch();
		const auto end = reader.getResBatch();
		check(first.size() == 2 && first[0].Query == "q1" && first[1].Sbj == "b" && second.size() == 1 &&
			second[0].Query == "q2" && end.empty() && reader.queriesRead() == 2,
			"the mapping reader must group contiguous query blocks across comments, headers and CRLF");
		check(quiet.sink.str().find("line 5") != std::string::npos, "a skipped record must be reported by line number");
	}
	const std::string alpha = "k__Bacteria;p__P;c__C;o__O;f__F;g__Genus;s__Genus alpha";
	const std::string shallow = "k__Bacteria;p__P";
	for (bool shallowFirst : {true, false}) {
		TempFile tax("# reference taxonomy\n" +
			(shallowFirst ? "dup\t" + shallow + "\ndup\t" + alpha : "dup\t" + alpha + "\ndup\t" + shallow) +
			"\nother\t" + alpha + "\textra column\n");
		Silence quiet;
		RefTax reference(tax.path, testNames, 7, false, true);
		const RefTax::Lineage* dup = reference.find("dup");
		const RefTax::Lineage* other = reference.find("other");
		check(dup != nullptr && dup->depth == 7, "duplicate identifiers must keep the more complete lineage");
		check(other != nullptr && reference.rankName(*other, 6) == "Genus alpha",
			"columns after the lineage must not become part of rank names");
		check(reference.find("absent") == nullptr && reference.size() == 2, "unknown identifiers must not be found");
	}
}

void orderedPipeline() {
	currentGroup = "ordered pipeline";
	for (int threads : {1, 4}) {
		int next = 0;
		std::vector<int> consumed;
		auto produce = [&next](int& item) {
			if (next == 200) { return false; }
			item = next++;
			return true;
		};
		auto slowSquare = [](int& item, int& result) {
			std::this_thread::sleep_for(std::chrono::microseconds((item * 7) % 5 * 40));
			result = item * item;
		};
		auto consume = [&consumed](int& result) { consumed.push_back(result); };
		runOrdered<int, int>(threads, produce, slowSquare, consume);
		bool ordered = consumed.size() == 200;
		for (size_t i = 0; ordered && i < consumed.size(); ++i) { ordered = consumed[i] == static_cast<int>(i * i); }
		check(ordered, "results must be consumed in production order");

		next = 0;
		consumed.clear();
		int code = 0;
		try {
			runOrdered<int, int>(threads, produce, [](int& item, int& result) {
				if (item == 57) { throw FatalError{7, "work"}; }
				result = item;
			}, consume);
		} catch (const FatalError& error) { code = error.code; }
		check(code == 7 && consumed.size() == 57 && consumed.back() == 56,
			"a failing item must be reported after every earlier result, and no later one");

		next = 0;
		consumed.clear();
		code = 0;
		try {
			runOrdered<int, int>(threads, [&next](int& item) {
				if (next == 30) { throw FatalError{8, "produce"}; }
				item = next++;
				return true;
			}, [](int& item, int& result) { result = item; }, consume);
		} catch (const FatalError& error) { code = error.code; }
		check(code == 8 && consumed.size() == 30, "a failing producer must still deliver the items produced before it");
	}
}

void taxonomyLoading() {
	currentGroup = "taxonomy loading";
	std::string text = "# reference taxonomy\r\n\n";
	for (int i = 0; i < 1500; ++i) {
		const std::string genus = "Genus" + std::to_string(i % 37);
		text += "ref" + std::to_string(i) + "\tk__Bacteria;p__P" + std::to_string(i % 5) + ";c__C" + std::to_string(i % 11) +
			";g__" + genus + ";s__" + genus + " species" + std::to_string(i % 13) + " strain" +
			(i % 7 == 0 ? "\textra" : "") + (i % 10 == 0 ? "\r" : "") + "\n";
		if (i % 97 == 0) { text += "malformed line without a tab\n"; }
		if (i % 101 == 0) { text += "ref" + std::to_string(i / 2) + "\tk__Bacteria;p__Shallow\n"; }
		if (i % 103 == 0) {
			text += "ref" + std::to_string(i / 3) + "\tk__Bacteria;p__P0;c__C0;o__Deeper;f__F;g__Genus0;s__Genus0 zeta\n";
		}
	}
	text += "last\tk__Bacteria;p__NoNewline";
	TempFile file(text);

	auto render = [](const RefTax& tax, const std::string& id) {
		const RefTax::Lineage* lineage = tax.find(id);
		if (lineage == nullptr) { return std::string("missing"); }
		std::string rendered = std::to_string(lineage->depth) + (lineage->speciesUncertain ? "u" : "c");
		for (int rank = 0; rank < 7; ++rank) { rendered += ";" + tax.rankName(*lineage, rank); }
		return rendered;
	};
	std::vector<std::string> expected;
	std::string expectedLog;
	const std::pair<int, size_t> configurations[] = {{1, size_t(4) << 20}, {1, 64}, {4, 64}, {4, 7}, {3, 1000}};
	for (const auto& [threads, blockSize] : configurations) {
		std::vector<std::string> lineages;
		std::string log;
		{
			Silence quiet;
			NamePool names; // one per load, so distinctNames counts this load only
			RefTax tax(file.path, names, 7, false, true, threads, blockSize);
			lineages.push_back(render(tax, "last"));
			for (int i = 0; i < 1500; ++i) { lineages.push_back(render(tax, "ref" + std::to_string(i))); }
			lineages.push_back(std::to_string(tax.size()) + "/" + std::to_string(tax.distinctNames()));
			if (expected.empty()) {
				TaxObj parsed("k__Bacteria;p__P1;c__C1;g__Genus1;s__Genus1 species1 strain", 7, false, false, testNames);
				const RefTax::Lineage* ref1 = tax.find("ref1");
				bool same = ref1 != nullptr && ref1->depth == parsed.depth;
				for (int rank = 0; same && rank < 7; ++rank) { same = tax.rankName(*ref1, rank) == parsed.get(rank); }
				check(same, "stored lineages must match the TaxObj parser");
				check(tax.rankName(*tax.find("ref0"), 3) == "Deeper" && tax.rankName(*tax.find("ref34"), 3) == "Deeper" &&
					tax.find("ref50")->depth == 7, "duplicates must keep the lineage with more known ranks");
				check(tax.rankName(*tax.find("ref70"), 6) == "Genus33 species5" && tax.rankName(*tax.find("last"), 1) == "NoNewline",
					"CRLF, extra columns and a final line without newline must be handled");
				check(tax.size() == 1501 && tax.find("malformed line without a tab") == nullptr,
					"every valid identifier must be stored once");
			}
			log = quiet.sink.str();
		}
		if (expected.empty()) {
			expected = lineages;
			expectedLog = log;
			check(log.find("skipped 16 taxonomy line(s)") != std::string::npos && log.find("first at line 4)") != std::string::npos &&
				log.find("on 215 taxonomy line(s)") != std::string::npos && log.find("first at line 3)") != std::string::npos,
				"warnings must count lines and name the first one");
		} else {
			check(lineages == expected, "lineages and duplicates must not depend on threads or block size");
			check(log == expectedLog, "reports must not depend on threads or block size");
		}
	}

	TempFile small("a\tk__Bacteria;p__P;c__C;o__O;f__F;g__G;s__G alpha\nu\tk__Bacteria;g__G;s__G sp. X\n");
	Silence quiet;
	RefTax tax(small.path, testNames, 7, false, true);
	HitLineage limited(tax, *tax.find("a"));
	limited.limitByIdentity(96.0, makeOptions().idThr);
	HitLineage masked(tax, *tax.find("a"));
	masked.speciesMasked = true;
	check(limited.depth == 6 && limited.get(6) == __unkwnTax && limited.get(5) == "G" &&
		masked.get(6) == __unkwnTax && masked.get(5) == "G" && masked.depth == 7,
		"hit lineages must hide ranks beyond their identity and a suppressed species");
	check(HitLineage(tax, *tax.find("u")).speciesUncertain() && !limited.speciesUncertain(),
		"hit lineages must report the reference's species uncertainty");
}

// A taxonomy read from its index must be the one read from its text, and the index must be rebuilt when the file
// or the settings that shape the taxonomy change.
void taxonomyIndex() {
	currentGroup = "taxonomy index";
	std::string text = "# reference taxonomy\n";
	for (int i = 0; i < 400; ++i) {
		const std::string genus = "Genus" + std::to_string(i % 23);
		text += "ref" + std::to_string(i % 350) + "\tk__Bacteria;p__P" + std::to_string(i % 4) + ";g__" +
			(i % 11 == 0 ? "uncultured" : genus) + ";s__" + genus + (i % 9 == 0 ? " sp. X" : " species" + std::to_string(i % 7)) +
			(i % 50 == 0 ? "\textra" : "") + "\n";
		if (i % 97 == 0) { text += "malformed line\n"; }
	}
	TempFile file(text);
	const std::string index = RefTax::indexPath(file.path);
	struct Remove {
		std::string path;
		~Remove() { std::error_code ignored; std::filesystem::remove_all(path, ignored); }
	} removeIndex{index};

	auto render = [](const RefTax& tax) {
		std::string rendered = std::to_string(tax.size()) + "/" + std::to_string(tax.distinctNames());
		for (int i = 0; i <= 350; ++i) {
			const std::string id = i < 350 ? "ref" + std::to_string(i) : "new";
			const RefTax::Lineage* lineage = tax.find(id);
			rendered += "\n" + id + (lineage == nullptr ? " missing" :
				" " + std::to_string(lineage->depth) + (lineage->speciesUncertain ? "u" : "c"));
			for (int rank = 0; lineage != nullptr && rank < lineage->depth; ++rank) { rendered += ";" + tax.rankName(*lineage, rank); }
		}
		return rendered;
	};
	// A pool of its own for each load, holding `before` first, as a pool shared with an earlier database would.
	auto load = [&file, &render](int depth, bool checkUnknown, bool useIndex, std::string& log,
		std::initializer_list<const char*> before = {}) {
		NamePool pool;
		for (const char* name : before) { pool.intern(name); }
		Silence quiet;
		const RefTax tax(file.path, pool, depth, false, checkUnknown, 2, 1000, useIndex);
		log = quiet.sink.str();
		return render(tax);
	};
	auto withoutIndex = [](const std::string& log) {
		std::istringstream lines(log);
		std::string kept, line;
		while (std::getline(lines, line)) {
			if (line.find("axonomy index") == std::string::npos) { kept += line + "\n"; }
		}
		return kept;
	};
	auto loadsLike = [&](int depth, bool checkUnknown, const char* expectedLine, const char* message,
		std::initializer_list<const char*> before = {}) {
		std::string plainLog, log;
		const std::string plain = load(depth, checkUnknown, false, plainLog, before);
		const std::string indexed = load(depth, checkUnknown, true, log, before);
		check(indexed == plain && withoutIndex(log) == plainLog && log.find(expectedLine) != std::string::npos, message);
		if (log.find(expectedLine) == std::string::npos) { std::fprintf(stderr, "  log:\n%s", log.c_str()); }
	};
	auto indexBytes = [&index] {
		std::ifstream in(index, std::ios::binary);
		return std::string(std::istreambuf_iterator<char>(in), {});
	};

	std::string log;
	load(7, true, false, log);
	check(!std::filesystem::exists(index), "no index may be written without useIndex");
	loadsLike(7, true, "Wrote taxonomy index", "the first load must write the index and give the text's taxonomy and reports");
	check(std::filesystem::is_regular_file(index), "the index must be written beside the taxonomy");
	loadsLike(7, true, "Read taxonomy index", "the next load must read the index and give the same taxonomy and reports");
	// Names that a pool holds already have other IDs there than in a pool of their own.
	const std::initializer_list<const char*> other = {"Genus5 species3", "Other", "P2", "Bacteria"};
	loadsLike(7, true, "Read taxonomy index", "an index must be read into a pool that holds other names", other);
	const std::string fresh = indexBytes();
	std::filesystem::remove(index);
	loadsLike(7, true, "Wrote taxonomy index", "an index must be written from a pool that holds other names", other);
	check(indexBytes() == fresh, "the index must not depend on the other names in the pool");
	loadsLike(8, true, "does not match", "another depth must rebuild the index");
	loadsLike(8, false, "does not match", "-no_taxDB_filter must rebuild the index");
	loadsLike(8, false, "Read taxonomy index", "the rebuilt index must be read");

	std::filesystem::last_write_time(file.path, std::filesystem::last_write_time(file.path) + std::chrono::seconds(2));
	loadsLike(8, false, "does not match", "a newer taxonomy file must rebuild the index");
	std::ofstream(file.path, std::ios::binary | std::ios::app) << "new\tk__Bacteria;p__New\n";
	loadsLike(8, false, "does not match", "a changed taxonomy file must rebuild the index");
	check(load(8, false, true, log).find("\nnew 2c;Bacteria;New") != std::string::npos,
		"the rebuilt index must hold the changed taxonomy");

	// A damaged index is rebuilt, never read.
	const std::string good = indexBytes();
	auto damaged = [&](std::string content, const char* message) {
		std::ofstream(index, std::ios::binary | std::ios::trunc) << content;
		loadsLike(8, false, "does not match", message);
	};
	damaged(good.substr(0, good.size() - 1), "a truncated index must be rebuilt");
	damaged(good + "x", "an index with trailing bytes must be rebuilt");
	std::string changed = good;
	changed[48] = static_cast<char>(changed[48] + 1); // the count of names
	damaged(changed, "an index whose counts do not add up must be rebuilt");
	changed = good;
	changed[changed.size() - 13] = '\xff'; // the last slot's lineage, before the 8-byte end mark
	damaged(changed, "an index with a slot past the lineages must be rebuilt");
	damaged("", "an empty index must be rebuilt");

	// An index that cannot be written leaves the taxonomy as read from text, and no temporary file behind.
	std::filesystem::remove(index);
	std::filesystem::create_directory(index);
	loadsLike(8, false, "could not write the taxonomy index", "an index that cannot be written must be reported");
	bool leftovers = false;
	for (const auto& entry : std::filesystem::directory_iterator(std::filesystem::path(index).parent_path())) {
		leftovers = leftovers || entry.path().string().rfind(index + ".tmp", 0) == 0;
	}
	check(!leftovers, "a failed index write must remove its temporary file");
}

void endToEnd() {
	currentGroup = "LCA driver";
	TempFile tax("aa\tk__Bacteria;p__P;c__C;o__O;f__F;g__Genus;s__Genus beta\n"
		"bb\tk__Bacteria;p__P;c__C;o__O;f__F;g__Genus;s__Genus alpha\n"
		"cc\tk__Bacteria;p__P;c__C;o__O;f__F;g__Genus;s__Genus alpha\n"
		"sp\tk__Bacteria;p__P;c__C;o__O;f__F;g__Genus;s__Genus sp. X\n");
	std::unique_ptr<RefTax> reference;
	{
		Silence quiet;
		reference = std::make_unique<RefTax>(tax.path, testNames, 7, false, true);
	}
	options opt = makeOptions({"-showHitRead", "-reportID"});
	auto run = [&](std::vector<BlastRes> hits) {
		TaxObj* result = LCA(hits, reference.get(), &opt);
		std::string text = result ? result->Subj + "\t" + result->getWriteString(opt.idThr) : "";
		delete result;
		return text;
	};
	check(run({hit("q", "aa", 99, 100, 100, 100), hit("q", "bb", 99, 100, 100, 100)}) ==
		"q\tBacteria\tP\tC\tO\tF\tGenus\t?\t?\t99.000000", "disagreeing species must stop at genus with Hit2DB ?");
	check(run({hit("q", "bb", 99, 100, 100, 100)}) == "q\tBacteria\tP\tC\tO\tF\tGenus\tGenus alpha\tbb\t99.000000",
		"a single subject must be named in Hit2DB");
	check(run({hit("q", "sp", 99, 100, 100, 100), hit("q", "bb", 99, 100, 100, 100)}) ==
		"q\tBacteria\tP\tC\tO\tF\tGenus\tGenus alpha\t?\t99.000000", "a certain species must outvote an uncertain one");
	check(run({hit("q", "sp", 99, 100, 100, 100)}) == "q\tBacteria\tP\tC\tO\tF\tGenus\tGenus sp. X\tsp\t99.000000",
		"an uncertain species alone must be kept");
	check(run({hit("q", "bb", 97.2, 100, 100, 100), hit("q", "cc", 96.6, 100, 100, 100)}) ==
		"q\tBacteria\tP\tC\tO\tF\tGenus\t?\t?\t96.900000", "the mean identity must limit the reported depth");
	check(run({hit("q", "bb", 50, 10, 10, 100)}).empty(), "a query without eligible hits must have no result");
	std::vector<BlastRes> single = {hit("q", "bb", 99, 100, 100, 100)};
	std::unique_ptr<TaxObj> held(LCA(single, reference.get(), &opt));
	check(&held->get(6) == nameOf("Genus alpha"),
		"assignments must point into the name pool, so that they outlive their RefTax with several databases");

	opt = makeOptions({"-showHitRead", "-reportBestHit"});
	const std::vector<BlastRes> tied = {hit("q", "cc", 99, 100, 100, 100), hit("q", "aa", 99, 100, 100, 100),
		hit("q", "bb", 99, 100, 100, 100), hit("q", "zz-worse", 98, 100, 100, 100)};
	check(run(tied) == "q\tBacteria\tP\tC\tO\tF\tGenus\t?\taa", "tied best hits of different species must give genus");
	opt.LCAfract = 0.6;
	check(run(tied) == "q\tBacteria\tP\tC\tO\tF\tGenus\tGenus alpha\tbb",
		"Hit2DB must name the first tied subject that supports the consensus");
}

void matrixCounts() {
	currentGroup = "abundance matrix";
	options opt = makeOptions();
	Matrix matrix(7, opt.Taxlvls, true);
	TaxObj* first = lineageOf({"B", "P", "C", "O", "F", "G", "S"});
	first->setHitDB("ref1");
	TaxObj* second = lineageOf({"B", "P", "C", "O", "F", "G"});
	second->setHitDB(__unkwnTax);
	matrix.add(first);
	matrix.add(first);
	matrix.add(second);
	check(matrix.rowIDs[0].at("B") == 3 && matrix.rowIDs[6].at("B;P;C;O;F;G;S") == 2 &&
		matrix.rowIDs[6].at("B;P;C;O;F;G;?") == 1 && matrix.rowIDs[7].at("B;P;C;O;F;G;S;ref1") == 2,
		"lineage counts must accumulate per rank, with ? for missing names");
	const auto paths = matrix.outputPaths("out.tsv");
	check(paths.size() == 8 && paths.front() == "out.tsv_Domain" && paths.back() == "out.tsv_Hit2DB",
		"one count file per rank plus Hit2DB must be named after the output path");
	delete first;
	delete second;
}

void emptyConversion() {
	currentGroup = "hit conversion";
	options opt = makeOptions();
	std::vector<BlastRes> empty;
	double meanIdentity = 123.0;
	const auto noTax = BlastToTax(empty, nullptr, &opt, meanIdentity);
	check(noTax.empty() && meanIdentity == 0.0, "empty hit conversion must reset the mean without dividing by zero");
}

lca::Taxonomy mockTaxonomy() {
	lca::Taxonomy t;
	t.AddLineage("d__Bacteria;p__Mockota;c__Mockia;o__Mockales;f__Mockaceae;g__Mockella;s__Mockella alpha");  // 0
	t.AddLineage("d__Bacteria;p__Mockota;c__Mockia;o__Mockales;f__Mockaceae;g__Mockella;s__Mockella beta");   // 1
	t.AddLineage("d__Bacteria;p__Mockota;c__Mockia;o__Mockales;f__Mockaceae;g__Otherella;s__Otherella one");  // 2
	t.AddLineage("d__Bacteria;p__Mockota;c__Mockia;o__Mockales;f__Mockaceae;g__Mockella");                    // 3: genus record
	t.AddLineage("k__Alveolata;p__Dinoflagellata");                                                             // 4: a decoy
	t.AddLineage("d__Bacteria;p__Mockota;c__Mockia;o__Mockales;f__Otheraceae;g__Farella;s__Farella far");      // 5
	t.AddLineage("d__Bacteria;p__Mockota;c__Mockia;f__Gapaceae;g__Gapella;s__Gapella one");                    // 6: no order
	t.AddLineage("d__Bacteria;p__Mockota;c__Mockia;f__Gapaceae;g__Gapella;s__Gapella two");                    // 7: no order
	t.AddLineage("d__Bacteria;p__Mockota;c__Mockia;o__Mockales;f__Mockaceae;g__Mockella;s__Mockella sp. X");   // 8
	return t;
}

lca::Hit libraryHit(uint32_t subject, double identity, uint32_t align = 150, uint32_t queryLength = 150, uint32_t covered = 150) {
	return {subject, identity, align, queryLength, covered};
}

std::string libraryLineage(const lca::Assignment& a, const lca::Taxonomy& t) {
	std::string text;
	for (int rank = 0; rank < a.depth; ++rank) { text += (rank == 0 ? "" : ";") + t.Name(a.name_ids[rank]); }
	return text;
}

// protal's cases for its stand-in (tests/test_RrnaLca.cpp there), with LCA's results where its rules differ.
void library() {
	currentGroup = "library";
	const lca::Taxonomy t = mockTaxonomy();
	const lca::Options opt;
	lca::Scratch s;
	auto assign = [&](std::vector<lca::Hit> hits, const lca::Options& o) { return lca::Assign(hits, t, o, s); };

	check(t.Size() == 9 && t.Name(t.NameId(0, 0)) == "Bacteria" && t.NameId(0, 0) == t.NameId(5, 0) && t.NameId(3, 6) == 0 &&
		t.Name(t.NameId(4, 0)) == "Alveolata" && t.NameId(4, 2) == 0 && t.Name(0) == "?" && t.DistinctNames() == 20,
		"lineages must be stored by index with interned names, k__ read as d__");
	check(t.SpeciesUncertain(8) && !t.SpeciesUncertain(0), "placeholder species must be marked uncertain");
	const std::string_view names[] = {"Bacteria", " P ", "unclassified", "", "F", "G", "G alpha strain 1"};
	lca::Taxonomy added;
	added.Add(names);
	check(added.NameId(0, 2) == 0 && added.Name(added.NameId(0, 1)) == "P" && added.Name(added.NameId(0, 6)) == "G alpha",
		"names added by rank must follow the taxonomy file's naming rules");

	lca::Assignment a = assign({libraryHit(0, 99.3)}, opt);
	check(a.depth == 7 && libraryLineage(a, t).ends_with(";Mockella alpha") && a.hits_used == 1 && a.mean_identity == 99.3 &&
		!a.species_uncertain, "one hit at high identity must give its species");

	check(assign({libraryHit(0, 96.0)}, opt).depth == 6 && assign({libraryHit(0, 93.5)}, opt).depth == 5 &&
		assign({libraryHit(0, 80.0)}, opt).depth == 2, "identity cutoffs must limit the depth");
	lca::Options yarza;
	yarza.min_identity = {98.7, 94.5, 86.5, 82, 78.5, 75, 0};
	check(assign({libraryHit(0, 93.5)}, yarza).depth == 5 && assign({libraryHit(0, 80.0)}, yarza).depth == 3 &&
		assign({libraryHit(0, 76.0)}, yarza).depth == 2 && assign({libraryHit(0, 98.7)}, yarza).depth == 7,
		"cutoffs must be read deepest first, inclusive");

	check(libraryLineage(assign({libraryHit(0, 99), libraryHit(1, 99)}, opt), t) ==
		"Bacteria;Mockota;Mockia;Mockales;Mockaceae;Mockella", "two species of one genus must give the genus");
	check(assign({libraryHit(0, 99), libraryHit(1, 99), libraryHit(2, 98.98)}, opt).depth == 5,
		"another genus inside the identity window must give the family");
	a = assign({libraryHit(0, 99), libraryHit(1, 99), libraryHit(2, 96)}, opt);
	check(a.depth == 6 && a.hits_used == 2, "hits below the identity window must be left out");
	check(assign({libraryHit(0, 99), libraryHit(2, 98.75)}, opt).depth == 5 && assign({libraryHit(0, 99), libraryHit(2, 98.74)}, opt).depth == 7,
		"the identity window at 99% must be 0.25 points, inclusive");
	check(lca::IdentityWindow(80) == 1.5 && lca::IdentityWindow(96.99) == 1.5 && lca::IdentityWindow(97) == 1.0 &&
		lca::IdentityWindow(98) == 0.75 && lca::IdentityWindow(99) == 0.25 && lca::IdentityWindow(99.5) == 0.15f &&
		lca::IdentityWindow(100) == 0.05f, "the identity window must narrow in steps toward 100%");
	lca::Taxonomy many;
	for (int i = 0; i < 9; ++i) { many.AddLineage("d__B;p__P;c__C;o__O;f__F;g__G;s__G same"); }
	many.AddLineage("d__B;p__P;c__C;o__O;f__F;g__G;s__G other");
	std::vector<lca::Hit> ten;
	for (uint32_t i = 0; i < 10; ++i) { ten.push_back(libraryHit(i, 99)); }
	lca::Options strict = opt;
	strict.lca_frac = 0.95;
	a = lca::Assign(ten, many, opt, s);
	check(a.depth == 7 && many.Name(a.name_ids[6]) == "G same" && lca::Assign(ten, many, strict, s).depth == 6,
		"a rank must need lca_frac of the known votes");
	lca::Taxonomy split;
	std::vector<lca::Hit> twentyFive;
	for (uint32_t i = 0; i < 25; ++i) {
		split.AddLineage(i < 14 ? "d__A" : "d__B");
		twentyFive.push_back(libraryHit(i, 99));
	}
	lca::Options exactShare = opt;
	exactShare.lca_frac = 0.56; // 25 * 0.56 is slightly more than 14 in floating point
	check(lca::Assign(twentyFive, split, exactShare, s).depth == 1, "an exact lca_frac share must be accepted despite rounding");

	check(libraryLineage(assign({libraryHit(0, 99), libraryHit(3, 99)}, opt), t).ends_with(";Mockella alpha"),
		"a hit without a species must abstain there");
	a = assign({libraryHit(4, 99.5), libraryHit(0, 97)}, opt);
	check(a.depth == 2 && libraryLineage(a, t) == "Alveolata;Dinoflagellata" && a.hits_used == 1,
		"a decoy as the best hit must be assigned alone");
	// Not in the stand-in: a rank no hit names stays unknown, and deeper ranks are still decided.
	a = assign({libraryHit(6, 99), libraryHit(7, 99)}, opt);
	check(libraryLineage(a, t) == "Bacteria;Mockota;Mockia;?;Gapaceae;Gapella" && a.name_ids[3] == 0,
		"ranks below one that no hit names must still be assigned");
	// Not in the stand-in: placeholder species abstain beside a certain one.
	a = assign({libraryHit(8, 99), libraryHit(0, 99)}, opt);
	check(libraryLineage(a, t).ends_with(";Mockella alpha") && !a.species_uncertain, "a certain species must outvote a placeholder");
	a = assign({libraryHit(8, 99)}, opt);
	check(libraryLineage(a, t).ends_with(";Mockella sp. X") && a.species_uncertain, "a placeholder species alone must be kept, as uncertain");

	a = assign({libraryHit(0, 99, 60, 150, 60), libraryHit(5, 95, 150, 400, 150)}, opt);
	check(a.depth == 0 && a.hits_used == 0, "short hits and hits covering too little of the query must be dropped");
	a = assign({libraryHit(0, 99), libraryHit(0, 98.99), libraryHit(1, 99)}, opt);
	check(a.hits_used == 2 && a.depth == 6 && a.mean_identity == 99, "one subject aligned twice must be one vote, its best");
	// 85 of 100 covered bases is kept; as fractions of 113 query bases, 85/113 < 100/113 * 0.85.
	check(assign({libraryHit(0, 99, 100, 113, 100), libraryHit(2, 99, 85, 113, 85)}, opt).depth == 5 &&
		assign({libraryHit(0, 99, 100, 113, 100), libraryHit(2, 99, 84, 113, 84)}, opt).depth == 7,
		"a hit covering exactly 85% of the anchor's query bases must be kept");
	lca::Options best = opt;
	best.best_hit = true;
	// The stand-in kept every hit at the anchor's identity; -reportBestHit keeps those also equal in coverage and length.
	a = assign({libraryHit(2, 99), libraryHit(0, 99, 150, 150, 140), libraryHit(1, 98.9)}, best);
	check(a.hits_used == 1 && libraryLineage(a, t).ends_with(";Otherella one"), "best_hit must keep the hits ranking equal to the anchor");
	a = assign({libraryHit(2, 99), libraryHit(0, 99), libraryHit(1, 98.9)}, best);
	check(a.hits_used == 2 && libraryLineage(a, t).ends_with(";Mockaceae"), "tied best hits must be combined by consensus");
	check(assign({libraryHit(99, 99.9), libraryHit(0, 99)}, opt).depth == 7, "hits of subjects outside the taxonomy must be ignored");

	std::vector<lca::Hit> large, small = {libraryHit(8, 99), libraryHit(0, 99)};
	for (uint32_t subject : {0u, 1u, 2u, 3u, 5u, 6u, 7u, 8u, 0u, 2u}) { large.push_back(libraryHit(subject, 98.5)); }
	lca::Scratch warm;
	lca::Assign(large, t, opt, warm);
	const size_t allocations = heapAllocations;
	for (int i = 0; i < 100; ++i) {
		lca::Assign(large, t, opt, warm);
		lca::Assign(small, t, best, warm);
	}
	check(heapAllocations == allocations, "Assign must not allocate once its scratch has grown");
}

// GTDB appends capital letters after an underscore (_A, _AB, _CF); they are never placeholder markers, while the
// lower case markers of PR2 and NCBI names are.
void gtdbSuffixes() {
	currentGroup = "GTDB suffixes";
	check(lca::kNamingRulesVersion == 2, "the naming rules version must be raised for the GTDB suffix rule");
	const char* certain[] = {"Streptococcus mitis_CF", "Streptococcus oralis_CF", "Bacillus_CF chungangensis",
		"Prevotella copri_A", "Bacillus_A cereus_AB", "Collinsella sp900762015", "Escherichia coli",
		"Foo bar_SP", "Foo_AFF baz", "Foo bar_SPP"};
	for (const char* name : certain) {
		TaxObj tax(std::string("k__Bacteria;g__G;s__") + name, 7, false, false, testNames);
		check(!lca::species_name_is_placeholder(name) && !tax.speciesUncertain,
			(std::string("a GTDB name must not be a placeholder: ") + name).c_str());
	}
	const char* placeholders[] = {"Streptococcus sp.", "Streptococcus sp. X1", "Lachnospiraceae bacterium",
		"Pichia sp. NRRL Y-12824", "Pichia_sp._NRRL_Y-12824", "Sordariomycetes_sp.", "Candida cf. albicans",
		"Candida_cf._albicans", "Mockella aff. alpha", "uncultured bacterium",
		"Foo_SP bar_sp", "Foo_CF_cf baz", "SP bar", "Foo_BACTERIUM"};
	for (const char* name : placeholders) {
		TaxObj tax(std::string("k__Bacteria;g__G;s__") + name, 7, false, false, testNames);
		check(lca::species_name_is_placeholder(name) && tax.speciesUncertain,
			(std::string("a placeholder name must stay a placeholder: ") + name).c_str());
	}
	check(!lca::species_name_is_placeholder("") && !lca::species_name_is_placeholder("_CF") &&
		!lca::species_name_is_placeholder("Genus_") && lca::species_name_is_placeholder("sp"),
		"empty names and a lone suffix must not be placeholders");

	const std::string family = "d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;";
	lca::Taxonomy t;
	t.AddLineage(family + "s__Streptococcus mitis_CF");  // 0
	t.AddLineage(family + "s__Streptococcus oralis");    // 1
	check(!t.SpeciesUncertain(0) && !t.SpeciesUncertain(1), "S. mitis_CF must be a certain species");
	const lca::Options opt;
	lca::Scratch s;
	std::vector<lca::Hit> hits = {libraryHit(0, 99), libraryHit(1, 99)};
	lca::Assignment a = lca::Assign(hits, t, opt, s);
	check(a.depth == 6 && libraryLineage(a, t).ends_with(";Streptococcus") && !a.species_uncertain,
		"two certain species at one identity must give the genus, not S. oralis");
	hits = {libraryHit(0, 99)};
	a = lca::Assign(hits, t, opt, s);
	check(a.depth == 7 && libraryLineage(a, t).ends_with(";Streptococcus mitis_CF") && !a.species_uncertain,
		"a single hit on S. mitis_CF must name it, as certain");
}

BlastRes subjectHit(int subject, double identity, int alignmentLength, int queryStop, int queryLength) {
	const std::string id = "s" + std::to_string(subject);
	if (queryLength > 0) { return hit("q", id.c_str(), identity, alignmentLength, queryStop, queryLength); }
	std::ostringstream line; // a 12-column BLAST record, without the query length
	line << "q\t" << id << '\t' << identity << '\t' << alignmentLength << "\t0\t0\t1\t" << queryStop << "\t1\t"
		<< alignmentLength << "\t1e-30\t100";
	return BlastRes(line.str(), 0);
}

// lca::Assign must give the program's assignment: random queries against a random taxonomy, in every mode.
void libraryMatchesProgram() {
	currentGroup = "library and program";
	std::mt19937 random(20261023);
	auto pick = [&random](int n) { return std::uniform_int_distribution<int>(0, n - 1)(random); };
	auto chance = [&random](double p) { return std::uniform_real_distribution<double>(0.0, 1.0)(random) < p; };

	// Names from a binary tree, so that subjects share upper ranks and differ below, with missing and unknown ranks,
	// shortened lineages, placeholder species and species names with extra words.
	const char* const prefixes[] = {"d__", "p__", "c__", "o__", "f__", "g__", "s__"};
	const int subjects = 80;
	std::string text;
	lca::Taxonomy library;
	std::vector<int> leaf(subjects);
	for (int subject = 0; subject < subjects; ++subject) {
		std::string lineage;
		const int ranks = chance(0.2) ? 1 + pick(7) : 7;
		int node = 0;
		for (int rank = 0; rank < 7; ++rank) {
			node = node * 2 + pick(2);
			if (rank >= ranks || chance(0.08)) { continue; }
			std::string name = "N" + std::to_string(rank) + "x" + std::to_string(node);
			if (rank == 6) {
				const std::string genus = "G" + std::to_string(node / 2);
				const int kind = pick(8);
				name = kind == 0 ? genus + " sp. " + std::to_string(subject) : kind == 1 ? genus + " bacterium" :
					kind == 2 ? genus + " e" + std::to_string(node) + " strain " + std::to_string(subject) :
					kind == 3 ? "Candidatus " + genus + " e" + std::to_string(node) + " x" : genus + " e" + std::to_string(node);
			}
			if (chance(0.05)) { name = chance(0.5) ? "?" : "uncultured"; }
			lineage += (lineage.empty() ? "" : ";") + std::string(prefixes[rank]) + name;
		}
		if (lineage.empty()) { lineage = "d__N0x0"; }
		leaf[subject] = node;
		text += "s" + std::to_string(subject) + "\t" + lineage + "\n";
		library.AddLineage(lineage);
	}
	TempFile file(text);
	std::unique_ptr<RefTax> program;
	{
		Silence quiet;
		program = std::make_unique<RefTax>(file.path, testNames, 7, false, true);
	}

	// Each query hits subjects next to each other in the tree, so that its hits share some ranks.
	std::vector<int> byLeaf(subjects);
	std::iota(byLeaf.begin(), byLeaf.end(), 0);
	std::sort(byLeaf.begin(), byLeaf.end(), [&leaf](int a, int b) { return leaf[a] < leaf[b]; });
	struct Query {
		std::vector<BlastRes> records;
		std::vector<lca::Hit> hits;
	};
	std::vector<Query> queries(4000);
	const double tops[] = {100, 99.5, 99.2, 99, 98.4, 98, 97.5, 97, 96, 94.5, 93, 90, 86.5, 82, 78.5, 76};
	const double steps[] = {0, 0, 0.05, 0.15, 0.25, 0.5, 0.75, 1, 1.5, 2, 3};
	for (Query& query : queries) {
		const int first = pick(subjects);
		const int count = 1 + pick(10);
		const bool legacy = chance(0.1);
		const int queryLength = 100 + pick(100);
		const double top = chance(0.7) ? tops[pick(16)] : 75 + pick(2500) / 100.0;
		for (int i = 0; i < count; ++i) {
			const int subject = byLeaf[std::min(subjects - 1, first + pick(6))];
			const double identity = i == 0 ? top : top - (chance(0.5) ? steps[pick(11)] : pick(300) / 100.0);
			const int alignmentLength = chance(0.1) ? 40 + pick(40) : 80 + pick(120);
			const int length = legacy ? 0 : chance(0.05) ? queryLength + 50 : queryLength;
			const int queryStop = legacy ? alignmentLength : chance(0.8) ? length - pick(20) : 1 + pick(length);
			const BlastRes record = subjectHit(subject, identity, alignmentLength, queryStop, length);
			if (record.fail) { continue; }
			query.records.push_back(record);
			query.hits.push_back({static_cast<uint32_t>(subject), record.perID, static_cast<uint32_t>(record.alLen),
				static_cast<uint32_t>(record.queryCoverageKnown ? record.queryLength : 0),
				static_cast<uint32_t>(record.queryCoverageKnown ? record.queryCovered : 0)});
		}
	}

	struct Mode {
		std::vector<std::string> args;
		lca::Options library;
	};
	lca::Options protal, half, unanimous, unfiltered, bestHit;
	protal.min_identity = {98.7, 94.5, 86.5, 82, 78.5, 75, 0};
	protal.min_align_len = 50;
	half.lca_frac = 0.5;
	unanimous.lca_frac = 1;
	unfiltered.no_filters = true;
	bestHit.best_hit = true;
	bestHit.min_identity.assign(7, 1.0);
	lca::Options bestUnfiltered = bestHit;
	bestUnfiltered.no_filters = true;
	const Mode modes[] = {
		{{}, lca::Options()},
		{{"-id", "98.7,94.5,86.5,82,78.5,75,0", "-minAlignLen", "50", "-cover", "0.5"}, protal},
		{{"-LCAfrac", "0.5"}, half},
		{{"-LCAfrac", "1"}, unanimous},
		{{"-no_bl_filter"}, unfiltered},
		{{"-reportBestHit"}, bestHit},
		{{"-reportBestHit", "-no_bl_filter"}, bestUnfiltered},
	};

	auto programLineage = [](const TaxObj* tax) {
		std::string text;
		for (int rank = 0; rank < 7; ++rank) { text += (rank == 0 ? "" : ";") + tax->get(rank); }
		return text;
	};
	auto libraryRanks = [&library](const lca::Assignment& a) {
		std::string text;
		for (int rank = 0; rank < 7; ++rank) {
			text += (rank == 0 ? "" : ";") + (rank < a.depth && a.name_ids[rank] != 0 ? library.Name(a.name_ids[rank]) : __unkwnTax);
		}
		return text;
	};
	lca::Scratch scratch;
	int mismatches = 0, placed = 0, species = 0, consensus = 0, gaps = 0, uncertain = 0;
	for (const Mode& mode : modes) {
		options opt = makeOptions(mode.args);
		for (const Query& query : queries) {
			std::vector<BlastRes> records = query.records; // LCA() leaves the kept subjects' hits
			std::unique_ptr<TaxObj> expected(LCA(records, program.get(), &opt));
			const lca::Assignment got = lca::Assign(query.hits, library, mode.library, scratch);
			const std::string want = expected ? programLineage(expected.get()) : "none";
			const std::string have = got.hits_used > 0 ? libraryRanks(got) : "none";
			bool same = want == have && records.size() == got.hits_used && (!expected || expected->perID == got.mean_identity);
			if (expected && records.size() == 1) {
				same = same && got.species_uncertain == (expected->speciesUncertain && expected->get(6) != __unkwnTax);
			}
			if (!same && mismatches++ < 5) {
				std::fprintf(stderr, "  mode %zu, %zu hits: program %s (%zu subjects, %f), library %s (%u subjects, %f)\n",
					static_cast<size_t>(&mode - modes), query.records.size(), want.c_str(), records.size(),
					expected ? expected->perID : 0.0, have.c_str(), got.hits_used, got.mean_identity);
			}
			placed += got.depth > 0;
			species += got.depth == 7;
			consensus += got.depth > 0 && got.hits_used > 1;
			uncertain += got.species_uncertain;
			for (int rank = 0; rank + 1 < got.depth; ++rank) {
				if (got.name_ids[rank] == 0) { ++gaps; break; }
			}
		}
	}
	check(mismatches == 0, "the library must assign every query as the program does");
	check(placed > 20000 && species > 5000 && consensus > 8000 && gaps > 3000 && uncertain > 1000,
		"the random queries must reach every rule");
}

// Subjects 0 to 11 for every naming rule, then a tree of subjects whose neighbours share their upper ranks, with
// missing ranks, shortened lineages and placeholder species. Depth 8 at most.
lca::Taxonomy namingRulesTaxonomy(int depth) {
	lca::Taxonomy t(depth);
	t.AddLineage("d__Bacteria;p__Firmicutes;c__Clostridia;o__Eubacteriales;f__Lachnospiraceae;g__Blautia;"
		"s__Blautia obeum ATCC 29174");                                                            // 0: shortened
	t.AddLineage("k__Bacteria;p__Firmicutes;c__uncultured;o__?;f__;g__Blautia;s__Blautia sp. X"); // 1: unknown, placeholder
	t.AddLineage("Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;"
		"Escherichia coli O157:H7");                                                               // 2: unprefixed
	t.AddLineage("d__Bacteria;p__Proteobacteria;s__Candidatus Liberibacter asiaticus strain Y");  // 3: Candidatus
	t.AddLineage("d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria");                         // 4: shorter than the depth
	t.AddLineage("d__Bacteria;g__Prevotella;s__uncultured bacterium clone X");                     // 5: placeholder
	t.AddLineage("D__Archaea;P__Euryarchaeota;s__Methanobrevibacter smithii;t__strain PS");        // 6: a strain
	t.AddLineage(" d__ Bacteria ; p__ Bacteroidota ; Bacteroidia ;; uncultured bacterium ;");       // 7: padded, mixed
	t.AddLineage("");                                                                              // 8: no names
	const std::string_view byRank[] = {"Bacteria", " Bacillota ", "unclassified", "", "Lachnospiraceae", "Unknown",
		"Roseburia intestinalis L1-82"};
	t.Add(byRank);                                                                                 // 9
	const std::string_view shortRanks[] = {"Archaea", "?"};
	t.Add(shortRanks);                                                                             // 10
	const std::string_view placeholder[] = {"Bacteria", "Bacillota", "Bacilli", "Lactobacillales", "Lactobacillaceae",
		"Lactobacillus", "Lactobacillus cf. casei"};
	t.Add(placeholder);                                                                            // 11

	std::mt19937 random(static_cast<unsigned>(depth));
	for (int subject = 0; subject < 300; ++subject) {
		const unsigned leaf = static_cast<unsigned>(subject) * 256 / 300;
		const int ranks = random() % 5 == 0 ? 1 + static_cast<int>(random() % depth) : depth;
		std::string lineage;
		for (int rank = 0; rank < ranks; ++rank) {
			const unsigned node = leaf >> (7 - rank);
			std::string name = random() % 20 == 0 ? "" : "N" + std::to_string(rank) + "x" + std::to_string(node);
			if (rank == lca::kSpeciesRank && !name.empty()) {
				const std::string genus = "G" + std::to_string(node / 2);
				name = random() % 3 == 0 ? genus + " sp. " + std::to_string(subject) : genus + " e" + std::to_string(node);
			}
			lineage += std::string(rank == 0 ? "" : ";") + "dpcofgst"[rank] + "__" + name;
		}
		t.AddLineage(lineage);
	}
	return t;
}

bool sameTaxonomy(const lca::Taxonomy& a, const lca::Taxonomy& b) {
	if (a.Depth() != b.Depth() || a.Size() != b.Size() || a.DistinctNames() != b.DistinctNames()) { return false; }
	for (uint32_t id = 0; id <= a.DistinctNames(); ++id) {
		if (a.Name(id) != b.Name(id)) { return false; }
	}
	for (uint32_t subject = 0; subject < a.Size(); ++subject) {
		if (a.SpeciesUncertain(subject) != b.SpeciesUncertain(subject)) { return false; }
		for (int rank = 0; rank < a.Depth(); ++rank) {
			if (a.NameId(subject, rank) != b.NameId(subject, rank)) { return false; }
		}
	}
	return true;
}

bool sameParts(const lca::TaxonomyParts& a, const lca::TaxonomyParts& b) {
	return a.depth == b.depth && a.names == b.names && a.ids == b.ids && a.species_uncertain == b.species_uncertain;
}

// A Taxonomy rebuilt from its parts must be the one they were taken from, for Assign and for Add.
void taxonomyParts() {
	currentGroup = "taxonomy parts";
	std::mt19937 random(20261024);
	const double identities[] = {100, 99.5, 99.2, 99, 98.4, 97, 96, 94, 92, 88, 80};
	lca::Scratch scratch;
	for (int depth : {7, 4, 8}) {
		lca::Taxonomy original = namingRulesTaxonomy(depth);
		auto name = [&original](uint32_t subject, int rank) { return original.Name(original.NameId(subject, rank)); };
		if (depth >= 7) {
			check(name(0, 6) == "Blautia obeum" && original.NameId(1, 2) == 0 && original.NameId(1, 3) == 0 &&
				original.NameId(1, 4) == 0 && original.SpeciesUncertain(1) && name(1, 6) == "Blautia sp. X" &&
				name(2, 6) == "Escherichia coli" && original.NameId(2, 0) == original.NameId(0, 0) &&
				name(3, 6) == "Candidatus Liberibacter asiaticus" && original.SpeciesUncertain(5) &&
				name(5, 6) == "uncultured bacterium clone X" && name(7, 2) == "Bacteroidia" && original.NameId(7, 4) == 0 &&
				name(9, 1) == "Bacillota" && original.NameId(9, 5) == 0 && name(9, 6) == "Roseburia intestinalis" &&
				original.NameId(10, 1) == 0 && original.SpeciesUncertain(11) && (depth == 7 || name(6, 7) == "strain PS"),
				"the round trip's subjects must reach every naming rule");
		}

		const lca::TaxonomyParts parts = original.Parts();
		bool layout = parts.depth == depth && parts.names.size() == original.DistinctNames() + 1 &&
			parts.ids.size() == original.Size() * depth && parts.species_uncertain.size() == original.Size();
		for (uint32_t id = 0; layout && id < parts.names.size(); ++id) { layout = parts.names[id] == original.Name(id); }
		for (uint32_t subject = 0; layout && subject < original.Size(); ++subject) {
			layout = parts.species_uncertain[subject] == (original.SpeciesUncertain(subject) ? 1 : 0);
			for (int rank = 0; layout && rank < depth; ++rank) {
				layout = parts.ids[subject * depth + rank] == original.NameId(subject, rank);
			}
		}
		check(layout && parts.names[0] == "?", "Parts must hold the names by id and depth ids per subject");

		std::string error;
		std::optional<lca::Taxonomy> built = lca::Taxonomy::FromParts(lca::TaxonomyParts(parts), error);
		check(built.has_value(), "the parts of a Taxonomy must be accepted");
		if (!built) {
			std::fprintf(stderr, "  depth %d: %s\n", depth, error.c_str());
			continue;
		}
		lca::Taxonomy rebuilt = std::move(*built);
		check(sameTaxonomy(original, rebuilt), "FromParts(Parts()) must give the same names, ids and placeholder marks");

		lca::Options opt;
		opt.depth = depth;
		opt.min_identity = {99, 97, 95, 93, 91, 88, 78, 0};
		int differences = 0, placed = 0, uncertain = 0;
		for (int set = 0; set < 300; ++set) {
			std::vector<lca::Hit> hits;
			const uint32_t first = static_cast<uint32_t>(random() % (original.Size() - 8));
			const uint32_t queryLength = random() % 5 == 0 ? 0 : 150;
			for (int i = 0, count = 1 + static_cast<int>(random() % 8); i < count; ++i) {
				const uint32_t align = static_cast<uint32_t>(60 + random() % 3 * 45);
				const uint32_t covered = queryLength == 0 ? 0 : std::min(align, queryLength) - static_cast<uint32_t>(random() % 2 * 10);
				hits.push_back({first + static_cast<uint32_t>(random() % 8), identities[random() % 11], align, queryLength, covered});
			}
			for (bool noFilters : {false, true}) {
				for (bool bestHit : {false, true}) {
					opt.no_filters = noFilters;
					opt.best_hit = bestHit;
					const lca::Assignment a = lca::Assign(hits, original, opt, scratch);
					const lca::Assignment b = lca::Assign(hits, rebuilt, opt, scratch);
					differences += a.depth != b.depth || a.name_ids != b.name_ids || a.species_uncertain != b.species_uncertain ||
						a.hits_used != b.hits_used || a.mean_identity != b.mean_identity;
					placed += a.depth >= 3;
					uncertain += a.species_uncertain;
				}
			}
		}
		check(differences == 0, "Assign must give the same assignments on the rebuilt Taxonomy");
		check(placed > 200 && (depth < 7 || uncertain > 20), "the random hit sets must reach deep ranks and placeholder species");

		const uint32_t known = static_cast<uint32_t>(original.DistinctNames());
		const std::string_view added[] = {"Newdomain", " Firmicutes"};
		for (lca::Taxonomy* t : {&original, &rebuilt}) {
			t->AddLineage("d__Bacteria;p__Newphylum;c__Clostridia");
			t->Add(added);
		}
		const uint32_t s = static_cast<uint32_t>(rebuilt.Size()) - 2;
		check(rebuilt.NameId(s, 0) == rebuilt.NameId(0, 0) && rebuilt.NameId(s, 1) == known + 1 &&
			rebuilt.NameId(s, 2) == rebuilt.NameId(0, 2) && rebuilt.NameId(s + 1, 0) == known + 2 &&
			rebuilt.NameId(s + 1, 1) == rebuilt.NameId(0, 1) && rebuilt.DistinctNames() == known + 2,
			"Add on a rebuilt Taxonomy must reuse the ids of its names and give new names new ids");
		check(sameTaxonomy(original, rebuilt), "Add must give the same ids on the rebuilt Taxonomy as on the original");
	}

	const lca::TaxonomyParts valid = namingRulesTaxonomy(7).Parts();
	auto refusal = [&valid](auto change, std::initializer_list<const char*> words, const char* message) {
		lca::TaxonomyParts parts = valid;
		change(parts);
		const lca::TaxonomyParts given = parts;
		std::string error;
		bool refused = !lca::Taxonomy::FromParts(std::move(parts), error);
		for (const char* word : words) { refused = refused && error.find(word) != std::string::npos; }
		check(refused && sameParts(parts, given), message); // refused parts are left unchanged
		if (!refused) { std::fprintf(stderr, "  error: \"%s\"\n", error.c_str()); }
	};
	using Parts = lca::TaxonomyParts;
	refusal([](Parts& p) { p.depth = 0; }, {"depth 0"}, "depth 0 must be refused");
	refusal([](Parts& p) { p.depth = lca::kMaxDepth + 1; }, {"depth 65"}, "a depth past kMaxDepth must be refused");
	refusal([](Parts& p) { p.names.clear(); }, {"name 0"}, "parts without names must be refused");
	refusal([](Parts& p) { p.names[0] = "unknown"; }, {"name 0", "unknown"}, "a name 0 other than ? must be refused");
	refusal([](Parts& p) { p.names[2].clear(); }, {"name 2 is empty"}, "an empty name must be refused");
	refusal([](Parts& p) { p.names[3] = p.names[1]; }, {"name 3", "repeats name 1"}, "a repeated name must be refused");
	refusal([](Parts& p) { p.names.back() = "?"; }, {"repeats name 0"}, "a second ? must be refused");
	refusal([](Parts& p) { p.ids.pop_back(); }, {"ids for"}, "ids that are not depth per subject must be refused");
	refusal([](Parts& p) { p.species_uncertain.push_back(0); }, {"ids for"}, "a subject without ids must be refused");
	refusal([](Parts& p) { p.ids[9] = static_cast<uint32_t>(p.names.size()); }, {"index 9"},
		"an id past the names must be refused");

	// About ribotal's reference (1.77 M subjects, 274 K names): FromParts must take a small part of the time of Add.
	const uint32_t pool[] = {3, 60, 600, 2500, 12000, 45000, 140000};
	std::vector<std::vector<std::string>> names(7);
	for (int rank = 0; rank < 7; ++rank) {
		for (uint32_t i = 0; i < pool[rank]; ++i) {
			names[rank].push_back(rank == lca::kSpeciesRank ? "Genus" + std::to_string(i % 45000) + " species" + std::to_string(i) :
				"Rank" + std::to_string(rank) + "name" + std::to_string(i));
		}
	}
	lca::Taxonomy large;
	std::array<std::string_view, 7> lineage;
	const auto addStart = std::chrono::steady_clock::now();
	for (uint32_t subject = 0; subject < 1000000; ++subject) {
		for (int rank = 0; rank < 7; ++rank) {
			uint32_t mixed = subject * 2654435761u + static_cast<uint32_t>(rank) * 0x9E3779B9u;
			mixed ^= mixed >> 15;
			lineage[rank] = names[rank][mixed % pool[rank]];
		}
		large.Add(lineage);
	}
	const double addSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - addStart).count();
	lca::TaxonomyParts largeParts = large.Parts();
	std::string error;
	const auto rebuildStart = std::chrono::steady_clock::now();
	std::optional<lca::Taxonomy> rebuilt = lca::Taxonomy::FromParts(std::move(largeParts), error);
	const double rebuildSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - rebuildStart).count();
	std::printf("  %zu subjects, %zu names: Add %.3f s, FromParts %.3f s\n", large.Size(), large.DistinctNames(),
		addSeconds, rebuildSeconds);
	check(rebuilt && sameTaxonomy(large, *rebuilt), "a large Taxonomy must be rebuilt unchanged");
	check(rebuildSeconds * 5 < addSeconds, "FromParts must take a small part of the time that Add takes");
}
}

int main() {
	taxonomyParsing();
	speciesNames();
	mappingRecords();
	hitRanking();
	primaryFilter();
	bestHitMode();
	consensus();
	identityLimits();
	outputFormatting();
	databaseSelection();
	optionParsing();
	readers();
	orderedPipeline();
	taxonomyLoading();
	taxonomyIndex();
	endToEnd();
	matrixCounts();
	emptyConversion();
	library();
	gtdbSuffixes();
	libraryMatchesProgram();
	taxonomyParts();

	if (failures != 0) {
		std::cerr << failures << " regression check(s) failed\n";
		return 1;
	}
	std::cout << "All regression checks passed\n";
	return 0;
}
