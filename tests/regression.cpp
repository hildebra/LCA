#include "../LCAimpl.h"
#include "../Matrix.h"
#include "../Pipeline.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <memory>

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

TaxObj* lineageOf(std::initializer_list<const char*> ranks) {
	TaxObj* tax = new TaxObj(7);
	int rank = 0;
	for (const char* name : ranks) {
		if (*name != '\0') { tax->set(rank, name); }
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
	TaxObj prefixed("s__Genus sp.;k__Bacteria;g__Genus", 7, false, false);
	check(prefixed.get(0) == "Bacteria", "k__ must map to Domain even after s__");
	check(prefixed.get(5) == "Genus", "g__ must map to Genus even after s__");
	check(prefixed.get(6) == "Genus sp.", "s__ must map to Species");
	check(prefixed.speciesUncertain, "a species containing sp. must be marked uncertain");

	TaxObj unprefixed("Bacteria;P;;O;F;G;Genus species", 7, false, false);
	check(unprefixed.get(2) == __unkwnTax && unprefixed.get(3) == "O" && unprefixed.depth == 7,
		"empty positions in unprefixed lineages must keep later ranks in place");
	TaxObj padded(" k__ Bacteria ; p__ P ;", 7, false, false);
	check(padded.get(0) == "Bacteria" && padded.get(1) == "P" && padded.depth == 2,
		"whitespace around names and a trailing separator must be ignored");
	TaxObj upper("D__Bacteria;P__Firmicutes", 7, false, false);
	check(upper.get(0) == "Bacteria" && upper.get(1) == "Firmicutes", "rank prefixes must be case-insensitive");
	TaxObj kept("k__Bacteria;g__uncultured", 7, false, true);
	TaxObj filtered("k__Bacteria;g__uncultured", 7, false, false);
	check(kept.get(5) == "uncultured" && filtered.get(5) == __unkwnTax && filtered.depth == 1,
		"-no_taxDB_filter must keep names that are otherwise treated as unknown");
	TaxObj unknownGenus("k__Bacteria;p__P;g__Unknown;s__Unidentified", 7, false, false);
	check(unknownGenus.get(5) == __unkwnTax && unknownGenus.get(6) == __unkwnTax && unknownGenus.depth == 2,
		"unknown and unidentified rank names must be treated as missing");
	TaxObj strain8("k__Bacteria;s__Genus species;t__strain X", 8, false, false);
	TaxObj strain7("k__Bacteria;s__Genus species;t__strain X", 7, false, false);
	check(strain8.get(7) == "strain X" && strain7.depth == 7, "t__ must be read only when Strain is configured");
}

void speciesNames() {
	currentGroup = "species names";
	TaxObj candidatus("k__Bacteria;p__P;c__C;o__O;f__F;g__Liberibacter;s__Candidatus Liberibacter asiaticus strain X",
		7, false, false);
	check(candidatus.get(6) == "Candidatus Liberibacter asiaticus",
		"Candidatus species names must retain their three-part species name");
	TaxObj serotype("k__Bacteria;s__Escherichia coli O157:H7", 7, false, false);
	TaxObj subspecies("k__Bacteria;s__Lactobacillus casei subsp. casei", 7, false, false);
	check(serotype.get(6) == "Escherichia coli" && subspecies.get(6) == "Lactobacillus casei",
		"certain species names must be shortened to genus and species");
	const char* placeholders[] = {"Lachnospiraceae bacterium 3_1_46FAA", "Methanomassiliicoccales archaeon RumEn M1",
		"human gut metagenome", "Blautia sp", "Bacteroides_sp._XY1", "Lactobacillus cf. casei",
		"Streptococcus aff. mitis", "unidentified Prevotella", "Ruminococcus spp."};
	for (const char* name : placeholders) {
		TaxObj tax(std::string("k__Bacteria;g__G;s__") + name, 7, false, false);
		check(tax.speciesUncertain && tax.get(6) == name, "placeholder species names must be marked uncertain, unshortened");
	}
	const char* species[] = {"Mycobacterium tuberculosis", "Bacteroides_A sp000123456", "Blautia obeum",
		"Spirochaeta spectabilis"};
	for (const char* name : species) {
		TaxObj tax(std::string("k__Bacteria;g__G;s__") + name, 7, false, false);
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
	parentTest[9]->set(1, "child-a");
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
	votes[8]->set(0, "B");
	lineage = LCAcore(votes, singular, 0.9, 7);
	check(lineage->depth == 0, "less than the -LCAfrac share of known votes must stop the LCA");
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

	TaxObj uncertain("k__Bacteria;s__Genus sp. X", 7, false, false);
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
		for (int i = 0; i < 7; ++i) { tax.set(i, "rank" + std::to_string(i)); }
		tax.limitByIdentity(identity, thresholds);
		check(tax.depth == depth, "each rank must require at least its identity cutoff");
	}
	TaxObj shallow(7);
	shallow.set(0, "B"); shallow.set(1, "P");
	shallow.limitByIdentity(100.0, thresholds);
	check(shallow.depth == 2, "identity must not deepen a lineage beyond its known ranks");
}

void outputFormatting() {
	currentGroup = "output formatting";
	TaxObj output(7);
	for (int i = 0; i < 7; ++i) { output.set(i, "rank" + std::to_string(i)); }
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
	partial.set(0, "B");
	partial.set(2, "C");
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
	for (int i = 0; i < 7; ++i) { truncated.set(i, "rank" + std::to_string(i)); }
	truncated.limitByIdentity(94, makeOptions().idThr);
	TaxObj moreKnown(7);
	for (int i = 0; i < 6; ++i) { moreKnown.set(i, "rank" + std::to_string(i)); }
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
		RefTax reference(tax.path, 7, false, true);
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
			RefTax tax(file.path, 7, false, true, threads, blockSize);
			lineages.push_back(render(tax, "last"));
			for (int i = 0; i < 1500; ++i) { lineages.push_back(render(tax, "ref" + std::to_string(i))); }
			lineages.push_back(std::to_string(tax.size()) + "/" + std::to_string(tax.distinctNames()));
			if (expected.empty()) {
				TaxObj parsed("k__Bacteria;p__P1;c__C1;g__Genus1;s__Genus1 species1 strain", 7, false, false);
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
	RefTax tax(small.path, 7, false, true);
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

void endToEnd() {
	currentGroup = "LCA driver";
	TempFile tax("aa\tk__Bacteria;p__P;c__C;o__O;f__F;g__Genus;s__Genus beta\n"
		"bb\tk__Bacteria;p__P;c__C;o__O;f__F;g__Genus;s__Genus alpha\n"
		"cc\tk__Bacteria;p__P;c__C;o__O;f__F;g__Genus;s__Genus alpha\n"
		"sp\tk__Bacteria;p__P;c__C;o__O;f__F;g__Genus;s__Genus sp. X\n");
	std::unique_ptr<RefTax> reference;
	{
		Silence quiet;
		reference = std::make_unique<RefTax>(tax.path, 7, false, true);
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
	endToEnd();
	matrixCounts();
	emptyConversion();

	if (failures != 0) {
		std::cerr << failures << " regression check(s) failed\n";
		return 1;
	}
	std::cout << "All regression checks passed\n";
	return 0;
}
