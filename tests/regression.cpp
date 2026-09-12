#include "../LCAimpl.h"

#include <cmath>
#include <iostream>

namespace {
int failures = 0;

void check(bool condition, const char* message) {
	if (!condition) {
		std::cerr << "FAIL: " << message << '\n';
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

options defaultOptions() {
	char a0[] = "regression";
	char a1[] = "-i";
	char a2[] = "input.m8";
	char a3[] = "-r";
	char a4[] = "reference.tax";
	char a5[] = "-o";
	char a6[] = "output.tsv";
	char* argv[] = {a0, a1, a2, a3, a4, a5, a6};
	return options(7, argv, 7);
}
}

int main() {
	TaxObj prefixed("s__Genus sp.;k__Bacteria;g__Genus", 7, false, false);
	check(prefixed.get(0) == "Bacteria", "k__ must map to Domain even after s__");
	check(prefixed.get(5) == "Genus", "g__ must map to Genus even after s__");
	check(prefixed.get(6) == "Genus sp.", "s__ must map to Species");
	check(prefixed.speciesUncertain, "a species containing sp. must be marked uncertain");
	TaxObj candidatus("k__Bacteria;p__P;c__C;o__O;f__F;g__Liberibacter;s__Candidatus Liberibacter asiaticus strain X",
		7, false, false);
	check(candidatus.get(6) == "Candidatus Liberibacter asiaticus",
		"Candidatus species names must retain their three-part species name");

	BlastRes crlf("q\ts\t99\t100\t0\t0\t100\t1\t200\t101\t100\r", 0);
	check(!crlf.fail && crlf.Query == "q" && crlf.Sbj == "s", "CRLF m8 records must parse");
	check(std::fabs(crlf.Qcoverage - 1.0f) < 1e-6f, "reverse query coordinates must produce valid coverage");
	check(BlastRes::isColumnHeader(
		"qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tqlen\r\n"),
		"the documented custom m8 header must be recognized");
	check(BlastRes::isColumnHeader(
		"qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"),
		"the legacy 12-column BLAST header must be recognized");
	BlastRes spaced("q s 99 100 0 0 1 100 1 100 100", 0);
	check(spaced.fail,
		"ordinary spaces must not be interpreted as column separators");
	BlastRes legacy(
		"HSQ-700358:311:H23FFBCXX:1:2210:10499:32978\tAACY020299484\t96.70\t91\t3\t0\t203\t293\t1912\t2002\t1e-35\t 152",
		0);
	check(!legacy.fail && !legacy.queryCoverageKnown && legacy.eval == 1e-35 && legacy.score == 152.0,
		"legacy qseqid-to-bitscore 12-column BLAST records must parse");
	BlastRes extra("q\ts\t99\t100\t0\t0\t1\t100\t1\t100\t100\t1e-20\t50", 0);
	check(extra.fail, "custom m8 records with an invalid column count must be rejected");

	BlastRes gapped("q\ts\t99\t100\t0\t1\t1\t90\t1\t100\t100", 0);
	check(!gapped.fail, "valid gapped custom m8 record must parse");
	check(std::fabs(gapped.Qcoverage - 0.9f) < 1e-6f,
		"coverage must use query-coordinate span, not gapped alignment length");
	BlastRes malformed("q\ts\t99\t100\tbad\t0\t1\t90\t1\t100\t100", 0);
	check(malformed.fail, "numeric custom m8 columns must be validated");

	options opt = defaultOptions();
	opt.minCover = 1.0f;
	std::vector<BlastRes> legacyHits = {legacy};
	double legacyBest = 0.0;
	filterBlastPrimary(legacyHits, &opt, legacyBest);
	check(legacyHits.size() == 1,
		"coverage filtering must not reject legacy records whose query length is unavailable");
	opt.minCover = 0.5f;
	std::vector<BlastRes> forward = {
		hit("q", "high_identity", 100.0, 100, 100, 200),
		hit("q", "longer", 95.0, 120, 120, 200)
	};
	std::vector<BlastRes> reverse = {forward[1], forward[0]};
	double forwardBest = 0.0, reverseBest = 0.0;
	filterBlastPrimary(forward, &opt, forwardBest);
	filterBlastPrimary(reverse, &opt, reverseBest);
	check(forward.size() == 1 && reverse.size() == 1, "primary filter should retain one anchor hit in this case");
	check(forward.front().Sbj == "longer" && reverse.front().Sbj == "longer",
		"primary filtering must be independent of input hit order");
	check(forwardBest == reverseBest, "order must not change the anchor identity");

	std::vector<TaxObj*> parentTest;
	for (int i = 0; i < 8; ++i) {
		TaxObj* tax = new TaxObj(2);
		tax->set(0, "A");
		parentTest.push_back(tax);
	}
	for (int i = 0; i < 2; ++i) {
		TaxObj* tax = new TaxObj(2);
		tax->set(1, i == 0 ? "invalid-child-a" : "invalid-child-b");
		parentTest.push_back(tax);
	}
	bool singular = false;
	TaxObj* lineage = LCAcore(parentTest, singular, 0.9, 2);
	check(lineage->depth == 1 && lineage->get(0) == "A" && lineage->get(1) == __unkwnTax,
		"children without the accepted parent must not influence deeper LCA ranks");
	delete lineage;
	for (TaxObj* tax : parentTest) { delete tax; }

	TaxObj output(7);
	for (int i = 0; i < 7; ++i) { output.set(i, "rank" + std::to_string(i)); }
	output.setHitDB("subject-read");
	output.setRepID(true);
	output.perID = 98.5f;
	const std::string rendered = output.getWriteString(std::vector<double>(7, 0.0));
	check(rendered == "rank0\trank1\trank2\trank3\trank4\trank5\trank6\tsubject-read\t98.500000",
		"hit read and identity must follow all normal taxonomy columns");
	TaxObj tieA(&output);
	TaxObj tieB(&output);
	tieA.setHitDB("z-subject");
	tieB.setHitDB("a-subject");
	check(tieA.evalAcpyTax(&tieB) && tieA.getHitDB() == "a-subject",
		"equal multi-database assignments must select hit reads deterministically");

	char b0[] = "regression";
	char b1[] = "-i";
	char b2[] = "input.m8";
	char b3[] = "-r";
	char b4[] = "reference.tax";
	char b5[] = "-o";
	char b6[] = "output.tsv";
	char b7[] = "-tdep";
	char b8[] = "8";
	char* depthArgv[] = {b0, b1, b2, b3, b4, b5, b6, b7, b8};
	options depthOptions(9, depthArgv, 7);
	check(depthOptions.taxDepth == 8 && depthOptions.idThr.size() == 8 &&
		depthOptions.Taxlvls.size() == 8 && depthOptions.Taxlvls[7] == "Strain",
		"-tdep must resize taxonomy thresholds and labels consistently");

	std::string queryToken = "stale";
	check(!BlastRes::extractQueryToken(" \t\r\n", queryToken),
		"blank lines must not contain a query token");
	check(BlastRes::extractQueryToken("  query with barcode \tsubject", queryToken) &&
		queryToken == "  query with barcode ", "query extraction must preserve the complete first tab-separated field");
	BlastRes padded("  q description \t s description \t 99 \t 100 \t0\t0\t1\t100\t1\t100\t 100 \r\n", 0);
	check(!padded.fail && padded.Query == "  q description " && padded.Sbj == " s description ",
		"numeric padding and CRLF must be accepted without trimming identifiers");

	const std::string illumina = "LH00409:413:22WTJGLT4:6:2481:9713:19847 1:N:0:CGTATCTC+CTCGAACA";
	const std::string reportedLine = illumina + "\tOBEP010038528\t 95.3 \t150\t7\t0\t1\t150\t795\t944\t150";
	BlastRes reported(reportedLine, 0);
	check(!reported.fail && reported.Query == illumina && reported.Sbj == "OBEP010038528" &&
		reported.perID == 95.3 && reported.queryCoverageKnown && reported.Qcoverage == 1.0f,
		"the reported Usearch row must preserve its complete Illumina query header");
	check(BlastRes::supportedColumnCount(reportedLine) == 11,
		"an embedded identifier space must not change the detected layout");
	check(BlastRes::extractQueryToken(reportedLine, queryToken) && queryToken == illumina,
		"query extraction must agree with parsing for Illumina metadata");
	check(!BlastRes::isColumnHeader("qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen"),
		"headers must use tabs as column separators");

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
	check(BlastRes("q\ts\t9 9\t100\t0\t0\t1\t100\t1\t100\t100", 0).fail,
		"spaces within a number must not be accepted as numeric padding");
	BlastRes zeroCoordinate("q\ts\t99\t100\t0\t0\t0\t99\t1\t100\t100", 0);
	check(zeroCoordinate.fail, "BLAST coordinates must be positive and one-based");
	BlastRes extraFields("q\ts\t99\t100\t0\t0\t1\t100\t1\t100\t100\t1\t2\t3", 0);
	check(extraFields.fail, "extra fields must not be silently discarded");

	std::vector<BlastRes> boundary = {hit("q", "a", 97, 100, 100, 100), hit("q", "b", 97, 85, 85, 100)};
	double bestID = 0.0;
	filterBlastPrimary(boundary, &opt, bestID);
	check(boundary.size() == 2, "an alignment at the exact 85% length cutoff must be retained");

	std::vector<BlastRes> duplicates = {hit("q", "a", 100, 50, 50, 100), hit("q", "a", 99, 100, 100, 100)};
	filterBlastPrimary(duplicates, &opt, bestID);
	check(duplicates.size() == 1 && duplicates[0].perID == 99,
		"a short duplicate must not hide an eligible alignment for the same subject");
	duplicates = {hit("q", "a", 100, 100, 100, 300), hit("q", "a", 99, 100, 100, 100)};
	filterBlastPrimary(duplicates, &opt, bestID);
	check(duplicates.size() == 1 && duplicates[0].perID == 99,
		"a low-coverage duplicate must not hide an eligible alignment for the same subject");

	opt.reportBestHit = true;
	for (bool filtered : {true, false}) {
		opt.BLfilter = filtered;
		std::vector<BlastRes> tied = {hit("q", "b", 99, 100, 100, 100), hit("q", "a", 99, 100, 100, 100)};
		filterBlastPrimary(tied, &opt, bestID);
		check(tied.size() == 1 && tied[0].Sbj == "a",
			"best-hit mode must select one deterministic subject even without primary filtering");
	}

	TaxObj truncated(7);
	for (int i = 0; i < 7; ++i) { truncated.set(i, "rank" + std::to_string(i)); }
	truncated.limitByIdentity(94, defaultOptions().idThr);
	TaxObj moreKnown(7);
	for (int i = 0; i < 6; ++i) { moreKnown.set(i, "rank" + std::to_string(i)); }
	check(truncated.depth == 5 && truncated.evalAcpyTax(&moreKnown) && truncated.depth == 6,
		"taxonomy beyond the identity cutoff must not count toward database selection");

	std::vector<BlastRes> empty;
	double meanIdentity = 123.0;
	const auto noTax = BlastToTax(empty, nullptr, &opt, meanIdentity);
	check(noTax.empty() && meanIdentity == 0.0f,
		"empty hit conversion must reset the mean without dividing by zero");

	if (failures != 0) {
		std::cerr << failures << " regression check(s) failed\n";
		return 1;
	}
	std::cout << "All regression checks passed\n";
	return 0;
}
