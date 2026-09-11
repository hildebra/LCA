#include "options.h"
#include <cerrno>
#include <climits>

namespace {
bool parseStrictDouble(const char* text, double& value) {
	if (text == nullptr || *text == '\0') { return false; }
	char* end = nullptr;
	errno = 0;
	value = std::strtod(text, &end);
	return end != text && *end == '\0' && errno != ERANGE && std::isfinite(value);
}

bool parseStrictInt(const char* text, int& value) {
	if (text == nullptr || *text == '\0') { return false; }
	char* end = nullptr;
	errno = 0;
	const long parsed = std::strtol(text, &end, 10);
	if (end == text || *end != '\0' || errno == ERANGE || parsed < INT_MIN || parsed > INT_MAX) { return false; }
	value = static_cast<int>(parsed);
	return true;
}

[[noreturn]] void invalidOptionValue(const char* option, const string& value, const string& expected) {
	cerr << "Invalid value '" << value << "' for " << option << "; expected " << expected << ".\n";
	exit(6);
}
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
	size_t begin = 0;
	while (begin <= s.size()) {
		const size_t end = s.find(delim, begin);
		elems.push_back(s.substr(begin, end == string::npos ? string::npos : end - begin));
		if (end == string::npos) { break; }
		begin = end + 1;
	}
	return elems;
}

void self_help() {
	cout << "LCA help\nUsage: ./LCA [[ optional_args ]] -i [blast m8 output] -r [taxonomy database] -o [output file]\n";
	cout << "Required arguments:\n";
	cout << "  -i custom 11-column m8 ending in qlen, or legacy 12-column BLAST ending in evalue bitscore\n";
	cout << "  -r taxonomy file with entries corresponding to sequences in ref database, that was mapped against\n";
	cout << "  -o output file containing the sequence name and the assigned taxonomy against the ref database\n";
	cout << "Optional arguments:\n";
	cout << "  -matHigh           calculate abundance of reads at different taxonomic levels. An extra file (derriving from -o) per tax level is written\n";
	cout << "  -showHitRead       append the uniquely assigned database entry after the taxonomy (or ? when the LCA used multiple entries)\n";
	cout << "  -no_bl_filter      use only, if custom scripts were used to pre-filter filter -i file and in-built filter should be deactivated\n";
	cout << "  -readInput         Are the inputs miTags? Default: off (assummes OTUs) \n";
	cout << "  -LCAfrac           (0-1] fraction of matching known taxonomies required at each level. Default=\"0.9\"\n";
	cout << "  -t                 retained for compatibility; only -t 1 is supported by this single-threaded build\n";
	cout << "  -tdep              [int] number of taxonomy levels to read and report. Default=7\n";
	cout << "  -id                comma-separated min %identity values from the deepest configured rank back to Domain; supply exactly -tdep values. Default=\"97,95,93,91,88,78,0\"\n";
	cout << "  -cover             [0-1] query coverage required when qlen is present; skipped for legacy 12-column input. Default=0.5\n";
	cout << "  -minAlignLen       [int] min num basepairs of hit to accept reported hit. Default=\"75\"\n";
	cout << "  -SLVfmt            retained for compatibility; rank prefixes are auto-detected. Default: off\n";
	cout << "  -reportID          append the mean identity of retained LCA hits\n";
	cout << "  -reportBestHit     report the deterministic best hit and use its complete eligible taxonomy\n";
	cout << "  -reportHitPattern  [file] Write a report of best search hits to reproduce Grant et al. (2023) hit patterns for LCA search \n";
	cout << endl;
	exit(0);
}


options::options(int argc, char **argv,int defDep):
	RefTaxFile(""), blastres(""), outF(""), input_format("bl8"), repHitPattern(""),
	BLfilter(true), calcHighMats(false), hitRD(false), isReads(false),
	nativeSlVdb(false), reportID(false), reportBestHit(false), checkTaxoUnkw(true),
	numThr(1), taxDepth(defDep), LCAfract(0.9), minCover(0.5f), minAliLen(75), idThr(),
	blFiles(0), refDBs(0), Taxlvls(), version(false)
{
	if (defDep <= 0) {
		cerr << "Configured taxonomy depth must be positive.\n";
		exit(6);
	}
	bool newIDthrs = false; string newIDthStr("");
	auto requireValue = [&](int& idx, const char* optName) -> const char* {
		if (idx + 1 >= argc) {
			cerr << "Missing value for option " << optName << "\n";
			exit(6);
		}
		return argv[++idx];
	};

	for (int i = 1; i < argc; i++)
	{
		if (!strcmp(argv[i], "-i"))
           blastres = requireValue(i, "-i");
		else if (!strcmp(argv[i], "-h"))
			self_help();
		else if (!strcmp(argv[i], "-r"))
         RefTaxFile = requireValue(i, "-r");
		else if (!strcmp(argv[i], "-f"))//input format
           input_format = requireValue(i, "-f");
		else if (!strcmp(argv[i], "-v"))
			version = true;
		else if (!strcmp(argv[i], "-o"))
           outF = requireValue(i, "-o");
		else if (!strcmp(argv[i], "-matHigh"))
			calcHighMats = true;
		else if (!strcmp(argv[i], "-showHitRead"))
			hitRD = true;
		else if (!strcmp(argv[i], "-no_bl_filter"))
			BLfilter = false;
		else if (!strcmp(argv[i], "-no_taxDB_filter"))
			checkTaxoUnkw = false;
		else if (!strcmp(argv[i], "-readInput"))
			isReads = true;
		else if (!strcmp(argv[i], "-minAlignLen"))
		{
			const char* value = requireValue(i, "-minAlignLen");
			int parsed = 0;
			if (!parseStrictInt(value, parsed) || parsed < 0) {
				invalidOptionValue("-minAlignLen", value, "an integer from 0 to INT_MAX");
			}
			minAliLen = parsed;
		}
		else if (!strcmp(argv[i], "-reportID"))
			reportID = true;
		else if (!strcmp(argv[i], "-reportBestHit"))
			reportBestHit = true;
		else if (!strcmp(argv[i], "-SLVfmt"))
			nativeSlVdb = true;
		else if (!strcmp(argv[i], "-reportHitPattern"))
			repHitPattern = requireValue(i, "-reportHitPattern");
		else if (!strcmp(argv[i], "-t"))
		{
			const char* value = requireValue(i, "-t");
			if (!parseStrictInt(value, numThr) || numThr != 1) {
				invalidOptionValue("-t", value, "1 (threading is not enabled in this build)");
			}
		}
		else if (!strcmp(argv[i], "-tdep"))
		{
			const char* value = requireValue(i, "-tdep");
			if (!parseStrictInt(value, taxDepth) || taxDepth <= 0 || taxDepth > 64) {
				invalidOptionValue("-tdep", value, "an integer from 1 to 64");
			}
		}
		else if (!strcmp(argv[i], "-cover"))
		{
			const char* value = requireValue(i, "-cover");
			double parsed = 0.0;
			if (!parseStrictDouble(value, parsed) || parsed < 0.0 || parsed > 1.0) {
				invalidOptionValue("-cover", value, "a finite number from 0 to 1");
			}
			minCover = static_cast<float>(parsed);
		}
		else if (!strcmp(argv[i], "-LCAfrac"))
		{
			const char* value = requireValue(i, "-LCAfrac");
			if (!parseStrictDouble(value, LCAfract) || LCAfract <= 0.0 || LCAfract > 1.0) {
				invalidOptionValue("-LCAfrac", value, "a finite number greater than 0 and at most 1");
			}
		}
		else if (!strcmp(argv[i], "-id")) {
			newIDthrs = true; newIDthStr = requireValue(i, "-id");
		}
		else {
			cerr << "Unknown option: " << argv[i] << "\nUse \"./LCA -h\" to get full help.\n";
			exit(6);
		}
	}
	if (version) { return; }

	const double defaultThresholds[] = {0, 78, 88, 91, 93, 95, 97};
	idThr.resize(taxDepth);
	for (int i = 0; i < taxDepth; ++i) {
		idThr[i] = defaultThresholds[std::min(i, 6)];
	}
	Taxlvls.resize(taxDepth);

	if (input_format != "bl8") {
		cerr << "Unsupported input format '" << input_format
			<< "'. Only the custom bl8 format is supported.\n";
		exit(6);
	}
	split(blastres, ',', blFiles);
	split(RefTaxFile, ',', refDBs);

	if (blFiles.size() != refDBs.size()) {
		cerr << "Unequal number of blast and refDB files!\n"; exit(24);
	}

	//check that all args were given
	bool isErr(false);
	if (blastres == "") { cerr << "Input file invalid or missing (-i)\n"; isErr = true; }
	if (RefTaxFile == "") { cerr << "RefDb file invalid or missing (-r)\n"; isErr = true; }
	if (outF == "") { cerr << "Output file invalid or missing (-o)\n"; isErr = true; }
	for (const auto& path : blFiles) {
		if (path.empty()) { cerr << "An empty blast input path was supplied to -i\n"; isErr = true; }
	}
	for (const auto& path : refDBs) {
		if (path.empty()) { cerr << "An empty taxonomy database path was supplied to -r\n"; isErr = true; }
	}
	//if (blastres == "") { cerr << "Input file invalid (-f)"; isErr = true; }
	if (isErr) { cerr << "Use \"./LCA -h\" to get full help.\nError in command line args.. exiting\n"; exit(5); }

	if (newIDthrs) {
		vector<string> idthrsrev;
		split(newIDthStr, ',',idthrsrev);
		if (idthrsrev.size() != static_cast<size_t>(taxDepth)) {
			cerr << "Wrong number of identity threshold levels (needs to be " << taxDepth << ").\nAborting..\n";
			exit(39);
		}
		for (size_t i = 0; i < idthrsrev.size(); i++) {
			double parsed = 0.0;
			const string& text = idthrsrev[idthrsrev.size() - 1 - i];
			if (!parseStrictDouble(text.c_str(), parsed) || parsed < 0.0 || parsed > 100.0) {
				invalidOptionValue("-id", text, "comma-separated finite percentages from 0 to 100");
			}
			idThr[i] = parsed;
		}
		for (size_t i = 1; i < idThr.size(); ++i) {
			if (idThr[i] < idThr[i - 1]) {
				cerr << "Identity thresholds must not decrease from Domain toward deeper ranks.\n";
				exit(39);
			}
		}
	}

	//simply overwrite to low values..
	if (reportBestHit) {
		std::fill(idThr.begin(), idThr.end(), 1.0);
	}


	vector<string> defTLvls(8, "");
	defTLvls[0] = "Domain"; defTLvls[1] = "Phylum"; defTLvls[2] = "Class";
	defTLvls[3] = "Order"; defTLvls[4] = "Family";  defTLvls[5] = "Genus";
	defTLvls[6] = "Species"; defTLvls[7] = "Strain";
	for (size_t i = 0; i < (size_t)taxDepth; i++) {
		Taxlvls[i] = i < defTLvls.size() ? defTLvls[i] : "Level" + to_string(i + 1);
	}

}
const string options::TaxLvl2string() {
	string ret = Taxlvls[0];
	for (size_t i = 1; i < Taxlvls.size(); i++) {
		ret += __defaultTaxSep + Taxlvls[i];
	}
	return ret;
}


options::~options()
{
}
