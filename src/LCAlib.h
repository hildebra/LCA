// LCAlib.h - LCA's assignment of one query, as a header-only library (since 0.33).
// SPDX-License-Identifier: GPL-2.0-only. Part of LCA, github.com/hildebra/LCA; the license is LCA's LICENSE file.
//
// The LCA program reads alignment and taxonomy files; this header gives the same assignment to programs that hold a
// query's hits in memory, such as protal placing reads on its rRNA reference. It needs only the C++20 standard library.
// A Taxonomy is built once and is then read-only, so threads may share it; each thread passes its own Scratch, and
// Assign allocates nothing once the Scratch has grown to the largest query. For one query, Assign applies:
//   1. hits shorter than min_align_len, or covering less than min_cover of the query, are dropped;
//   2. the anchor is the hit of the highest identity; query coverage, then alignment length, break ties;
//   3. kept are the hits whose identity is within IdentityWindow(anchor's identity) of the anchor's and that cover at
//      least kRelativeCover as much of the query as the anchor (best_hit: the hits ranking equal to the anchor instead);
//   4. several alignments to one subject are one vote, its best alignment's;
//   5. each hit votes down to the deepest rank whose min_identity it reaches; when one hit names a certain species,
//      the placeholder species ("Genus sp. X") of the others abstain;
//   6. from the domain down, a rank is assigned when one name has at least lca_frac of the known votes there (ties:
//      the alphabetically first name); a hit naming another taxon at an assigned rank stops voting, a hit without a
//      name there abstains, and a rank that no hit names stays unknown while deeper ranks are still decided;
//   7. the mean identity of the kept hits limits the depth as a single hit's identity would.
// no_filters skips 1 and 3 (best_hit's choice still applies). LCAimpl.cpp applies the same rules to the program's
// records, and tests/regression.cpp checks that both give the same assignments. The naming rules for taxonomy names
// (parseLineage) are the program's own.
// A Taxonomy can be saved (Parts, since 0.34) and rebuilt without interning every subject's names again (FromParts);
// saved parts are valid only under the kNamingRulesVersion they were made with.
#pragma once

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <deque>
#include <iterator>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace lca {
inline constexpr int kDefaultDepth = 7; // domain, phylum, class, order, family, genus, species
inline constexpr int kMaxDepth = 64;    // as -tdep
inline constexpr int kSpeciesRank = 6;
inline constexpr double kRelativeCover = 0.85; // rule 3: of the anchor's query coverage
// Raise it whenever a subject's names could come out differently for the same input: is_unknown_taxon_name's list,
// species_name_is_placeholder, storeRankName (species shortened to two words, Candidatus to three), Add's trimming,
// parseLineage's rank prefixes. Saved TaxonomyParts of another version must be rebuilt from the lineage text.
inline constexpr int kNamingRulesVersion = 1;

//*******************************************************
//        Taxonomy names
//*******************************************************

inline bool iequals_ascii(std::string_view value, const char* txt) {
	size_t i = 0;
	for (; txt[i] != '\0'; ++i) {
		if (i >= value.size()) { return false; }
		if (std::tolower(static_cast<unsigned char>(value[i])) != std::tolower(static_cast<unsigned char>(txt[i]))) {
			return false;
		}
	}
	return i == value.size();
}

inline bool istarts_with_ascii(std::string_view value, const char* txt) {
	for (size_t i = 0; txt[i] != '\0'; ++i) {
		if (i >= value.size()) { return false; }
		if (std::tolower(static_cast<unsigned char>(value[i])) != std::tolower(static_cast<unsigned char>(txt[i]))) {
			return false;
		}
	}
	return true;
}

inline std::string_view trim_view(std::string_view value, const char* whitespace) {
	const size_t begin = value.find_first_not_of(whitespace);
	if (begin == std::string_view::npos) { return std::string_view(); }
	return value.substr(begin, value.find_last_not_of(whitespace) - begin + 1);
}

inline int tax_rank_from_prefix(char rank) {
	switch (std::tolower(static_cast<unsigned char>(rank))) {
	case 'd':
	case 'k': return 0;
	case 'p': return 1;
	case 'c': return 2;
	case 'o': return 3;
	case 'f': return 4;
	case 'g': return 5;
	case 's': return 6;
	case 't': return 7;
	default: return -1;
	}
}

inline bool is_name_separator(char c) {
	return c == '_' || std::isspace(static_cast<unsigned char>(c));
}

// Species names that do not identify a species, e.g. "Genus sp. X", "Lachnospiraceae bacterium 3_1",
// "gut metagenome". Markers must be whole words; spaces and underscores separate words.
inline bool species_name_is_placeholder(std::string_view value) {
	static const char* const markers[] = {
		"sp", "sp.", "spp", "spp.", "cf", "cf.", "aff", "aff.", "bacterium", "archaeon", "metagenome",
		"unidentified", "unknown", "environmental", "uncultured", "unclassified", "incertae"
	};
	if (istarts_with_ascii(value, "uncultured") || istarts_with_ascii(value, "unclassified")) {
		return true;
	}
	size_t begin = 0;
	while (begin < value.size()) {
		while (begin < value.size() && is_name_separator(value[begin])) { ++begin; }
		size_t end = begin;
		while (end < value.size() && !is_name_separator(value[end])) { ++end; }
		if (end > begin) {
			const std::string_view word = value.substr(begin, end - begin);
			for (const char* marker : markers) {
				if (iequals_ascii(word, marker)) { return true; }
			}
		}
		begin = end;
	}
	return false;
}

inline bool is_unknown_taxon_name(std::string_view token) {
	static const char* const unknownNames[] = {
		"unclassified", "uncultured", "uncultured bacterium", "unknown", "unidentified", "unassigned"
	};
	if (token.empty() || token == "?") { return true; }
	for (const char* name : unknownNames) {
		if (iequals_ascii(token, name)) { return true; }
	}
	return false;
}

// Rank names of one lineage as views into its text; an empty view is an unknown rank.
struct LineageViews {
	std::vector<std::string_view> ranks;
	int depth = 0;
	bool speciesUncertain = false;
};

// Stores the trimmed name `token` at `rank` (below the configured depth) of `out`. Unknown names are left out unless
// doNotCheckTax; a species name is shortened to genus and species (a Candidatus name keeps a third word) unless it
// is a placeholder, which marks the lineage's species uncertain.
inline void storeRankName(std::string_view token, int rank, bool doNotCheckTax, LineageViews& out) {
	const bool taxKnown = !is_unknown_taxon_name(token);
	if (rank == kSpeciesRank && !token.empty()) {
		const bool currentSpeciesUncertain = species_name_is_placeholder(token);
		out.speciesUncertain = currentSpeciesUncertain;
		if (taxKnown && !currentSpeciesUncertain) {
			const size_t firstSpace = token.find(' ');
			if (firstSpace != std::string_view::npos) {
				size_t truncateAt = token.find(' ', firstSpace + 1);
				if (istarts_with_ascii(token, "candidatus ") && truncateAt != std::string_view::npos) {
					truncateAt = token.find(' ', truncateAt + 1);
				}
				if (truncateAt != std::string_view::npos) { token = token.substr(0, truncateAt); }
			}
		}
	}
	if (doNotCheckTax || taxKnown) {
		out.ranks[rank] = token;
		out.depth = std::max(out.depth, rank + 1);
	}
}

// A semicolon-separated lineage of a taxonomy file: rank prefixes (d__ or k__, p__, ..., s__, t__) place names at
// their rank, names without one follow the previous rank; names at rank `d` or deeper are ignored.
inline void parseLineage(std::string_view X, int d, bool doNotCheckTax, LineageViews& out) {
	out.ranks.assign(static_cast<size_t>(d), std::string_view());
	out.depth = 0;
	out.speciesUncertain = false;
	int sequentialRank = 0;
	size_t tokenBegin = 0;
	while (tokenBegin <= X.size()) {
		const size_t separator = X.find(';', tokenBegin);
		const size_t tokenEnd = separator == std::string_view::npos ? X.size() : separator;
		std::string_view token = trim_view(X.substr(tokenBegin, tokenEnd - tokenBegin), " \t\r\n");

		int rank = sequentialRank;
		if (token.size() >= 3 && token[1] == '_' && token[2] == '_') {
			const int prefixedRank = tax_rank_from_prefix(token[0]);
			if (prefixedRank >= 0) {
				rank = prefixedRank;
				token = trim_view(token.substr(3), " \t\r\n");
			}
		}
		sequentialRank = std::max(sequentialRank + 1, rank + 1);
		if (rank >= 0 && rank < d) { storeRankName(token, rank, doNotCheckTax, out); }

		if (separator == std::string_view::npos) { break; }
		tokenBegin = separator + 1;
	}
}

//*******************************************************
//        Assignment
//*******************************************************

// A Taxonomy's contents, to save beside kNamingRulesVersion and rebuild with Taxonomy::FromParts.
struct TaxonomyParts {
	int depth = kDefaultDepth;
	std::vector<std::string> names;         // by id; names[0] is "?" (unknown)
	std::vector<uint32_t> ids;              // depth per subject, subject-major
	std::vector<uint8_t> species_uncertain; // per subject, 0 or 1
};

// Reference lineages by subject index, names interned (id 0: unknown, named "?"), under the program's naming rules.
class Taxonomy
{
public:
	explicit Taxonomy(int depth = kDefaultDepth) : m_depth(std::clamp(depth, 1, kMaxDepth)) { m_names.emplace_back("?"); }
	Taxonomy(const Taxonomy&) = delete;            // the name index refers into this object's names
	Taxonomy& operator=(const Taxonomy&) = delete;
	Taxonomy(Taxonomy&&) = default;                // a deque's elements stay in place, and so the index stays valid
	Taxonomy& operator=(Taxonomy&&) = default;

	// Adds the next subject: names per rank from the domain down, without rank prefixes; ranks after them are unknown.
	// Returns the subject's index.
	uint32_t Add(std::span<const std::string_view> names) {
		m_views.ranks.assign(static_cast<size_t>(m_depth), std::string_view());
		m_views.depth = 0;
		m_views.speciesUncertain = false;
		for (size_t rank = 0; rank < names.size() && rank < static_cast<size_t>(m_depth); ++rank) {
			storeRankName(trim_view(names[rank], " \t\r\n"), static_cast<int>(rank), false, m_views);
		}
		return Store();
	}
	// Adds the next subject from a lineage as taxonomy files have it, e.g. "d__Bacteria;p__P;...;s__Genus alpha".
	uint32_t AddLineage(std::string_view lineage) {
		parseLineage(lineage, m_depth, false, m_views);
		return Store();
	}

	size_t Size() const { return m_species_uncertain.size(); }
	int Depth() const { return m_depth; }
	uint32_t NameId(uint32_t subject, int rank) const {
		return rank >= 0 && rank < m_depth ? m_ids[static_cast<size_t>(subject) * m_depth + rank] : 0;
	}
	const std::string& Name(uint32_t id) const { return m_names[id]; }
	size_t DistinctNames() const { return m_names.size() - 1; }
	// The subject's species name is a placeholder such as "Genus sp. X".
	bool SpeciesUncertain(uint32_t subject) const { return m_species_uncertain[subject] != 0; }

	// A copy of the contents, to save beside kNamingRulesVersion.
	TaxonomyParts Parts() const {
		return {m_depth, std::vector<std::string>(m_names.begin(), m_names.end()), m_ids, m_species_uncertain};
	}
	// Rebuilds a Taxonomy from Parts made under this kNamingRulesVersion, moving the vectors in. The names are taken as
	// they are, and only the name index is rebuilt (one insert per name), so that Add reuses their ids. Inconsistent
	// parts give nullopt and the reason in `error`, and are left unchanged.
	static std::optional<Taxonomy> FromParts(TaxonomyParts&& parts, std::string& error) {
		const auto refuse = [&error](std::string reason) {
			error = std::move(reason);
			return std::optional<Taxonomy>();
		};
		if (parts.depth < 1 || parts.depth > kMaxDepth) {
			return refuse("depth " + std::to_string(parts.depth) + " is outside 1 to " + std::to_string(kMaxDepth));
		}
		if (parts.names.empty()) { return refuse("no names; name 0 must be \"?\""); }
		if (parts.names[0] != "?") { return refuse("name 0 is \"" + parts.names[0] + "\", not \"?\""); }
		for (size_t id = 1; id < parts.names.size(); ++id) {
			if (parts.names[id].empty()) { return refuse("name " + std::to_string(id) + " is empty"); }
		}
		const size_t depth = static_cast<size_t>(parts.depth), subjects = parts.species_uncertain.size();
		if (parts.ids.size() != depth * subjects) {
			return refuse(std::to_string(parts.ids.size()) + " ids for " + std::to_string(subjects) + " subjects at depth " +
				std::to_string(depth) + ", not " + std::to_string(depth * subjects));
		}
		const size_t names = parts.names.size();
		const auto past = std::find_if(parts.ids.begin(), parts.ids.end(), [names](uint32_t id) { return id >= names; });
		if (past != parts.ids.end()) {
			const size_t i = static_cast<size_t>(past - parts.ids.begin());
			return refuse("id " + std::to_string(*past) + " at index " + std::to_string(i) + " (subject " +
				std::to_string(i / depth) + ", rank " + std::to_string(i % depth) + ") is past the " +
				std::to_string(names) + " names");
		}

		Taxonomy tax(parts.depth);
		tax.m_names.assign(std::make_move_iterator(parts.names.begin()), std::make_move_iterator(parts.names.end()));
		tax.m_name_ids.reserve(names);
		for (size_t id = 1; id < names; ++id) {
			const std::string& name = tax.m_names[id];
			const auto [entry, added] = tax.m_name_ids.emplace(std::string_view(name), static_cast<uint32_t>(id));
			if (!added || name == tax.m_names[0]) { // a name's votes would be split between its ids
				std::string reason = "name " + std::to_string(id) + " (\"" + name + "\") repeats name " +
					std::to_string(added ? 0 : entry->second);
				std::move(tax.m_names.begin(), tax.m_names.end(), parts.names.begin());
				return refuse(std::move(reason));
			}
		}
		tax.m_ids = std::move(parts.ids);
		tax.m_species_uncertain = std::move(parts.species_uncertain);
		return tax;
	}

private:
	uint32_t Store() {
		for (int rank = 0; rank < m_depth; ++rank) {
			m_ids.push_back(rank < m_views.depth ? Intern(m_views.ranks[rank]) : 0);
		}
		m_species_uncertain.push_back(m_views.speciesUncertain);
		return static_cast<uint32_t>(m_species_uncertain.size() - 1);
	}
	uint32_t Intern(std::string_view name) {
		if (name.empty()) { return 0; }
		const auto found = m_name_ids.find(name);
		if (found != m_name_ids.end()) { return found->second; }
		m_names.emplace_back(name);
		const uint32_t id = static_cast<uint32_t>(m_names.size() - 1);
		m_name_ids.emplace(std::string_view(m_names.back()), id); // deque elements do not move
		return id;
	}

	int m_depth;
	std::deque<std::string> m_names;                     // by id
	std::unordered_map<std::string_view, uint32_t> m_name_ids;
	std::vector<uint32_t> m_ids;                         // m_depth per subject
	std::vector<uint8_t> m_species_uncertain;            // per subject
	LineageViews m_views;                                // of the lineage being added
};

struct Hit {
	uint32_t subject = 0;       // index into the Taxonomy; hits of subjects outside it are ignored
	double identity = 0;        // percent
	uint32_t align_len = 0;     // alignment length
	uint32_t query_len = 0;     // the query's length; 0: unknown, which skips the coverage filter (rule 1)
	uint32_t query_covered = 0; // query bases the alignment covers

	uint32_t Covered() const { return std::min(query_covered, query_len); }
	double Cover() const { return query_len != 0 ? static_cast<double>(Covered()) / query_len : 1.0; }
};

// The program's options of the same names in brackets.
struct Options {
	// Minimum identity per rank in percent, inclusive, from the deepest rank up to the domain, one per rank (-id).
	std::vector<double> min_identity = {97, 95, 93, 91, 88, 78, 0};
	double lca_frac = 0.9;       // -LCAfrac
	uint32_t min_align_len = 75; // -minAlignLen
	double min_cover = 0.5;      // -cover
	bool best_hit = false;       // -reportBestHit's choice of hits; the program also sets every min_identity to 1
	bool no_filters = false;     // -no_bl_filter: rules 1 and 3 are skipped, except best_hit's choice
	int depth = kDefaultDepth;   // -tdep; the Taxonomy's depth if that is smaller

	// The cutoff of a rank, 0 the domain; 0 for a rank without one.
	double CutoffOf(int rank) const {
		const size_t i = static_cast<size_t>(depth - 1 - rank);
		return i < min_identity.size() ? min_identity[i] : 0.0;
	}
};

struct Assignment {
	int depth = 0;                              // ranks assigned from the domain down (0: none)
	std::array<uint32_t, kMaxDepth> name_ids{}; // per rank; 0 past depth, or at a rank inside it that no hit names
	double mean_identity = 0;                   // of the kept hits
	uint32_t hits_used = 0;                     // kept subjects (0: no hit passed the filters)
	bool species_uncertain = false;             // the assigned species is a placeholder name
};

// Per-thread buffers that Assign reuses.
struct Scratch {
	struct Voter {
		uint32_t subject;
		int eligible;       // ranks the hit's identity reaches
		bool speciesMasked; // a placeholder species beside a certain one
		bool active;        // no conflict with an assigned rank
	};
	struct Vote {
		uint32_t id;
		uint32_t count;
		size_t voter; // the first naming it
	};
	std::vector<Hit> kept;
	std::vector<Voter> voters;
	std::vector<Vote> votes;
};

// Rule 3's window below the anchor's identity, in percentage points. The float literals are the program's values
// since 0.31 (as doubles, 0.05f and 0.15f are slightly larger than 0.05 and 0.15).
inline double IdentityWindow(double anchorIdentity) {
	if (anchorIdentity >= 100) { return 0.05f; }
	if (anchorIdentity >= 99.5) { return 0.15f; }
	if (anchorIdentity >= 99) { return 0.25f; }
	if (anchorIdentity >= 98) { return 0.75f; }
	if (anchorIdentity >= 97) { return 1.0f; }
	return 1.5;
}

// Negative when a ranks above b: identity, then query coverage (a known coverage first), then alignment length.
inline int CompareRank(const Hit& a, const Hit& b) {
	if (a.identity != b.identity) { return a.identity > b.identity ? -1 : 1; }
	const bool aKnown = a.query_len != 0, bKnown = b.query_len != 0;
	if (aKnown != bKnown) { return aKnown ? -1 : 1; }
	if (aKnown && a.Cover() != b.Cover()) { return a.Cover() > b.Cover() ? -1 : 1; }
	if (a.align_len != b.align_len) { return a.align_len > b.align_len ? -1 : 1; }
	return 0;
}

// Rule 3's coverage: of the anchor's covered query bases, or of its alignment length if a query length is unknown.
inline bool CoversEnough(const Hit& hit, const Hit& anchor) {
	if (hit.query_len != 0 && anchor.query_len != 0) {
		if (hit.query_len == anchor.query_len) {
			return hit.Covered() >= std::ceil(static_cast<double>(anchor.Covered()) * kRelativeCover);
		}
		return hit.Cover() >= anchor.Cover() * kRelativeCover;
	}
	return hit.align_len >= std::ceil(static_cast<double>(anchor.align_len) * kRelativeCover);
}

inline Assignment Assign(std::span<const Hit> hits, const Taxonomy& tax, const Options& opt, Scratch& scratch) {
	Assignment result;
	const int depth = std::clamp(std::min(opt.depth, tax.Depth()), 0, kMaxDepth);
	std::vector<Hit>& kept = scratch.kept;
	kept.clear();
	for (const Hit& hit : hits) {
		if (hit.subject >= tax.Size()) { continue; }
		if (!opt.no_filters &&
			(hit.align_len < opt.min_align_len || (hit.query_len != 0 && hit.Cover() < opt.min_cover))) { continue; }
		kept.push_back(hit);
	}
	if (kept.empty()) { return result; }

	Hit anchor = kept.front();
	for (const Hit& hit : kept) {
		if (CompareRank(hit, anchor) < 0) { anchor = hit; }
	}
	if (opt.best_hit) {
		std::erase_if(kept, [&anchor](const Hit& hit) { return CompareRank(hit, anchor) != 0; });
	} else if (!opt.no_filters) {
		const double window = IdentityWindow(anchor.identity);
		std::erase_if(kept, [&anchor, window](const Hit& hit) {
			return hit.identity + window < anchor.identity || !CoversEnough(hit, anchor);
		});
	}
	// One hit per subject, where the subject first occurs; queries have few hits, so a linear search will do.
	size_t subjects = 0;
	for (size_t i = 0; i < kept.size(); ++i) {
		const Hit hit = kept[i];
		size_t j = 0;
		while (j < subjects && kept[j].subject != hit.subject) { ++j; }
		if (j == subjects) { kept[subjects++] = hit; }
		else if (CompareRank(hit, kept[j]) < 0) { kept[j] = hit; }
	}
	kept.resize(subjects);

	long double identitySum = 0.0L; // summed as the program does, so the means are identical
	for (const Hit& hit : kept) { identitySum += hit.identity; }
	result.hits_used = static_cast<uint32_t>(kept.size());
	result.mean_identity = static_cast<double>(identitySum / kept.size());
	auto eligibleDepth = [&opt, depth](double identity) {
		int eligible = 0;
		while (eligible < depth && opt.CutoffOf(eligible) <= identity) { ++eligible; }
		return eligible;
	};

	std::vector<Scratch::Voter>& voters = scratch.voters;
	voters.clear();
	bool certainSpecies = false;
	for (const Hit& hit : kept) {
		const int eligible = eligibleDepth(hit.identity);
		voters.push_back({hit.subject, eligible, false, true});
		if (eligible > kSpeciesRank && tax.NameId(hit.subject, kSpeciesRank) != 0 && !tax.SpeciesUncertain(hit.subject)) {
			certainSpecies = true;
		}
	}
	if (certainSpecies && voters.size() > 1) {
		for (Scratch::Voter& voter : voters) {
			voter.speciesMasked = voter.eligible > kSpeciesRank && tax.SpeciesUncertain(voter.subject);
		}
	}
	auto voteOf = [&tax](const Scratch::Voter& voter, int rank) -> uint32_t {
		if (!voter.active || rank >= voter.eligible || (voter.speciesMasked && rank == kSpeciesRank)) { return 0; }
		return tax.NameId(voter.subject, rank);
	};

	size_t speciesVoter = 0; // a hit naming the assigned species
	if (voters.size() == 1) {
		for (int rank = 0; rank < voters.front().eligible; ++rank) { result.name_ids[rank] = voteOf(voters.front(), rank); }
		result.depth = voters.front().eligible;
	} else {
		std::vector<Scratch::Vote>& votes = scratch.votes;
		for (int rank = 0; rank < depth; ++rank) {
			votes.clear();
			uint32_t known = 0;
			for (size_t i = 0; i < voters.size(); ++i) {
				const uint32_t id = voteOf(voters[i], rank);
				if (id == 0) { continue; }
				++known;
				// Distinct names at a rank are few, so a linear scan beats hashing.
				auto vote = std::find_if(votes.begin(), votes.end(), [id](const Scratch::Vote& v) { return v.id == id; });
				if (vote == votes.end()) { votes.push_back({id, 1, i}); }
				else { ++vote->count; }
			}
			if (known == 0) { continue; }

			const Scratch::Vote* best = &votes.front();
			for (const Scratch::Vote& vote : votes) {
				if (vote.count > best->count || (vote.count == best->count && tax.Name(vote.id) < tax.Name(best->id))) {
					best = &vote;
				}
			}
			if (best->count < static_cast<double>(known) * opt.lca_frac) { break; }
			const uint32_t consensus = best->id;
			result.name_ids[rank] = consensus;
			result.depth = rank + 1;
			if (rank == kSpeciesRank) { speciesVoter = best->voter; }
			for (Scratch::Voter& voter : voters) {
				const uint32_t id = voteOf(voter, rank);
				if (id != 0 && id != consensus) { voter.active = false; }
			}
		}
	}

	result.depth = std::min(result.depth, eligibleDepth(result.mean_identity));
	while (result.depth > 0 && result.name_ids[result.depth - 1] == 0) { --result.depth; }
	std::fill(result.name_ids.begin() + result.depth, result.name_ids.begin() + depth, 0u);
	result.species_uncertain = result.depth > kSpeciesRank && result.name_ids[kSpeciesRank] != 0 &&
		tax.SpeciesUncertain(voters[speciesVoter].subject);
	return result;
}
}
