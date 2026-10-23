#include "LCAimpl.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>

namespace {
// Keeps one alignment per subject: the best one, at the position where the subject first occurs.
void deduplicateBlastHits(vector<BlastRes>& BR) {
	const size_t count = BR.size();
	if (count < 2) { return; }
	thread_local vector<uint32_t> order;
	order.resize(count);
	std::iota(order.begin(), order.end(), 0);
	std::sort(order.begin(), order.end(), [&BR](uint32_t a, uint32_t b) {
		const int bySubject = BR[a].Sbj.compare(BR[b].Sbj);
		return bySubject != 0 ? bySubject < 0 : a < b;
	});
	bool duplicates = false;
	for (size_t i = 1; i < count && !duplicates; ++i) {
		duplicates = BR[order[i]].Sbj == BR[order[i - 1]].Sbj;
	}
	if (!duplicates) { return; }

	thread_local vector<std::pair<uint32_t, uint32_t> > kept; // first occurrence, best alignment
	kept.clear();
	for (size_t i = 0; i < count;) {
		uint32_t best = order[i];
		size_t j = i + 1;
		for (; j < count && BR[order[j]].Sbj == BR[order[i]].Sbj; ++j) {
			if (BR[order[j]].compareRank(BR[best]) < 0) { best = order[j]; }
		}
		kept.emplace_back(order[i], best);
		i = j;
	}
	std::sort(kept.begin(), kept.end());
	vector<BlastRes> unique;
	unique.reserve(kept.size());
	for (const auto& entry : kept) { unique.push_back(std::move(BR[entry.second])); }
	BR.swap(unique);
}

// The fields that rank and filter hits, without copying the identifier strings.
BlastRes rankingFields(const BlastRes& hit) {
	BlastRes key;
	key.alLen = hit.alLen;
	key.perID = hit.perID;
	key.Qcoverage = hit.Qcoverage;
	key.queryCovered = hit.queryCovered;
	key.queryLength = hit.queryLength;
	key.queryCoverageKnown = hit.queryCoverageKnown;
	key.fail = hit.fail;
	return key;
}

// Query coverage relative to the anchor; alignment length only when qlen is unavailable.
bool coversEnoughRelativeTo(const BlastRes& hit, const BlastRes& anchor, double fraction) {
	if (hit.queryCoverageKnown && anchor.queryCoverageKnown) {
		if (hit.queryLength == anchor.queryLength) {
			return hit.queryCovered >= static_cast<int>(std::ceil(static_cast<double>(anchor.queryCovered) * fraction));
		}
		return hit.Qcoverage >= anchor.Qcoverage * fraction;
	}
	return hit.alLen >= static_cast<int>(std::ceil(static_cast<double>(anchor.alLen) * fraction));
}

// True when the subject's lineage contains every rank the consensus reports.
bool supportsLineage(const HitLineage& subject, const TaxObj& consensus) {
	for (int rank = 0; rank < consensus.depth; ++rank) {
		const string& value = consensus.get(rank);
		if (value != __unkwnTax && subject.get(rank) != value) { return false; }
	}
	return true;
}

const string& rankName(const TaxObj* lineage, int rank) { return lineage->get(rank); }
const string& rankName(const HitLineage& lineage, int rank) { return lineage.get(rank); }
int lineageDepth(const TaxObj* lineage) { return lineage->depth; }
int lineageDepth(const HitLineage& lineage) { return lineage.depth; }
bool lineageSpeciesUncertain(const TaxObj* lineage) { return lineage->speciesUncertain; }
bool lineageSpeciesUncertain(const HitLineage& lineage) { return lineage.speciesUncertain(); }
// Vote identity: interned name IDs for reference lineages, the name itself for TaxObj inputs.
uint32_t voteKey(const HitLineage& lineage, int rank) { return lineage.nameId(rank); }
std::string_view voteKey(const TaxObj* lineage, int rank) { return lineage->get(rank); }
bool unknownKey(uint32_t key) { return key == 0; }
bool unknownKey(std::string_view key) { return key.empty(); }

template <class Lineage>
TaxObj* consensusLineage(const vector<Lineage>& TO, bool& hitRd, double LCAfrac, int tdepth) {
	TaxObj* ret = new TaxObj(tdepth);
	hitRd = false;
	const size_t count = TO.size();
	if (count == 0) { return ret; }
	if (count == 1) {
		const int depth = lineageDepth(TO.front());
		ret->SavedTaxs.assign(depth, __unkwnTax);
		for (int rank = 0; rank < depth; ++rank) { ret->SavedTaxs[rank] = rankName(TO.front(), rank); }
		ret->depth = depth;
		ret->speciesUncertain = lineageSpeciesUncertain(TO.front());
		hitRd = true;
		return ret;
	}

	using Key = decltype(voteKey(TO.front(), 0));
	struct Vote {
		Key key;
		size_t example; // a hit naming this taxon
		int count;
	};
	vector<char> active(count, 1);
	// Distinct names at a rank are few, so a linear scan beats hashing.
	vector<Vote> votes;
	votes.reserve(8);
	for (int DL = 0; DL < tdepth; DL++) {
		votes.clear();
		int knownCount = 0;
		for (size_t pos = 0; pos < count; ++pos) {
			if (!active[pos]) { continue; }
			const Key key = voteKey(TO[pos], DL);
			// Hits without a name at this rank abstain but keep voting at deeper ranks.
			if (unknownKey(key)) { continue; }
			++knownCount;
			auto vote = std::find_if(votes.begin(), votes.end(), [&key](const Vote& v) { return v.key == key; });
			if (vote == votes.end()) { votes.push_back(Vote{key, pos, 1}); }
			else { ++vote->count; }
		}
		// A rank missing from every remaining lineage stays unknown; matching deeper names still count.
		if (knownCount == 0) { continue; }

		const Vote* best = &votes.front();
		for (const Vote& vote : votes) {
			if (vote.count > best->count || (vote.count == best->count &&
				rankName(TO[vote.example], DL) < rankName(TO[best->example], DL))) {
				best = &vote;
			}
		}
		if (best->count < static_cast<double>(knownCount) * LCAfrac) { break; }
		const Key consensus = best->key;
		ret->set(DL, rankName(TO[best->example], DL));

		// Hits naming a different taxon at an accepted rank cannot vote at deeper ranks.
		for (size_t pos = 0; pos < count; ++pos) {
			if (!active[pos]) { continue; }
			const Key key = voteKey(TO[pos], DL);
			if (!unknownKey(key) && key != consensus) { active[pos] = 0; }
		}
	}
	return ret;
}
}

TaxObj* LCA(vector<BlastRes>& BR, RefTax* RT, options* opt) {
	if (BR.size() == 0) { return NULL; }

	double bestPerID(0.0);
	filterBlastPrimary(BR, opt, bestPerID);
	if (BR.size() == 0) { return NULL; }

	double avgPerID(0.0);
	const vector<HitLineage> allTax = BlastToTax(BR, RT, opt, avgPerID);

	bool singularHit(false);
	TaxObj* ret = LCAcore(allTax, singularHit, opt->LCAfract, opt->taxDepth);
	ret->Subj = BR.front().Query;
	ret->perID = avgPerID;
	ret->limitByIdentity(avgPerID, opt->idThr);

	if (opt->hitRD) {
		string subject = __unkwnTax;
		if (singularHit) {
			subject = BR.front().Sbj;
		} else if (opt->reportBestHit) {
			// Equally ranked best hits: name the top-ranked subject whose lineage supports the consensus.
			vector<size_t> order(BR.size());
			std::iota(order.begin(), order.end(), 0);
			std::sort(order.begin(), order.end(), [&BR](size_t a, size_t b) { return BR[a].betterThan(BR[b]); });
			subject = BR[order.front()].Sbj;
			for (size_t idx : order) {
				if (supportsLineage(allTax[idx], *ret)) {
					subject = BR[idx].Sbj;
					break;
				}
			}
		}
		ret->setHitDB(subject);
	}
	ret->setRepID(opt->reportID);
	return ret;
}

TaxObj* LCAcore(const vector<HitLineage>& TO, bool& hitRd, double LCAfrac, int tdepth) {
	return consensusLineage(TO, hitRd, LCAfrac, tdepth);
}

TaxObj* LCAcore(const vector<TaxObj*>& TO, bool& hitRd, double LCAfrac, int tdepth) {
	return consensusLineage(TO, hitRd, LCAfrac, tdepth);
}

double filterBlastPrimary(vector<BlastRes>& BR, options* opt, double& bestID) {
	bestID = 0.0;
	if (BR.empty()) { return 0.0; }
	const double minCov = opt->minCover;
	auto passesCoverage = [minCov](const BlastRes& hit) {
		return !hit.queryCoverageKnown || hit.Qcoverage >= minCov;
	};
	const int minAliLen = static_cast<int>(std::ceil(opt->minAliLen));

	if (opt->BLfilter) {
		BR.erase(std::remove_if(BR.begin(), BR.end(), [&](const BlastRes& hit) {
			return hit.alLen < minAliLen || !passesCoverage(hit);
		}), BR.end());
	}
	if (BR.empty()) { return 0.0; }

	// Identity decides the anchor; query coverage, then alignment length, only break identity ties.
	const BlastRes anchor = rankingFields(*std::min_element(BR.begin(), BR.end(),
		[](const BlastRes& a, const BlastRes& b) { return a.betterThan(b); }));
	bestID = anchor.perID;
	const double score = bestID * static_cast<double>(anchor.alLen);

	if (opt->reportBestHit) {
		// All equally ranked best hits are kept, so ties resolve by consensus rather than subject name.
		BR.erase(std::remove_if(BR.begin(), BR.end(), [&anchor](const BlastRes& hit) {
			return hit.compareRank(anchor) != 0;
		}), BR.end());
		deduplicateBlastHits(BR);
		return score;
	}
	if (!opt->BLfilter) {
		deduplicateBlastHits(BR);
		return score;
	}

	//filter parameters
	const double relativeCoverage(0.85);
	double tolerance(1.5);
	if (bestID >= 100) { tolerance = 0.05f; }
	else if (bestID >= 99.5) { tolerance = 0.15f; }
	else if (bestID >= 99) { tolerance = 0.25f; }
	else if (bestID >= 98) { tolerance = 0.75f; }
	else if (bestID >= 97) { tolerance = 1.0f; }

	BR.erase(std::remove_if(BR.begin(), BR.end(), [&](const BlastRes& cur) {
		return (cur.perID + tolerance) < bestID || !coversEnoughRelativeTo(cur, anchor, relativeCoverage);
	}), BR.end());
	// Only surviving alignments may represent a subject in the LCA vote.
	deduplicateBlastHits(BR);
	return score;
}


vector<HitLineage> BlastToTax(const vector<BlastRes>& BR, RefTax* RT, options* opt, double& consPerID) {
	vector<HitLineage> ret;
	ret.reserve(BR.size());
	const vector<double>& thr = opt->idThr;

	consPerID = 0.0;
	if (BR.empty()) { return ret; }
	long double identitySum = 0.0L;
	bool anySpeciesCertain(false);
	for (const BlastRes& hit : BR) {
		identitySum += hit.perID;
		const RefTax::Lineage* reference = RT->find(hit.Sbj);
		if (reference == nullptr) {
			throw FatalError{74, "Could not find tax for Subject " + hit.Sbj + "\n"};
		}
		HitLineage lineage(*RT, *reference);
		lineage.limitByIdentity(hit.perID, thr);
		if (lineage.get(6) != __unkwnTax && !lineage.speciesUncertain()) { anySpeciesCertain = true; }
		ret.push_back(lineage);
	}
	consPerID = static_cast<double>(identitySum / BR.size());

	//remove uncertain species, in case of enough good hits
	if (anySpeciesCertain && ret.size() > 1) {
		for (HitLineage& lineage : ret) {
			if (lineage.depth > 6 && lineage.speciesUncertain()) { lineage.speciesMasked = true; }
		}
	}
	return ret;
}
