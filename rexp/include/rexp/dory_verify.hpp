#pragma once

#include "rexp/dory.hpp"

#include <cstddef>
#include <optional>

namespace rexp {

struct VerifyMetrics {
    std::size_t terminal_pairings = 0;
    std::size_t symbolic_atom_insertions = 0;
    std::size_t coalesced_duplicate_bases = 0;
    std::size_t zero_coefficients_removed = 0;
    std::size_t identity_bases_removed = 0;
    std::size_t actual_gt_bases = 0;
    std::size_t nonzero_gt_coefficients = 0;
    std::size_t coefficient_one_bases = 0;
    std::size_t gt_multiexp_calls = 0;
    const char* gt_multiexp_backend = "mcl::bn::GT::powVec";
    double transcript_ms = 0;
    double batch_inversion_ms = 0;
    double symbolic_replay_ms = 0;
    double gt_multiexp_ms = 0;
    double terminal_pairing_ms = 0;
    double total_ms = 0;
};

struct VerifyAuditTrace {
    std::vector<GT> gt_bases;
    std::vector<Fr> gt_exponents;
    G1 terminal_left_g1;
    G2 terminal_right_g2;
    GT evaluated_lhs;
    GT terminal_pairing_rhs;
};



class PreparedVerifier {
public:
    const DoryCRS& crs() const { return *crs_; }
    const DoryPrecomputation& precomp() const { return *precomp_; }

private:
    PreparedVerifier(
        const DoryCRS& crs,
        const DoryPrecomputation& precomp)
        : crs_(&crs), precomp_(&precomp) {}
    const DoryCRS* crs_;
    const DoryPrecomputation* precomp_;
    friend std::optional<PreparedVerifier> PrepareVerifier(
        const DoryCRS&, const DoryPrecomputation&);
};

std::optional<PreparedVerifier> PrepareVerifier(
    const DoryCRS& crs,
    const DoryPrecomputation& precomp);

bool Verify(
    const DoryCRS& crs,
    const DoryPrecomputation& precomp,
    const DoryStatement& statement,
    const DoryProof& proof,
    VerifyMetrics* metrics = nullptr);

bool VerifyBatch(
    const DoryCRS& crs,
    const DoryPrecomputation& precomp,
    const std::vector<DoryStatement>& statements,
    const DoryBatchProof& proof,
    VerifyMetrics* metrics = nullptr);

bool VerifyPrepared(
    const PreparedVerifier& prepared,
    const DoryStatement& statement,
    const DoryProof& proof,
    VerifyMetrics* metrics = nullptr);

bool VerifyEmbeddedDeferred(
    const DoryCRS& crs,
    const DoryPrecomputation& precomp,
    const DoryStatement& statement,
    const DoryProof& proof,
    const Digest& transcript_in,
    Digest* transcript_end = nullptr,
    VerifyMetrics* metrics = nullptr);



bool VerifyEmbeddedReference(
    const DoryCRS& crs,
    const DoryPrecomputation& precomp,
    const DoryStatement& statement,
    const DoryProof& proof,
    const Digest& transcript_in,
    Digest* transcript_end = nullptr);

bool VerifyBatchPrepared(
    const PreparedVerifier& prepared,
    const std::vector<DoryStatement>& statements,
    const DoryBatchProof& proof,
    VerifyMetrics* metrics = nullptr);



bool VerifyPreparedAudit(
    const PreparedVerifier& prepared,
    const DoryStatement& statement,
    const DoryProof& proof,
    VerifyMetrics& metrics,
    VerifyAuditTrace& trace);

bool VerifyBatchPreparedAudit(
    const PreparedVerifier& prepared,
    const std::vector<DoryStatement>& statements,
    const DoryBatchProof& proof,
    VerifyMetrics& metrics,
    VerifyAuditTrace& trace);

}
