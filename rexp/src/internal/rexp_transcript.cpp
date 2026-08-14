#include "rexp_transcript.hpp"

#include "crypto.hpp"

namespace rexp::internal {

Digest rexp_crs_digest(const RawRexpCRS& crs) {
    Bytes bytes;
    frame(bytes, "REXP-G1-CRS-BN254-V1");
    frame(bytes, "BN254");
    frame(bytes, "G1-G2-GT-PAIRING");
    append_u64(bytes, crs.d);
    append_u64(bytes, crs.n);
    frame(bytes, "LEFT-PREFIX-CHAIN-V1");
    for (const auto& point : crs.Gamma) frame_element(bytes, point);
    for (const auto& point : crs.Lambda) frame_element(bytes, point);
    return sha256(bytes);
}

Digest rexp_statement_digest(
    const PreparedPublicParameters& params,
    const std::vector<G1>& points,
    const GT& d1,
    const GT& e,
    const GT& f,
    const GT& tl,
    const GT& tr) {
    Bytes bytes;
    frame(bytes, "REXP-G1-STATEMENT-BN254-V1");
    frame_digest(bytes, params.digest());
    append_u64(bytes, params.d());
    append_u64(bytes, params.n());
    for (const auto& point : points) frame_element(bytes, point);
    frame_element(bytes, d1);
    frame_element(bytes, e);
    frame_element(bytes, f);
    frame_element(bytes, tl);
    frame_element(bytes, tr);
    return sha256(bytes);
}

Digest rexp_initial_transcript(
    const PreparedPublicParameters& params,
    const PreparedStatement& statement) {
    Bytes bytes;
    frame(bytes, "REXP-G1-FS-v1");
    frame(bytes, "BN254");
    frame(bytes, "G1");
    frame_digest(bytes, params.digest());
    frame_digest(bytes, statement.digest());
    append_u64(bytes, params.d());
    append_u64(bytes, params.n());
    return sha256(bytes);
}

Digest rexp_absorb_round(
    const Digest& transcript,
    std::size_t round,
    std::size_t dimension,
    const RexpRoundMessage& message) {
    Bytes bytes;
    frame(bytes, "REXP-G1-ROUND-MESSAGE-V1");
    frame_digest(bytes, transcript);
    append_u64(bytes, round);
    append_u64(bytes, dimension);
    frame_element(bytes, message.E);
    frame_element(bytes, message.F);
    frame_element(bytes, message.TL);
    frame_element(bytes, message.TR);
    return sha256(bytes);
}

Digest rexp_enter_dory(
    const Digest& transcript,
    std::size_t round,
    std::size_t half_dimension,
    const DoryStatement& statement) {
    Bytes bytes;
    frame(bytes, "REXP-G1-EMBEDDED-DORY-ENTER-V1");
    frame_digest(bytes, transcript);
    append_u64(bytes, round);
    append_u64(bytes, round + 1);
    append_u64(bytes, half_dimension);
    frame_element(bytes, statement.D0);
    frame_element(bytes, statement.D1);
    frame_element(bytes, statement.D2);
    return sha256(bytes);
}

Digest rexp_leave_dory(const Digest& transcript, std::size_t round) {
    Bytes bytes;
    frame(bytes, "REXP-G1-EMBEDDED-DORY-EXIT-V1");
    frame_digest(bytes, transcript);
    append_u64(bytes, round);
    return sha256(bytes);
}

Digest rexp_absorb_final(const Digest& transcript, const G1& result) {
    Bytes bytes;
    frame(bytes, "REXP-G1-FINAL-OUTPUT-V1");
    frame_digest(bytes, transcript);
    frame_element(bytes, result);
    return sha256(bytes);
}

RexpRoundMessage rexp_initial_round_message(
    const PreparedStatement& statement) {
    return {statement.E0(), statement.F0(), statement.TL0(), statement.TR0()};
}

} 
