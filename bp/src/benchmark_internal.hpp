#pragma once

#include "bp/bp.hpp"

namespace bp::benchmark_internal {

struct Profile {
  double total_ms{};
  double public_parameter_validation_ms{};
  double proof_parsing_validation_ms{};
  double initial_allocation_copy_ms{};
  double msm_input_copy_ms{};
  double proof_message_msm_ms{};
  double generator_folding_ms{};
  double witness_folding_ms{};
  double transcript_challenge_ms{};
  double verifier_z_update_ms{};
  double terminal_verification_ms{};

  Profile& operator+=(const Profile& other);
  Profile& operator/=(double divisor);
  double accounted_ms() const;
  double residual_overhead_ms() const;
};

struct PrevalidatedVerification {
  const PublicParams* pp{};
  const Group* Z{};
  const Proof* proof{};
};

Proof ProveProfiled(const PublicParams& pp, const Group& Z,
                    std::span<const Scalar> x, std::span<const Scalar> y,
                    Profile& profile);
bool VerifySerializedProfiled(const PublicParams& pp, const Group& Z,
                              std::span<const std::uint8_t> proof_bytes,
                              Profile& profile);
bool PrevalidateVerification(const PublicParams& pp, const Group& Z,
                             const Proof& proof,
                             PrevalidatedVerification& output) noexcept;
bool VerifyPrevalidatedProfiled(const PrevalidatedVerification& inputs,
                                Profile& profile) noexcept;

}
