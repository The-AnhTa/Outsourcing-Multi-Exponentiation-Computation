#include "hp_vme_internal.hpp"

#include <mcl/fp.hpp>
#include <mcl/gmp_util.hpp>
#include <algorithm>
#include <array>
#include <chrono>
#include <compare>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <string_view>

namespace bp::hp_internal {
namespace {

using Clock = std::chrono::steady_clock;

constexpr Digest kRejectionLimit = {
    0xde,0xd4,0x5b,0x0d,0x80,0x00,0x00,0x0a,
    0x5d,0x39,0xd1,0x00,0x00,0x00,0x00,0x2f,
    0xfd,0xbd,0x00,0x00,0x00,0x00,0x00,0x63,
    0xc6,0x00,0x00,0x00,0x00,0x00,0x00,0x4e};

double ms(Clock::time_point a, Clock::time_point b = Clock::now()) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}

Bytes u64be(std::uint64_t value) {
  Bytes out;
  for (int shift = 56; shift >= 0; shift -= 8)
    out.push_back(static_cast<std::uint8_t>(value >> shift));
  return out;
}

void raw(Bytes& out, std::span<const std::uint8_t> field) {
  out.insert(out.end(), field.begin(), field.end());
}

void frame(Bytes& out, std::span<const std::uint8_t> field) {
  raw(out, u64be(field.size()));
  raw(out, field);
}

void frame(Bytes& out, std::string_view field) {
  frame(out, {reinterpret_cast<const std::uint8_t*>(field.data()), field.size()});
}

template<class T>
Bytes encode(const T& value) {
  Bytes out(1024);
  const auto written = value.serialize(out.data(), out.size());
  if (!written) throw std::runtime_error("MCL serialization failed");
  out.resize(written);
  return out;
}

bool power_of_two(std::size_t n) { return n && !(n & (n - 1)); }

std::size_t exact_log2(std::size_t n) {
  if (!power_of_two(n)) throw std::invalid_argument("dimension is not a power of two");
  std::size_t d = 0;
  while (n > 1) { n >>= 1; ++d; }
  return d;
}

Scalar fadd(const Scalar& a, const Scalar& b) {
  Scalar out; Scalar::add(out, a, b); return out;
}
Scalar fmul(const Scalar& a, const Scalar& b) {
  Scalar out; Scalar::mul(out, a, b); return out;
}
Scalar finv(const Scalar& a) {
  if (a.isZero()) throw std::invalid_argument("zero scalar inversion");
  Scalar out; Scalar::inv(out, a); return out;
}
Scalar fneg(const Scalar& a) {
  Scalar out; Scalar::neg(out, a); return out;
}
G1 g1add(const G1& a, const G1& b) {
  G1 out; G1::add(out, a, b); return out;
}
G1 g1mul(const G1& p, const Scalar& x) {
  G1 out; G1::mul(out, p, x); return out;
}
Group g2add(const Group& a, const Group& b) {
  Group out; Group::add(out, a, b); return out;
}
Group g2mul(const Group& p, const Scalar& x) {
  Group out; Group::mul(out, p, x); return out;
}
GT gtmul(const GT& a, const GT& b) {
  GT out; GT::mul(out, a, b); return out;
}
GT gtpow(const GT& a, const Scalar& x) {
  GT out; GT::pow(out, a, x); return out;
}
GT gtinv(const GT& a) {
  GT out; GT::inv(out, a); return out;
}
GT product(std::initializer_list<GT> values) {
  GT out; out.setOne();
  for (const auto& value : values) out = gtmul(out, value);
  return out;
}

bool valid_g1(const G1& p, bool nonidentity = false) {
  return p.isValid() && p.isValidOrder() && (!nonidentity || !p.isZero());
}
bool valid_g2(const Group& p, bool nonidentity = false) {
  return p.isValid() && p.isValidOrder() && (!nonidentity || !p.isZero());
}

bool valid_gt(const GT& value) {
  try {
    const Bytes bytes = serialize_gt(value);
    GT decoded;
    if (decoded.deserialize(bytes.data(), bytes.size()) != bytes.size() ||
        decoded != value)
      return false;
    mpz_class order;
    mcl::gmp::setStr(order, Scalar::getModulo(), 10);
    GT powered;
    GT::pow(powered, value, order);
    GT one; one.setOne();
    return powered == one;
  } catch (...) { return false; }
}

GT pairing(const G1& a, const Group& b) {
  GT out; mcl::bn::pairing(out, a, b); return out;
}

GT pairing_product(std::span<const G1> a, std::span<const Group> b) {
  if (a.size() != b.size()) throw std::invalid_argument("pairing length mismatch");
  GT out;
  if (a.empty()) { out.setOne(); return out; }
  GT miller;
  mcl::bn::millerLoopVec(miller, a.data(), b.data(), a.size(), true);
  mcl::bn::finalExp(out, miller);
  return out;
}

Scalar inner(std::span<const Scalar> a, std::span<const Scalar> b) {
  if (a.size() != b.size()) throw std::invalid_argument("inner-product length mismatch");
  Scalar out; out.clear();
  for (std::size_t i = 0; i < a.size(); ++i)
    out = fadd(out, fmul(a[i], b[i]));
  return out;
}

std::vector<Scalar> tensor_vector(std::span<const Scalar> challenges) {
  std::vector<Scalar> out(1);
  out[0] = 1;
  for (std::size_t j = challenges.size(); j > 0; --j) {
    std::vector<Scalar> next;
    next.reserve(out.size() * 2);
    for (const auto& value : out) {
      next.push_back(value);
      next.push_back(fmul(value, challenges[j - 1]));
    }
    out = std::move(next);
  }
  return out;
}

Digest digest_scalars(std::span<const Scalar> values) {
  Bytes input;
  frame(input, "BPVME/HP/VME/Z-DIGEST/v1");
  for (const auto& value : values) frame(input, serialize(value));
  return sha256(input);
}

Digest digest_points(std::span<const Group> points) {
  Bytes input;
  frame(input, "BPVME/HP/VME/P-DIGEST/v1");
  for (const auto& point : points) frame(input, serialize(point));
  return sha256(input);
}

Digest statement_digest(const VmePublicParams& pp, const Digest& context,
                        const Group& X, std::span<const Scalar> z) {
  Bytes input;
  frame(input, "BPVME/HP/VME/STATEMENT/v1");
  frame(input, pp.transcript_domain);
  frame(input, pp.digest);
  frame(input, context);
  frame(input, u64be(pp.dimension));
  frame(input, digest_points(pp.fixed_P));
  frame(input, serialize(X));
  frame(input, digest_scalars(z));
  return sha256(input);
}

class VmeTranscript {
 public:
  explicit VmeTranscript(const Digest& statement) {
    Bytes input;
    frame(input, "BPVME/HP/VME-FS/v1");
    frame(input, statement);
    state_ = sha256(input);
  }

  void absorb(std::string_view label, std::span<const Bytes> fields) {
    Bytes input;
    frame(input, "BPVME/HP/VME-ABSORB/v1");
    frame(input, state_);
    frame(input, label);
    for (const auto& field : fields) frame(input, field);
    state_ = sha256(input);
  }

  void absorb(std::string_view label, const Bytes& field) {
    absorb(label, std::span<const Bytes>(&field, 1));
  }

  Scalar challenge(std::string_view label, std::uint64_t index) {
    for (std::uint64_t counter = 0;; ++counter) {
      Bytes input;
      frame(input, "BPVME/HP/VME-CHALLENGE/v1");
      frame(input, state_);
      frame(input, label);
      frame(input, u64be(index));
      frame(input, u64be(counter));
      const Digest hash = sha256(input);
      if (!std::lexicographical_compare(
              hash.begin(), hash.end(), kRejectionLimit.begin(), kRejectionLimit.end()))
        continue;
      Scalar out;
      out.setBigEndianMod(hash.data(), hash.size());
      if (out.isZero()) continue;
      absorb("challenge-value", serialize(out));
      return out;
    }
  }

 private:
  Digest state_{};
};

G1 derive_g1(std::string_view domain, std::span<const std::uint8_t> seed,
             std::size_t index) {
  Bytes input;
  frame(input, domain);
  frame(input, seed);
  frame(input, u64be(index));
  G1 out;
  mcl::bn::hashAndMapToG1(out, input.data(), input.size());
  if (!valid_g1(out, true)) throw std::runtime_error("invalid derived G1");
  return out;
}

Group derive_g2(std::string_view domain, std::span<const std::uint8_t> seed) {
  Bytes input;
  frame(input, domain);
  frame(input, seed);
  Group out;
  mcl::bn::hashAndMapToG2(out, input.data(), input.size());
  if (!valid_g2(out, true)) throw std::runtime_error("invalid derived G2");
  return out;
}

Digest crs_digest(const VmePublicParams& pp) {
  Bytes input;
  frame(input, "BPVME/HP/VME/CRS/v1");
  frame(input, kGroupIdentifier);
  frame(input, kScalarModulus);
  frame(input, u64be(pp.dimension));
  frame(input, u64be(pp.log_dimension));
  frame(input, pp.transcript_domain);
  for (const auto& p : pp.auxiliary_G) frame(input, serialize_g1(p));
  for (const auto& p : pp.fixed_P) frame(input, serialize(p));
  frame(input, serialize_g1(pp.L));
  frame(input, serialize(pp.Lprime));
  return sha256(input);
}

Digest precomp_digest(const VmePublicParams& pp, const VmePrecomputation& pre) {
  Bytes input;
  frame(input, "BPVME/HP/VME/PRECOMPUTATION/v1");
  frame(input, pp.digest);
  for (const auto& x : pre.pairing_x) frame(input, serialize_gt(x));
  for (const auto& x : pre.delta1R) frame(input, serialize_gt(x));
  for (const auto& x : pre.delta2R) frame(input, serialize_gt(x));
  frame(input, serialize_gt(pre.pairing_LLprime));
  return sha256(input);
}

struct FreshInstance {
  GT D0, D1, D2;
  std::vector<G1> phi;
  std::vector<Group> theta;
};

struct Target { GT D0, D1, D2; };

void absorb_rexp(VmeTranscript& transcript, std::size_t j, std::size_t m,
                 const VmeRexpClaims& claims) {
  std::array<Bytes, 6> fields{u64be(j), u64be(m), serialize_gt(claims.E),
      serialize_gt(claims.F), serialize_gt(claims.TL), serialize_gt(claims.TR)};
  transcript.absorb("BPVME/HP/VME/REXP-CLAIMS/v1", fields);
}

VmeRexpClaims fixed_rexp_claim(const VmePrecomputation& pre) {
  return {pre.delta1R[0], pre.delta2R[0], pre.pairing_x[1], pre.delta1R[0]};
}

struct Phase1 {
  VmeStatement statement;
  std::vector<Scalar> rho;
  std::vector<Scalar> r;
  std::vector<VmeRexpClaims> claims;
  std::vector<FreshInstance> fresh;
  G1 R;
  Digest transcript_state{};
};

Phase1 prove_phase1(const VmePublicParams& pp, const VmePrecomputation& pre,
                    const Digest& context, const Group& X,
                    std::span<const Scalar> z) {
  Phase1 out;
  out.statement.application_context = context;
  out.statement.X = X;
  out.statement.z.assign(z.begin(), z.end());
  out.statement.digest = statement_digest(pp, context, X, z);
  VmeTranscript transcript(out.statement.digest);
  out.rho.resize(pp.log_dimension);
  out.fresh.resize(pp.log_dimension + 1);
  std::vector<G1> current = pp.auxiliary_G;
  GT outer = pre.pairing_x[0];
  for (std::size_t j = 0; j < pp.log_dimension; ++j) {
    const std::size_t m = pp.dimension >> j, half = m / 2, t = j + 1;
    VmeRexpClaims claims;
    if (j == 0) {
      claims = fixed_rexp_claim(pre);
    } else {
      claims.E = pairing_product(std::span(current).subspan(half, half),
                                 std::span(pp.fixed_P).first(half));
      claims.F = pairing_product(std::span(current).first(half),
                                 std::span(pp.fixed_P).subspan(half, half));
      claims.TL = pairing_product(std::span(current).first(half),
                                  std::span(pp.fixed_P).first(half));
      claims.TR = pairing_product(std::span(current).subspan(half, half),
                                  std::span(pp.fixed_P).first(half));
      out.claims.push_back(claims);
    }
    absorb_rexp(transcript, j, m, claims);
    const Scalar rho = transcript.challenge("BPVME/HP/VME/REXP-RHO/v1", j);
    const Scalar rho_inv = finv(rho);
    out.rho[j] = rho;
    std::vector<G1> next;
    std::vector<Group> theta;
    next.reserve(half); theta.reserve(half);
    for (std::size_t i = 0; i < half; ++i) {
      next.push_back(g1add(current[i], g1mul(current[half + i], rho)));
      theta.push_back(g2add(pp.fixed_P[i], g2mul(pp.fixed_P[half + i], rho_inv)));
    }
    auto& fresh = out.fresh[t];
    fresh.D0 = product({outer, gtpow(claims.E, rho), gtpow(claims.F, rho_inv)});
    fresh.D1 = gtmul(claims.TL, gtpow(claims.TR, rho));
    fresh.D2 = gtmul(pre.pairing_x[t], gtpow(pre.delta2R[j], rho_inv));
    fresh.phi = next;
    fresh.theta = std::move(theta);
    current = std::move(next);
    outer = fresh.D1;
  }
  out.R = current[0];
  transcript.absorb("BPVME/HP/VME/REXP-R/v1", serialize_g1(out.R));
  out.r.resize(pp.log_dimension);
  for (std::size_t j = 0; j < pp.log_dimension; ++j)
    out.r[pp.log_dimension - j - 1] = out.rho[j];
  return out;
}

VmeProof prove_phase2(const VmePublicParams& pp, const VmePrecomputation& pre,
                      const Phase1& phase1) {
  VmeProof proof;
  proof.rexp_claims = phase1.claims;
  proof.R = phase1.R;
  proof.dory_folds.reserve(pp.log_dimension);
  proof.batch_U.reserve(pp.log_dimension);
  VmeTranscript transcript(phase1.statement.digest);
  for (std::size_t j = 0; j < pp.log_dimension; ++j) {
    const VmeRexpClaims claims =
        j == 0 ? fixed_rexp_claim(pre) : phase1.claims[j - 1];
    absorb_rexp(transcript, j, pp.dimension >> j, claims);
    (void)transcript.challenge("BPVME/HP/VME/REXP-RHO/v1", j);
  }
  transcript.absorb("BPVME/HP/VME/REXP-R/v1", serialize_g1(phase1.R));
  const auto tensor = tensor_vector(phase1.r);
  const Scalar q = inner(tensor, phase1.statement.z);
  std::array<Bytes, 3> initial{u64be(pp.log_dimension), u64be(pp.dimension),
                               serialize(q)};
  transcript.absorb("BPVME/HP/VME/DORY-INITIAL/v1", initial);

  std::vector<G1> phi;
  std::vector<Group> theta;
  phi.reserve(pp.dimension); theta.reserve(pp.dimension);
  for (std::size_t i = 0; i < pp.dimension; ++i) {
    phi.push_back(g1mul(pp.L, phase1.statement.z[i]));
    theta.push_back(g2mul(pp.Lprime, tensor[i]));
  }
  Target aggregate{gtpow(pre.pairing_LLprime, q),
                   pairing(pp.L, phase1.statement.X),
                   pairing(phase1.R, pp.Lprime)};

  for (std::size_t t = 1; t <= pp.log_dimension; ++t) {
    const std::size_t k = t - 1, m = pp.dimension >> k, half = m / 2;
    VmeDoryFold fold;
    fold.D1L = pairing_product(std::span(phi).first(half),
                               std::span(pp.fixed_P).first(half));
    fold.D1R = pairing_product(std::span(phi).subspan(half, half),
                               std::span(pp.fixed_P).first(half));
    fold.D2L = pairing_product(std::span(pp.auxiliary_G).first(half),
                               std::span(theta).first(half));
    fold.D2R = pairing_product(std::span(pp.auxiliary_G).first(half),
                               std::span(theta).subspan(half, half));
    std::array<Bytes, 6> beta_fields{u64be(k), u64be(m), serialize_gt(fold.D1L),
        serialize_gt(fold.D1R), serialize_gt(fold.D2L), serialize_gt(fold.D2R)};
    transcript.absorb("BPVME/HP/VME/DORY-BETA-MESSAGE/v1", beta_fields);
    const Scalar beta = transcript.challenge("BPVME/HP/VME/DORY-BETA/v1", k);
    const Scalar beta_inv = finv(beta);

    std::vector<G1> phi_cross;
    std::vector<Group> theta_cross;
    phi_cross.reserve(m); theta_cross.reserve(m);
    for (std::size_t i = 0; i < m; ++i) {
      phi_cross.push_back(g1add(phi[i], g1mul(pp.auxiliary_G[i], beta)));
      theta_cross.push_back(g2add(theta[i], g2mul(pp.fixed_P[i], beta_inv)));
    }
    fold.W1 = pairing_product(std::span(phi_cross).first(half),
                              std::span(theta_cross).subspan(half, half));
    fold.W2 = pairing_product(std::span(phi_cross).subspan(half, half),
                              std::span(theta_cross).first(half));
    std::array<Bytes, 4> alpha_fields{u64be(k), u64be(m),
                                      serialize_gt(fold.W1), serialize_gt(fold.W2)};
    transcript.absorb("BPVME/HP/VME/DORY-ALPHA-MESSAGE/v1", alpha_fields);
    const Scalar alpha = transcript.challenge("BPVME/HP/VME/DORY-ALPHA/v1", k);
    const Scalar alpha_inv = finv(alpha);

    std::vector<G1> folded_phi;
    std::vector<Group> folded_theta;
    folded_phi.reserve(half); folded_theta.reserve(half);
    for (std::size_t i = 0; i < half; ++i) {
      folded_phi.push_back(g1add(g1mul(phi_cross[i], alpha), phi_cross[half + i]));
      folded_theta.push_back(
          g2add(g2mul(theta_cross[i], alpha_inv), theta_cross[half + i]));
    }
    Target bar;
    bar.D0 = product({aggregate.D0, pre.pairing_x[k],
        gtpow(aggregate.D1, beta_inv), gtpow(aggregate.D2, beta),
        gtpow(fold.W1, alpha), gtpow(fold.W2, alpha_inv)});
    bar.D1 = product({gtpow(fold.D1L, alpha), fold.D1R,
        gtpow(pre.pairing_x[k + 1], fmul(alpha, beta)),
        gtpow(pre.delta1R[k], beta)});
    bar.D2 = product({gtpow(fold.D2L, alpha_inv), fold.D2R,
        gtpow(pre.pairing_x[k + 1], fmul(alpha_inv, beta_inv)),
        gtpow(pre.delta2R[k], beta_inv)});

    const auto& fresh = phase1.fresh[t];
    const GT U = gtmul(pairing_product(folded_phi, fresh.theta),
                       pairing_product(fresh.phi, folded_theta));
    proof.batch_U.push_back(U);
    std::array<Bytes, 3> u_fields{u64be(t), u64be(half), serialize_gt(U)};
    transcript.absorb("BPVME/HP/VME/BATCH-U/v1", u_fields);
    const Scalar gamma = transcript.challenge("BPVME/HP/VME/DORY-GAMMA/v1", t);
    aggregate.D0 = product({gtpow(bar.D0, fmul(gamma, gamma)),
                            gtpow(U, gamma), fresh.D0});
    aggregate.D1 = gtmul(gtpow(bar.D1, gamma), fresh.D1);
    aggregate.D2 = gtmul(gtpow(bar.D2, gamma), fresh.D2);
    for (std::size_t i = 0; i < half; ++i) {
      folded_phi[i] = g1add(g1mul(folded_phi[i], gamma), fresh.phi[i]);
      folded_theta[i] = g2add(g2mul(folded_theta[i], gamma), fresh.theta[i]);
    }
    phi = std::move(folded_phi);
    theta = std::move(folded_theta);
    proof.dory_folds.push_back(fold);
  }
  proof.phi_final = phi[0];
  proof.theta_final = theta[0];
  return proof;
}

struct Replay {
  VmeTranscript transcript;
  std::vector<Scalar> rho, beta, alpha, gamma;
  Scalar epsilon, q;
  std::vector<Target> fresh;
  Target aggregate;
  GT outer;

  explicit Replay(const Digest& digest) : transcript(digest) {}
};

enum class GtAtomKind {
  PairingX, Delta1R, Delta2R, PairingLLprime,
  RexpE, RexpF, RexpTL, RexpTR,
  DoryD1L, DoryD1R, DoryD2L, DoryD2R, DoryW1, DoryW2, BatchU
};
struct GtAtomId {
  GtAtomKind kind{};
  std::size_t index{};
  auto operator<=>(const GtAtomId&) const = default;
};
enum class PairAtomKind {
  InitialLX, InitialRLprime, TerminalDory, TerminalRexp
};
struct PairAtomId {
  PairAtomKind kind{};
  auto operator<=>(const PairAtomId&) const = default;
};
struct SymbolicExpression {
  std::map<GtAtomId, Scalar> gt;
  std::map<PairAtomId, Scalar> pairs;
};
struct PairInputs { G1 first; Group second; };

Scalar fone() { Scalar out; out = 1; return out; }

template<class K>
void add_coefficient(std::map<K, Scalar>& terms, const K& id,
                     const Scalar& coefficient) {
  auto [it, inserted] = terms.try_emplace(id);
  if (inserted) it->second.clear();
  it->second = fadd(it->second, coefficient);
}

SymbolicExpression gt_atom(GtAtomKind kind, std::size_t index = 0) {
  SymbolicExpression out;
  out.gt[{kind, index}] = fone();
  return out;
}

SymbolicExpression pair_atom(PairAtomKind kind) {
  SymbolicExpression out;
  out.pairs[{kind}] = fone();
  return out;
}

void multiply(SymbolicExpression& target, const SymbolicExpression& source) {
  for (const auto& [id, coefficient] : source.gt)
    add_coefficient(target.gt, id, coefficient);
  for (const auto& [id, coefficient] : source.pairs)
    add_coefficient(target.pairs, id, coefficient);
}

SymbolicExpression powered(SymbolicExpression value, const Scalar& exponent) {
  for (auto& [id, coefficient] : value.gt)
    coefficient = fmul(coefficient, exponent);
  for (auto& [id, coefficient] : value.pairs)
    coefficient = fmul(coefficient, exponent);
  return value;
}

void multiply_powered(SymbolicExpression& target,
                      const SymbolicExpression& source,
                      const Scalar& exponent) {
  multiply(target, powered(source, exponent));
}

void normalize(SymbolicExpression& value) {
  for (auto it = value.gt.begin(); it != value.gt.end();)
    if (it->second.isZero()) it = value.gt.erase(it); else ++it;
  for (auto it = value.pairs.begin(); it != value.pairs.end();)
    if (it->second.isZero()) it = value.pairs.erase(it); else ++it;
}

std::vector<Scalar> batch_invert(std::span<const Scalar> values) {
  std::vector<Scalar> prefix(values.size()), result(values.size());
  Scalar accumulator = fone();
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (values[i].isZero()) throw std::invalid_argument("zero batch inverse");
    prefix[i] = accumulator;
    accumulator = fmul(accumulator, values[i]);
  }
  Scalar inverse = finv(accumulator);
  for (std::size_t i = values.size(); i-- > 0;) {
    result[i] = fmul(inverse, prefix[i]);
    inverse = fmul(inverse, values[i]);
  }
  return result;
}

unsigned hex_nibble(char c) {
  if (c >= '0' && c <= '9') return unsigned(c - '0');
  if (c >= 'a' && c <= 'f') return unsigned(c - 'a' + 10);
  if (c >= 'A' && c <= 'F') return unsigned(c - 'A' + 10);
  throw std::runtime_error("invalid scalar hex");
}

unsigned scalar_window_digit(const std::string& hex, std::size_t offset,
                             std::size_t width) {
  unsigned out = 0;
  for (std::size_t bit = 0; bit < width; ++bit) {
    const std::size_t absolute = offset + bit, nibble = absolute / 4;
    if (nibble < hex.size())
      out |= ((hex_nibble(hex[hex.size() - 1 - nibble]) >>
               (absolute % 4)) & 1U) << bit;
  }
  return out;
}

GT gt_multiexp(std::span<const GT> bases,
               std::span<const Scalar> exponents) {
  if (bases.size() != exponents.size())
    throw std::invalid_argument("GT multiexp length mismatch");
  GT one; one.setOne();
  if (bases.empty()) return one;
  const std::size_t width =
      bases.size() <= 4 ? 2 : bases.size() <= 16 ? 3 :
      bases.size() <= 160 ? 4 : bases.size() <= 256 ? 5 : 6;
  std::vector<std::string> hex;
  std::size_t bits = 0;
  hex.reserve(exponents.size());
  for (const auto& exponent : exponents) {
    hex.push_back(exponent.getStr(16));
    bits = std::max(bits, hex.back().size() * 4);
  }
  const std::size_t windows = (bits + width - 1) / width;
  const std::size_t bucket_count = std::size_t{1} << width;
  std::vector<GT> buckets(bucket_count);
  GT accumulator = one;
  for (std::size_t window = windows; window-- > 0;) {
    for (std::size_t i = 0; i < width; ++i)
      accumulator = gtmul(accumulator, accumulator);
    for (auto& bucket : buckets) bucket = one;
    for (std::size_t i = 0; i < bases.size(); ++i) {
      const unsigned digit =
          scalar_window_digit(hex[i], window * width, width);
      if (digit) buckets[digit] = gtmul(buckets[digit], bases[i]);
    }
    GT running = one;
    for (std::size_t i = bucket_count; i-- > 1;) {
      running = gtmul(running, buckets[i]);
      accumulator = gtmul(accumulator, running);
    }
  }
  return accumulator;
}

bool replay(const VmePublicParams& pp, const VmePrecomputation& pre,
            const VmeStatement& statement, const VmeProof& proof, Replay& out,
            VmeVerificationStats* stats) {
  const auto start = Clock::now();
  out.rho.resize(pp.log_dimension);
  out.fresh.resize(pp.log_dimension + 1);
  out.outer = pre.pairing_x[0];
  for (std::size_t j = 0; j < pp.log_dimension; ++j) {
    const VmeRexpClaims claims =
        j == 0 ? fixed_rexp_claim(pre) : proof.rexp_claims[j - 1];
    absorb_rexp(out.transcript, j, pp.dimension >> j, claims);
    const Scalar rho = out.transcript.challenge("BPVME/HP/VME/REXP-RHO/v1", j);
    const Scalar rho_inv = finv(rho);
    out.rho[j] = rho;
    auto& fresh = out.fresh[j + 1];
    fresh.D0 = product({out.outer, gtpow(claims.E, rho),
                        gtpow(claims.F, rho_inv)});
    fresh.D1 = gtmul(claims.TL, gtpow(claims.TR, rho));
    fresh.D2 = gtmul(pre.pairing_x[j + 1],
                     gtpow(pre.delta2R[j], rho_inv));
    out.outer = fresh.D1;
  }
  out.transcript.absorb("BPVME/HP/VME/REXP-R/v1", serialize_g1(proof.R));
  std::vector<Scalar> r(pp.log_dimension);
  for (std::size_t j = 0; j < pp.log_dimension; ++j)
    r[pp.log_dimension - j - 1] = out.rho[j];
  const auto tensor = tensor_vector(r);
  out.q = inner(tensor, statement.z);
  std::array<Bytes, 3> initial{u64be(pp.log_dimension), u64be(pp.dimension),
                               serialize(out.q)};
  out.transcript.absorb("BPVME/HP/VME/DORY-INITIAL/v1", initial);
  out.aggregate = {gtpow(pre.pairing_LLprime, out.q),
                   pairing(pp.L, statement.X), pairing(proof.R, pp.Lprime)};
  out.beta.resize(pp.log_dimension);
  out.alpha.resize(pp.log_dimension);
  out.gamma.resize(pp.log_dimension + 1);
  if (stats) stats->transcript_ms += ms(start);

  const auto recurrence_start = Clock::now();
  for (std::size_t t = 1; t <= pp.log_dimension; ++t) {
    const std::size_t k = t - 1, m = pp.dimension >> k, half = m / 2;
    const auto& fold = proof.dory_folds[k];
    std::array<Bytes, 6> beta_fields{u64be(k), u64be(m), serialize_gt(fold.D1L),
        serialize_gt(fold.D1R), serialize_gt(fold.D2L), serialize_gt(fold.D2R)};
    out.transcript.absorb("BPVME/HP/VME/DORY-BETA-MESSAGE/v1", beta_fields);
    const Scalar beta =
        out.transcript.challenge("BPVME/HP/VME/DORY-BETA/v1", k);
    const Scalar beta_inv = finv(beta);
    out.beta[k] = beta;
    std::array<Bytes, 4> alpha_fields{u64be(k), u64be(m),
                                      serialize_gt(fold.W1), serialize_gt(fold.W2)};
    out.transcript.absorb("BPVME/HP/VME/DORY-ALPHA-MESSAGE/v1", alpha_fields);
    const Scalar alpha =
        out.transcript.challenge("BPVME/HP/VME/DORY-ALPHA/v1", k);
    const Scalar alpha_inv = finv(alpha);
    out.alpha[k] = alpha;
    Target bar;
    bar.D0 = product({out.aggregate.D0, pre.pairing_x[k],
        gtpow(out.aggregate.D1, beta_inv), gtpow(out.aggregate.D2, beta),
        gtpow(fold.W1, alpha), gtpow(fold.W2, alpha_inv)});
    bar.D1 = product({gtpow(fold.D1L, alpha), fold.D1R,
        gtpow(pre.pairing_x[k + 1], fmul(alpha, beta)),
        gtpow(pre.delta1R[k], beta)});
    bar.D2 = product({gtpow(fold.D2L, alpha_inv), fold.D2R,
        gtpow(pre.pairing_x[k + 1], fmul(alpha_inv, beta_inv)),
        gtpow(pre.delta2R[k], beta_inv)});
    const GT& U = proof.batch_U[k];
    std::array<Bytes, 3> u_fields{u64be(t), u64be(half), serialize_gt(U)};
    out.transcript.absorb("BPVME/HP/VME/BATCH-U/v1", u_fields);
    const Scalar gamma =
        out.transcript.challenge("BPVME/HP/VME/DORY-GAMMA/v1", t);
    out.gamma[t] = gamma;
    out.aggregate.D0 = product({gtpow(bar.D0, fmul(gamma, gamma)),
                                gtpow(U, gamma), out.fresh[t].D0});
    out.aggregate.D1 = gtmul(gtpow(bar.D1, gamma), out.fresh[t].D1);
    out.aggregate.D2 = gtmul(gtpow(bar.D2, gamma), out.fresh[t].D2);
  }
  std::array<Bytes, 2> final_fields{serialize_g1(proof.phi_final),
                                    serialize(proof.theta_final)};
  out.transcript.absorb("BPVME/HP/VME/DORY-FINAL/v1", final_fields);
  out.epsilon =
      out.transcript.challenge("BPVME/HP/VME/DORY-EPSILON/v1", pp.log_dimension);
  if (stats) stats->recurrence_ms += ms(recurrence_start);
  return true;
}

bool verify_deferred(const VmePublicParams& pp, const VmePrecomputation& pre,
                     const Digest& context, const Group& X,
                     std::span<const Scalar> z, const VmeProof& proof,
                     bool inputs_validated, VmeVerificationStats* stats) {
  if (stats) *stats = {};
  const auto total_start = Clock::now();
  try {
    if ((!inputs_validated && !validate_vme_proof(pp, proof)) ||
        z.size() != pp.dimension || !valid_g2(X))
      return false;
    const Digest digest = statement_digest(pp, context, X, z);
    const auto transcript_start = Clock::now();
    VmeTranscript transcript(digest);
    const std::size_t D = pp.log_dimension;
    std::vector<Scalar> rho(D), beta(D), alpha(D), gamma(D + 1);
    for (std::size_t j = 0; j < D; ++j) {
      const VmeRexpClaims claims =
          j == 0 ? fixed_rexp_claim(pre) : proof.rexp_claims[j - 1];
      absorb_rexp(transcript, j, pp.dimension >> j, claims);
      rho[j] = transcript.challenge("BPVME/HP/VME/REXP-RHO/v1", j);
    }
    transcript.absorb("BPVME/HP/VME/REXP-R/v1", serialize_g1(proof.R));
    std::vector<Scalar> reversed_rho(D);
    for (std::size_t j = 0; j < D; ++j)
      reversed_rho[D - j - 1] = rho[j];
    const auto tensor = tensor_vector(reversed_rho);
    const Scalar q = inner(tensor, z);
    std::array<Bytes, 3> initial{u64be(D), u64be(pp.dimension), serialize(q)};
    transcript.absorb("BPVME/HP/VME/DORY-INITIAL/v1", initial);
    for (std::size_t t = 1; t <= D; ++t) {
      const std::size_t k = t - 1, m = pp.dimension >> k, half = m / 2;
      const auto& fold = proof.dory_folds[k];
      std::array<Bytes, 6> beta_fields{
          u64be(k), u64be(m), serialize_gt(fold.D1L),
          serialize_gt(fold.D1R), serialize_gt(fold.D2L),
          serialize_gt(fold.D2R)};
      transcript.absorb("BPVME/HP/VME/DORY-BETA-MESSAGE/v1", beta_fields);
      beta[k] = transcript.challenge("BPVME/HP/VME/DORY-BETA/v1", k);
      std::array<Bytes, 4> alpha_fields{
          u64be(k), u64be(m), serialize_gt(fold.W1), serialize_gt(fold.W2)};
      transcript.absorb("BPVME/HP/VME/DORY-ALPHA-MESSAGE/v1", alpha_fields);
      alpha[k] = transcript.challenge("BPVME/HP/VME/DORY-ALPHA/v1", k);
      std::array<Bytes, 3> u_fields{
          u64be(t), u64be(half), serialize_gt(proof.batch_U[k])};
      transcript.absorb("BPVME/HP/VME/BATCH-U/v1", u_fields);
      gamma[t] = transcript.challenge("BPVME/HP/VME/DORY-GAMMA/v1", t);
    }
    std::array<Bytes, 2> final_fields{serialize_g1(proof.phi_final),
                                      serialize(proof.theta_final)};
    transcript.absorb("BPVME/HP/VME/DORY-FINAL/v1", final_fields);
    const Scalar epsilon =
        transcript.challenge("BPVME/HP/VME/DORY-EPSILON/v1", D);
    const Scalar eta =
        transcript.challenge("BPVME/HP/VME/COMBINE-ETA/v1", D);
    if (stats) stats->transcript_ms = ms(transcript_start);

    const auto inversion_start = Clock::now();
    std::vector<Scalar> inverse_inputs;
    inverse_inputs.reserve(3 * D + 1);
    inverse_inputs.insert(inverse_inputs.end(), rho.begin(), rho.end());
    inverse_inputs.insert(inverse_inputs.end(), beta.begin(), beta.end());
    inverse_inputs.insert(inverse_inputs.end(), alpha.begin(), alpha.end());
    inverse_inputs.push_back(epsilon);
    const std::vector<Scalar> inverses = batch_invert(inverse_inputs);
    if (stats) stats->batch_inversion_ms = ms(inversion_start);
    const std::size_t beta_offset = D, alpha_offset = 2 * D,
                      epsilon_offset = 3 * D;

    const auto recurrence_start = Clock::now();
    std::vector<std::array<SymbolicExpression, 3>> fresh(D + 1);
    SymbolicExpression outer = gt_atom(GtAtomKind::PairingX, 0);
    for (std::size_t j = 0; j < D; ++j) {
      const SymbolicExpression E = j
          ? gt_atom(GtAtomKind::RexpE, j)
          : gt_atom(GtAtomKind::Delta1R, 0);
      const SymbolicExpression F = j
          ? gt_atom(GtAtomKind::RexpF, j)
          : gt_atom(GtAtomKind::Delta2R, 0);
      const SymbolicExpression TL = j
          ? gt_atom(GtAtomKind::RexpTL, j)
          : gt_atom(GtAtomKind::PairingX, 1);
      const SymbolicExpression TR = j
          ? gt_atom(GtAtomKind::RexpTR, j)
          : gt_atom(GtAtomKind::Delta1R, 0);
      fresh[j + 1][0] = outer;
      multiply_powered(fresh[j + 1][0], E, rho[j]);
      multiply_powered(fresh[j + 1][0], F, inverses[j]);
      fresh[j + 1][1] = TL;
      multiply_powered(fresh[j + 1][1], TR, rho[j]);
      fresh[j + 1][2] = gt_atom(GtAtomKind::PairingX, j + 1);
      multiply_powered(fresh[j + 1][2],
                       gt_atom(GtAtomKind::Delta2R, j), inverses[j]);
      outer = fresh[j + 1][1];
    }
    SymbolicExpression A0 =
        powered(gt_atom(GtAtomKind::PairingLLprime), q);
    SymbolicExpression A1 = pair_atom(PairAtomKind::InitialLX);
    SymbolicExpression A2 = pair_atom(PairAtomKind::InitialRLprime);
    for (std::size_t t = 1; t <= D; ++t) {
      const std::size_t k = t - 1;
      const Scalar& b = beta[k];
      const Scalar& b_inv = inverses[beta_offset + k];
      const Scalar& a = alpha[k];
      const Scalar& a_inv = inverses[alpha_offset + k];
      const Scalar& g = gamma[t];
      SymbolicExpression B0 = A0;
      multiply(B0, gt_atom(GtAtomKind::PairingX, k));
      multiply_powered(B0, A1, b_inv);
      multiply_powered(B0, A2, b);
      multiply_powered(B0, gt_atom(GtAtomKind::DoryW1, k), a);
      multiply_powered(B0, gt_atom(GtAtomKind::DoryW2, k), a_inv);
      SymbolicExpression B1 =
          powered(gt_atom(GtAtomKind::DoryD1L, k), a);
      multiply(B1, gt_atom(GtAtomKind::DoryD1R, k));
      multiply_powered(B1, gt_atom(GtAtomKind::PairingX, k + 1),
                       fmul(a, b));
      multiply_powered(B1, gt_atom(GtAtomKind::Delta1R, k), b);
      SymbolicExpression B2 =
          powered(gt_atom(GtAtomKind::DoryD2L, k), a_inv);
      multiply(B2, gt_atom(GtAtomKind::DoryD2R, k));
      multiply_powered(B2, gt_atom(GtAtomKind::PairingX, k + 1),
                       fmul(a_inv, b_inv));
      multiply_powered(B2, gt_atom(GtAtomKind::Delta2R, k), b_inv);
      A0 = powered(B0, fmul(g, g));
      multiply_powered(A0, gt_atom(GtAtomKind::BatchU, t), g);
      multiply(A0, fresh[t][0]);
      A1 = powered(B1, g);
      multiply(A1, fresh[t][1]);
      A2 = powered(B2, g);
      multiply(A2, fresh[t][2]);
      normalize(A0); normalize(A1); normalize(A2);
    }
    if (stats) stats->recurrence_ms = ms(recurrence_start);

    const auto terminal_start = Clock::now();
    const Scalar& epsilon_inv = inverses[epsilon_offset];
    const G1 terminal_g1 =
        g1add(proof.phi_final, g1mul(pp.auxiliary_G[0], epsilon));
    const Group terminal_g2 =
        g2add(proof.theta_final, g2mul(pp.fixed_P[0], epsilon_inv));
    SymbolicExpression dory = A0;
    multiply_powered(dory, A1, epsilon_inv);
    multiply_powered(dory, A2, epsilon);
    multiply(dory, gt_atom(GtAtomKind::PairingX, D));
    multiply_powered(dory, pair_atom(PairAtomKind::TerminalDory),
                     fneg(fone()));
    SymbolicExpression rexp = pair_atom(PairAtomKind::TerminalRexp);
    multiply_powered(rexp, outer, fneg(fone()));
    SymbolicExpression combined = dory;
    multiply_powered(combined, rexp, eta);
    normalize(combined);
    if (stats) {
      stats->terminal_ms = ms(terminal_start);
      stats->gt_terms = combined.gt.size();
      stats->pairing_terms = combined.pairs.size();
    }

    std::vector<GT> gt_bases;
    std::vector<Scalar> gt_exponents;
    gt_bases.reserve(combined.gt.size());
    gt_exponents.reserve(combined.gt.size());
    auto resolve_gt = [&](const GtAtomId& id) -> const GT& {
      switch (id.kind) {
        case GtAtomKind::PairingX: return pre.pairing_x[id.index];
        case GtAtomKind::Delta1R: return pre.delta1R[id.index];
        case GtAtomKind::Delta2R: return pre.delta2R[id.index];
        case GtAtomKind::PairingLLprime: return pre.pairing_LLprime;
        case GtAtomKind::RexpE: return proof.rexp_claims[id.index - 1].E;
        case GtAtomKind::RexpF: return proof.rexp_claims[id.index - 1].F;
        case GtAtomKind::RexpTL: return proof.rexp_claims[id.index - 1].TL;
        case GtAtomKind::RexpTR: return proof.rexp_claims[id.index - 1].TR;
        case GtAtomKind::DoryD1L: return proof.dory_folds[id.index].D1L;
        case GtAtomKind::DoryD1R: return proof.dory_folds[id.index].D1R;
        case GtAtomKind::DoryD2L: return proof.dory_folds[id.index].D2L;
        case GtAtomKind::DoryD2R: return proof.dory_folds[id.index].D2R;
        case GtAtomKind::DoryW1: return proof.dory_folds[id.index].W1;
        case GtAtomKind::DoryW2: return proof.dory_folds[id.index].W2;
        case GtAtomKind::BatchU: return proof.batch_U[id.index - 1];
      }
      throw std::logic_error("unknown GT atom");
    };
    for (const auto& [id, coefficient] : combined.gt) {
      gt_bases.push_back(resolve_gt(id));
      gt_exponents.push_back(coefficient);
    }
    const auto multiexp_start = Clock::now();
    const GT gt_value = gt_multiexp(gt_bases, gt_exponents);
    if (stats) stats->gt_multiexp_ms = ms(multiexp_start);

    auto resolve_pair = [&](const PairAtomId& id) -> PairInputs {
      switch (id.kind) {
        case PairAtomKind::InitialLX: return {pp.L, X};
        case PairAtomKind::InitialRLprime: return {proof.R, pp.Lprime};
        case PairAtomKind::TerminalDory: return {terminal_g1, terminal_g2};
        case PairAtomKind::TerminalRexp: return {proof.R, pp.fixed_P[0]};
      }
      throw std::logic_error("unknown pairing atom");
    };
    struct PairAccumulator { G1 first; Group second; };
    std::map<Bytes, PairAccumulator> pair_groups;
    const auto pairing_start = Clock::now();
    for (const auto& [id, coefficient] : combined.pairs) {
      const PairInputs inputs = resolve_pair(id);
      const G1 scaled = g1mul(inputs.first, coefficient);
      if (scaled.isZero()) continue;
      Bytes key = serialize(inputs.second);
      auto it = pair_groups.find(key);
      if (it == pair_groups.end())
        pair_groups.emplace(std::move(key),
                            PairAccumulator{scaled, inputs.second});
      else
        it->second.first = g1add(it->second.first, scaled);
    }
    std::vector<G1> left;
    std::vector<Group> right;
    for (const auto& [key, accumulator] : pair_groups)
      if (!accumulator.first.isZero()) {
        left.push_back(accumulator.first);
        right.push_back(accumulator.second);
      }
    const GT pair_value = pairing_product(left, right);
    if (stats) stats->multi_pairing_ms = ms(pairing_start);
    const GT combined_value = gtmul(gt_value, pair_value);
    GT one; one.setOne();
    const bool accepted = combined_value == one;
    if (stats) {
      stats->accepted = accepted;
      stats->total_ms = ms(total_start);
    }
    return accepted;
  } catch (...) {
    if (stats) stats->total_ms = ms(total_start);
    return false;
  }
}

bool verify_common(const VmePublicParams& pp, const VmePrecomputation& pre,
                   const Digest& context, const Group& X,
                   std::span<const Scalar> z, const VmeProof& proof,
                   bool optimized, bool inputs_validated,
                   VmeVerificationStats* stats) {
  if (stats) *stats = {};
  try {
    if ((!inputs_validated && !validate_vme_proof(pp, proof)) ||
        z.size() != pp.dimension || !valid_g2(X))
      return false;
    VmeStatement statement{context, X, std::vector<Scalar>(z.begin(), z.end()),
                           statement_digest(pp, context, X, z)};
    Replay replay_state(statement.digest);
    if (!replay(pp, pre, statement, proof, replay_state, stats)) return false;
    const auto terminal_start = Clock::now();
    const Scalar epsilon_inv = finv(replay_state.epsilon);
    const GT nonpair = product({replay_state.aggregate.D0,
        gtpow(replay_state.aggregate.D1, epsilon_inv),
        gtpow(replay_state.aggregate.D2, replay_state.epsilon),
        pre.pairing_x[pp.log_dimension]});
    const G1 terminal_g1 =
        g1add(proof.phi_final, g1mul(pp.auxiliary_G[0], replay_state.epsilon));
    const Group terminal_g2 =
        g2add(proof.theta_final, g2mul(pp.fixed_P[0], epsilon_inv));
    bool accepted = false;
    if (!optimized) {
      const GT dory_pair = pairing(terminal_g1, terminal_g2);
      const GT rexp_pair = pairing(proof.R, pp.fixed_P[0]);
      accepted = nonpair == dory_pair && rexp_pair == replay_state.outer;
    } else {
      const Scalar eta =
          replay_state.transcript.challenge("BPVME/HP/VME/COMBINE-ETA/v1",
                                            pp.log_dimension);
      GT gt_part = gtmul(nonpair, gtpow(gtinv(replay_state.outer), eta));
      Scalar minus_one;
      minus_one = -1;
      std::array<G1, 2> left{g1mul(terminal_g1, minus_one),
                             g1mul(proof.R, eta)};
      std::array<Group, 2> right{terminal_g2, pp.fixed_P[0]};
      const GT paired = pairing_product(left, right);
      GT combined = gtmul(gt_part, paired);
      GT one; one.setOne();
      accepted = combined == one;
    }
    if (stats) {
      stats->terminal_ms += ms(terminal_start);
      stats->accepted = accepted;
    }
    return accepted;
  } catch (...) { return false; }
}

}

Bytes serialize_g1(const G1& value) { return encode(value); }
Bytes serialize_gt(const GT& value) { return encode(value); }

std::size_t g1_bytes() {
  initialize();
  static const std::size_t size = [] {
    G1 p; mcl::bn::hashAndMapToG1(p, "hp-g1-size", 10);
    return serialize_g1(p).size();
  }();
  return size;
}

std::size_t gt_bytes() {
  initialize();
  static const std::size_t size = [] {
    G1 p; Group q;
    mcl::bn::hashAndMapToG1(p, "hp-gt-size-1", 12);
    mcl::bn::hashAndMapToG2(q, "hp-gt-size-2", 12);
    return serialize_gt(pairing(p, q)).size();
  }();
  return size;
}

bool deserialize_g1(std::span<const std::uint8_t> bytes, G1& out) {
  G1 value;
  if (bytes.size() != g1_bytes() ||
      value.deserialize(bytes.data(), bytes.size()) != bytes.size() ||
      !valid_g1(value) || serialize_g1(value) != Bytes(bytes.begin(), bytes.end()))
    return false;
  out = value;
  return true;
}

bool deserialize_gt(std::span<const std::uint8_t> bytes, GT& out) {
  GT value;
  if (bytes.size() != gt_bytes() ||
      value.deserialize(bytes.data(), bytes.size()) != bytes.size() ||
      !valid_gt(value) || serialize_gt(value) != Bytes(bytes.begin(), bytes.end()))
    return false;
  out = value;
  return true;
}

bool deserialize_gt_canonical_unchecked(
    std::span<const std::uint8_t> bytes, GT& out) {
  GT value;
  if (bytes.size() != gt_bytes() ||
      value.deserialize(bytes.data(), bytes.size()) != bytes.size() ||
      serialize_gt(value) != Bytes(bytes.begin(), bytes.end()))
    return false;
  out = value;
  return true;
}

VmePublicParams setup_vme(std::span<const Group> fixed_P,
                          std::span<const std::uint8_t> seed,
                          std::string_view transcript_domain) {
  initialize();
  VmePublicParams pp;
  pp.dimension = fixed_P.size();
  pp.log_dimension = exact_log2(pp.dimension);
  pp.fixed_P.assign(fixed_P.begin(), fixed_P.end());
  pp.transcript_domain = std::string(transcript_domain);
  pp.auxiliary_G.reserve(pp.dimension);
  for (std::size_t i = 0; i < pp.dimension; ++i)
    pp.auxiliary_G.push_back(derive_g1("BPVME/HP/VME/AUX-G1/v1", seed, i));
  pp.L = derive_g1("BPVME/HP/VME/L-G1/v1", seed, pp.dimension);
  pp.Lprime = derive_g2("BPVME/HP/VME/LPRIME-G2/v1", seed);
  pp.digest = crs_digest(pp);
  return pp;
}

VmePrecomputation precompute_vme(const VmePublicParams& pp) {
  if (!validate_vme_params(pp)) throw std::invalid_argument("invalid VME parameters");
  VmePrecomputation pre;
  pre.pairing_x.reserve(pp.log_dimension + 1);
  for (std::size_t k = 0; k <= pp.log_dimension; ++k) {
    const std::size_t m = pp.dimension >> k;
    pre.pairing_x.push_back(pairing_product(
        std::span(pp.auxiliary_G).first(m), std::span(pp.fixed_P).first(m)));
  }
  pre.delta1R.reserve(pp.log_dimension);
  pre.delta2R.reserve(pp.log_dimension);
  for (std::size_t k = 0; k < pp.log_dimension; ++k) {
    const std::size_t m = pp.dimension >> k, half = m / 2;
    pre.delta1R.push_back(pairing_product(
        std::span(pp.auxiliary_G).subspan(half, half),
        std::span(pp.fixed_P).first(half)));
    pre.delta2R.push_back(pairing_product(
        std::span(pp.auxiliary_G).first(half),
        std::span(pp.fixed_P).subspan(half, half)));
  }
  pre.pairing_LLprime = pairing(pp.L, pp.Lprime);
  pre.binding_digest = precomp_digest(pp, pre);
  return pre;
}

bool validate_vme_params(const VmePublicParams& pp) noexcept {
  try {
    if (!power_of_two(pp.dimension) || pp.transcript_domain.empty() ||
        pp.log_dimension != exact_log2(pp.dimension) ||
        pp.auxiliary_G.size() != pp.dimension ||
        pp.fixed_P.size() != pp.dimension)
      return false;
    for (const auto& p : pp.auxiliary_G) if (!valid_g1(p, true)) return false;
    for (const auto& p : pp.fixed_P) if (!valid_g2(p, true)) return false;
    return valid_g1(pp.L, true) && valid_g2(pp.Lprime, true) &&
           crs_digest(pp) == pp.digest;
  } catch (...) { return false; }
}

bool validate_vme_precomputation_binding(const VmePublicParams& pp,
                                         const VmePrecomputation& pre) noexcept {
  try {
    return pre.pairing_x.size() == pp.log_dimension + 1 &&
           pre.delta1R.size() == pp.log_dimension &&
           pre.delta2R.size() == pp.log_dimension &&
           precomp_digest(pp, pre) == pre.binding_digest;
  } catch (...) { return false; }
}

bool validate_vme_proof(const VmePublicParams& pp,
                        const VmeProof& proof) noexcept {
  try {
    if (proof.rexp_claims.size() + 1 != pp.log_dimension ||
        proof.dory_folds.size() != pp.log_dimension ||
        proof.batch_U.size() != pp.log_dimension ||
        !valid_g1(proof.R) || !valid_g1(proof.phi_final) ||
        !valid_g2(proof.theta_final))
      return false;
    for (const auto& c : proof.rexp_claims)
      if (!valid_gt(c.E) || !valid_gt(c.F) || !valid_gt(c.TL) || !valid_gt(c.TR))
        return false;
    for (const auto& f : proof.dory_folds)
      if (!valid_gt(f.D1L) || !valid_gt(f.D1R) || !valid_gt(f.D2L) ||
          !valid_gt(f.D2R) || !valid_gt(f.W1) || !valid_gt(f.W2))
        return false;
    for (const auto& u : proof.batch_U) if (!valid_gt(u)) return false;
    return true;
  } catch (...) { return false; }
}

bool validate_vme_proof_batched(const VmePublicParams& pp,
                                const VmeProof& proof) noexcept {
  try {
    if (proof.rexp_claims.size() + 1 != pp.log_dimension ||
        proof.dory_folds.size() != pp.log_dimension ||
        proof.batch_U.size() != pp.log_dimension ||
        !valid_g1(proof.R) || !valid_g1(proof.phi_final) ||
        !valid_g2(proof.theta_final))
      return false;
    std::vector<const GT*> fields;
    fields.reserve(4 * proof.rexp_claims.size() +
                   6 * proof.dory_folds.size() +
                   proof.batch_U.size());
    for (const auto& claims : proof.rexp_claims)
      for (const GT* field :
           {&claims.E, &claims.F, &claims.TL, &claims.TR})
        fields.push_back(field);
    for (const auto& fold : proof.dory_folds)
      for (const GT* field :
           {&fold.D1L, &fold.D1R, &fold.D2L, &fold.D2R,
            &fold.W1, &fold.W2})
        fields.push_back(field);
    for (const auto& U : proof.batch_U) fields.push_back(&U);
    Bytes seed;
    frame(seed, "BPVME/VME/GT-BATCH-MEMBERSHIP/v1");
    frame(seed, pp.digest);
    frame(seed, u64be(fields.size()));
    std::vector<GT> bases;
    bases.reserve(fields.size());
    for (const GT* field : fields) {
      const Bytes encoded = serialize_gt(*field);
      GT decoded;
      if (decoded.deserialize(encoded.data(), encoded.size()) !=
              encoded.size() ||
          decoded != *field)
        return false;
      frame(seed, encoded);
      bases.push_back(*field);
    }
    const Digest batch_seed = sha256(seed);
    std::vector<Scalar> coefficients;
    coefficients.reserve(fields.size());
    for (std::size_t i = 0; i < fields.size(); ++i) {
      Scalar coefficient;
      for (std::uint64_t counter = 0;; ++counter) {
        Bytes input;
        frame(input, "BPVME/VME/GT-BATCH-COEFFICIENT/v1");
        frame(input, batch_seed);
        frame(input, u64be(i));
        frame(input, u64be(counter));
        const Digest hash = sha256(input);
        coefficient.setBigEndianMod(hash.data(), hash.size());
        if (!coefficient.isZero()) break;
      }
      coefficients.push_back(coefficient);
    }
    const GT combined = gt_multiexp(bases, coefficients);
    mpz_class order;
    mcl::gmp::setStr(order, Scalar::getModulo(), 10);
    GT powered;
    GT::pow(powered, combined, order);
    GT one; one.setOne();
    return powered == one;
  } catch (...) { return false; }
}

VmeProof prove_vme(const VmePublicParams& pp, const VmePrecomputation& pre,
                   const Digest& context, const Group& X,
                   std::span<const Scalar> z) {
  if (!validate_vme_params(pp) || !validate_vme_precomputation_binding(pp, pre) ||
      z.size() != pp.dimension || !valid_g2(X))
    throw std::invalid_argument("invalid VME prove input");
  const Phase1 phase1 = prove_phase1(pp, pre, context, X, z);
  return prove_phase2(pp, pre, phase1);
}

bool verify_vme_reference(const VmePublicParams& pp,
                          const VmePrecomputation& pre, const Digest& context,
                          const Group& X, std::span<const Scalar> z,
                          const VmeProof& proof, VmeVerificationStats* stats) {
  if (!validate_vme_precomputation_binding(pp, pre)) return false;
  return verify_common(pp, pre, context, X, z, proof, false, false, stats);
}

bool verify_vme_optimized(const VmePublicParams& pp,
                          const VmePrecomputation& pre, const Digest& context,
                          const Group& X, std::span<const Scalar> z,
                          const VmeProof& proof, VmeVerificationStats* stats) {
  if (!validate_vme_precomputation_binding(pp, pre)) return false;
  return verify_deferred(pp, pre, context, X, z, proof, false, stats);
}

bool verify_vme_prevalidated_optimized(
    const VmePublicParams& pp, const VmePrecomputation& pre,
    const Digest& context, const Group& X, std::span<const Scalar> z,
    const VmeProof& proof, VmeVerificationStats* stats) {
  return verify_deferred(pp, pre, context, X, z, proof, true, stats);
}

}
