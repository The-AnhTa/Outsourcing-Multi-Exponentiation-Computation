#include "crypto.hpp"

#include <mcl/fp.hpp>

#include <limits>

namespace rexp::internal {

void append(Bytes& out, const void* data, std::size_t size) {
    const auto* first = static_cast<const std::uint8_t*>(data);
    out.insert(out.end(), first, first + size);
}

void append(Bytes& out, std::string_view text) {
    append(out, text.data(), text.size());
}

void append_u64(Bytes& out, std::uint64_t value) {
    for (int shift = 56; shift >= 0; shift -= 8) {
        out.push_back(static_cast<std::uint8_t>(value >> shift));
    }
}

void append_u64_be(Bytes& out, std::uint64_t value) {
    append_u64(out, value);
}

void frame(Bytes& out, const void* data, std::size_t size) {
    append_u64(out, static_cast<std::uint64_t>(size));
    append(out, data, size);
}

void frame(Bytes& out, std::string_view text) {
    frame(out, text.data(), text.size());
}

void append_frame(Bytes& out, const void* data, std::size_t size) {
    frame(out, data, size);
}

void append_frame(Bytes& out, std::string_view text) {
    frame(out, text);
}

void frame_digest(Bytes& out, const Digest& digest) {
    frame(out, digest.data(), digest.size());
}

Digest sha256(const Bytes& input) {
    Digest out{};
    const std::uint32_t written = mcl::fp::sha256(
        out.data(), static_cast<std::uint32_t>(out.size()),
        input.data(), static_cast<std::uint32_t>(input.size()));
    if (written != out.size()) throw std::runtime_error("SHA-256 failed");
    return out;
}

std::uint64_t Reader::readU64() {
    if (bytes_.size() - position_ < 8) {
        throw std::invalid_argument("truncated u64");
    }
    std::uint64_t value = 0;
    for (int i = 0; i < 8; ++i) {
        value = (value << 8) | bytes_[position_++];
    }
    return value;
}

std::size_t Reader::readSize() {
    const std::uint64_t value = readU64();
    if (value > std::numeric_limits<std::size_t>::max()) {
        throw std::invalid_argument("wire size does not fit size_t");
    }
    return static_cast<std::size_t>(value);
}

Bytes Reader::readFrame() {
    const std::uint64_t size = readU64();
    if (size > bytes_.size() - position_) {
        throw std::invalid_argument("truncated frame");
    }
    const std::size_t count = static_cast<std::size_t>(size);
    Bytes out(bytes_.begin() + position_, bytes_.begin() + position_ + count);
    position_ += count;
    return out;
}

void Reader::expect(std::string_view text) {
    const Bytes actual = readFrame();
    if (actual.size() != text.size()
        || !std::equal(actual.begin(), actual.end(), text.begin())) {
        throw std::invalid_argument("wire domain mismatch");
    }
}

void Reader::finish() const {
    if (position_ != bytes_.size()) {
        throw std::invalid_argument("trailing bytes");
    }
}

void Reader::requireFramedElementsFit(std::size_t count) const {
    // Every framed element needs at least its eight-byte length prefix. This
    // rejects impossible counts before reserve() can perform a large
    // attacker-controlled allocation.
    if (count > remaining() / 8) {
        throw std::invalid_argument("element count exceeds remaining wire data");
    }
}

void check_dimension(std::size_t d, std::size_t n, const char* protocol) {
    if (d >= std::numeric_limits<std::size_t>::digits) {
        throw std::invalid_argument(std::string(protocol) + " d is too large");
    }
    if (n != (std::size_t{1} << d)) {
        throw std::invalid_argument(std::string(protocol) + " n != 2^d");
    }
    if (n > (std::numeric_limits<std::size_t>::max() - 3U) / 4U) {
        throw std::invalid_argument(
            std::string(protocol) + " dimension overflows pairing count");
    }
}

GT pairing_product(
    const std::vector<G1>& left,
    std::size_t left_offset,
    const std::vector<G2>& right,
    std::size_t right_offset,
    std::size_t count) {
    if (left_offset > left.size() || right_offset > right.size()
        || count > left.size() - left_offset
        || count > right.size() - right_offset) {
        throw std::invalid_argument("pairing-product slice is out of range");
    }
    if (count == 0) {
        GT one;
        one.setOne();
        return one;
    }
    GT miller;
    mcl::bn::millerLoopVec(
        miller, left.data() + left_offset, right.data() + right_offset,
        count, true);
    GT out;
    mcl::bn::finalExp(out, miller);
    return out;
}

Fr inverse(const Fr& value) {
    if (value.isZero()) throw std::invalid_argument("cannot invert zero challenge");
    Fr out;
    Fr::inv(out, value);
    return out;
}

Fr fr_mul(const Fr& left, const Fr& right) {
    Fr out;
    Fr::mul(out, left, right);
    return out;
}

G1 g1_mul(const G1& point, const Fr& scalar) {
    G1 out;
    G1::mul(out, point, scalar);
    return out;
}

G2 g2_mul(const G2& point, const Fr& scalar) {
    G2 out;
    G2::mul(out, point, scalar);
    return out;
}

G1 g1_add(const G1& left, const G1& right) {
    G1 out;
    G1::add(out, left, right);
    return out;
}

G2 g2_add(const G2& left, const G2& right) {
    G2 out;
    G2::add(out, left, right);
    return out;
}

GT gt_mul(const GT& left, const GT& right) {
    GT out;
    GT::mul(out, left, right);
    return out;
}

GT gt_pow(const GT& value, const Fr& scalar) {
    GT out;
    GT::pow(out, value, scalar);
    return out;
}

} // namespace rexp::internal
