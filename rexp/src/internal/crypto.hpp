#pragma once

#include "rexp/dory_setup.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace rexp::internal {

using Bytes = std::vector<std::uint8_t>;

void append(Bytes& out, const void* data, std::size_t size);
void append(Bytes& out, std::string_view text);
void append_u64(Bytes& out, std::uint64_t value);
void append_u64_be(Bytes& out, std::uint64_t value);
void frame(Bytes& out, const void* data, std::size_t size);
void frame(Bytes& out, std::string_view text);
void append_frame(Bytes& out, const void* data, std::size_t size);
void append_frame(Bytes& out, std::string_view text);
void frame_digest(Bytes& out, const Digest& digest);

template<class T>
Bytes encode(const T& value) {
    Bytes out(2048);
    const std::size_t written = value.serialize(out.data(), out.size());
    if (written == 0) throw std::runtime_error("mcl element serialization failed");
    out.resize(written);
    return out;
}

template<class T>
Bytes serialize_element(const T& value) {
    return encode(value);
}

template<class T>
void frame_element(Bytes& out, const T& value) {
    const Bytes encoded = encode(value);
    frame(out, encoded.data(), encoded.size());
}

template<class T>
void append_framed_element(Bytes& out, const T& value) {
    frame_element(out, value);
}

template<class T>
void append_fixed_element(
    Bytes& out,
    const T& value,
    std::size_t expected_size,
    const char* description) {
    const Bytes encoded = encode(value);
    if (encoded.size() != expected_size) {
        throw std::runtime_error(
            std::string("unexpected fixed-width ") + description + " encoding");
    }
    append(out, encoded.data(), encoded.size());
}

Digest sha256(const Bytes& input);

class Reader {
public:
    explicit Reader(const Bytes& bytes) : bytes_(bytes) {}

    std::uint64_t readU64();
    std::size_t readSize();
    Bytes readFrame();
    void expect(std::string_view text);

    template<class T>
    T element() {
        const Bytes bytes = readFrame();
        T value;
        const std::size_t used = value.deserialize(bytes.data(), bytes.size());
        if (used != bytes.size() || encode(value) != bytes) {
            throw std::invalid_argument("noncanonical element encoding");
        }
        return value;
    }

    void finish() const;
    std::size_t remaining() const { return bytes_.size() - position_; }
    void requireFramedElementsFit(std::size_t count) const;

private:
    const Bytes& bytes_;
    std::size_t position_ = 0;
};

void check_dimension(std::size_t d, std::size_t n, const char* protocol);

template<class Point>
void require_point(const Point& point, bool nonzero, const char* name) {
    if (!point.isValid() || !point.isValidOrder()
        || (nonzero && point.isZero())) {
        throw std::invalid_argument(std::string("invalid ") + name);
    }
}

GT pairing_product(
    const std::vector<G1>& left,
    std::size_t left_offset,
    const std::vector<G2>& right,
    std::size_t right_offset,
    std::size_t count);

Fr inverse(const Fr& value);
Fr fr_mul(const Fr& left, const Fr& right);
G1 g1_mul(const G1& point, const Fr& scalar);
G2 g2_mul(const G2& point, const Fr& scalar);
G1 g1_add(const G1& left, const G1& right);
G2 g2_add(const G2& left, const G2& right);
GT gt_mul(const GT& left, const GT& right);
GT gt_pow(const GT& value, const Fr& scalar);

} // namespace rexp::internal
