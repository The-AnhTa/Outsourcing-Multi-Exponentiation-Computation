#include "rexp/rexp.hpp"

#include <mcl/fp.hpp>
#include <algorithm>
#include <chrono>
#include <limits>
#include <stdexcept>

namespace rexp {
namespace {
using Bytes = std::vector<std::uint8_t>;
using Clock = std::chrono::steady_clock;

void raw(Bytes& o, const void* p, std::size_t n) {
    const auto* b = static_cast<const std::uint8_t*>(p);
    o.insert(o.end(), b, b + n);
}
void u64(Bytes& o, std::uint64_t x) {
    for (int s = 56; s >= 0; s -= 8)
        o.push_back(static_cast<std::uint8_t>(x >> s));
}
void frame(Bytes& o, const void* p, std::size_t n) {
    u64(o, n); raw(o, p, n);
}
void frame(Bytes& o, std::string_view s) { frame(o, s.data(), s.size()); }
void fd(Bytes& o, const Digest& d) { frame(o, d.data(), d.size()); }
template<class T> Bytes enc(const T& x) {
    Bytes b(2048);
    const std::size_t n = x.serialize(b.data(), b.size());
    if (!n) throw std::runtime_error("mcl serialization failed");
    b.resize(n);
    return b;
}
template<class T> void fe(Bytes& o, const T& x) {
    const Bytes b = enc(x); frame(o, b.data(), b.size());
}
Digest sha(const Bytes& b) {
    Digest d{};
    if (mcl::fp::sha256(
            d.data(), 32, b.data(), static_cast<std::uint32_t>(b.size())) != 32)
        throw std::runtime_error("SHA-256 failed");
    return d;
}
class Reader {
public:
    explicit Reader(const Bytes& bytes):bytes_(bytes){}
    std::uint64_t readU64(){
        if(bytes_.size()-pos_<8)throw std::invalid_argument("truncated u64");
        std::uint64_t x=0;for(int i=0;i<8;++i)x=(x<<8)|bytes_[pos_++];
        return x;
    }
    Bytes readFrame(){
        const auto n=readU64();
        if(n>bytes_.size()-pos_)throw std::invalid_argument("truncated frame");
        Bytes out(bytes_.begin()+pos_,bytes_.begin()+pos_+static_cast<std::size_t>(n));
        pos_+=static_cast<std::size_t>(n);return out;
    }
    void expect(std::string_view text){
        const Bytes b=readFrame();
        if(b.size()!=text.size()||!std::equal(b.begin(),b.end(),text.begin()))
            throw std::invalid_argument("wire domain mismatch");
    }
    template<class T>T element(){
        const Bytes b=readFrame();T x;
        const std::size_t used=x.deserialize(b.data(),b.size());
        if(used!=b.size()||enc(x)!=b)
            throw std::invalid_argument("noncanonical element encoding");
        return x;
    }
    void finish(){if(pos_!=bytes_.size())throw std::invalid_argument("trailing bytes");}
private:
    const Bytes& bytes_;std::size_t pos_=0;
};
void check_dimension(std::size_t d, std::size_t n) {
    if (d >= std::numeric_limits<std::size_t>::digits)
        throw std::invalid_argument("Rexp d is too large");
    if (n != (std::size_t{1} << d))
        throw std::invalid_argument("Rexp n != 2^d");
    if (n > (std::numeric_limits<std::size_t>::max() - 3U) / 4U)
        throw std::invalid_argument("Rexp dimension overflows pairing count");
}
template<class Point>
void require_point(const Point& p, bool nonzero, const char* name) {
    if (!p.isValid() || !p.isValidOrder() || (nonzero && p.isZero()))
        throw std::invalid_argument(std::string("invalid ") + name);
}
GT pp(const std::vector<G1>& a, std::size_t ao,
      const std::vector<G2>& b, std::size_t bo, std::size_t n) {
    if (ao > a.size() || bo > b.size()
        || n > a.size() - ao || n > b.size() - bo)
        throw std::invalid_argument("pairing-product slice out of range");
    GT m, r;
    mcl::bn::millerLoopVec(m, a.data() + ao, b.data() + bo, n, true);
    mcl::bn::finalExp(r, m);
    return r;
}
GT mul(const GT&a,const GT&b){GT r;GT::mul(r,a,b);return r;}
GT power(const GT&a,const Fr&s){GT r;GT::pow(r,a,s);return r;}
G1 smul(const G1&a,const Fr&s){G1 r;G1::mul(r,a,s);return r;}
G2 smul(const G2&a,const Fr&s){G2 r;G2::mul(r,a,s);return r;}
G1 add(const G1&a,const G1&b){G1 r;G1::add(r,a,b);return r;}
G2 add(const G2&a,const G2&b){G2 r;G2::add(r,a,b);return r;}
Fr inverse(const Fr& x) {
    if (x.isZero()) throw std::invalid_argument("zero challenge");
    Fr r; Fr::inv(r,x); return r;
}

Digest crs_digest(const RawRexpCRS& c) {
    Bytes b;
    frame(b,"REXP-G1-CRS-BN254-V1"); frame(b,"BN254");
    frame(b,"G1-G2-GT-PAIRING"); u64(b,c.d); u64(b,c.n);
    frame(b,"LEFT-PREFIX-CHAIN-V1");
    for(const auto&x:c.Gamma)fe(b,x);
    for(const auto&x:c.Lambda)fe(b,x);
    return sha(b);
}
Digest statement_digest(
    const PreparedPublicParameters& p,
    const std::vector<G1>& H,
    const GT& d1,const GT& e,const GT& f,const GT& tl,const GT& tr) {
    Bytes b;
    frame(b,"REXP-G1-STATEMENT-BN254-V1"); fd(b,p.digest());
    u64(b,p.d()); u64(b,p.n());
    for(const auto&x:H)fe(b,x);
    fe(b,d1);fe(b,e);fe(b,f);fe(b,tl);fe(b,tr);
    return sha(b);
}
Digest initial(const PreparedPublicParameters&p,const PreparedStatement&s){
    Bytes b;frame(b,"REXP-G1-FS-v1");frame(b,"BN254");frame(b,"G1");
    fd(b,p.digest());fd(b,s.digest());
    u64(b,p.d());u64(b,p.n());return sha(b);
}
Digest round_digest(const Digest&t,std::size_t k,std::size_t m,
                    const RexpRoundMessage&r){
    Bytes b;frame(b,"REXP-G1-ROUND-MESSAGE-V1");fd(b,t);u64(b,k);u64(b,m);
    fe(b,r.E);fe(b,r.F);fe(b,r.TL);fe(b,r.TR);return sha(b);
}
Digest enter(const Digest&t,std::size_t k,std::size_t h,
             const DoryStatement&s){
    Bytes b;frame(b,"REXP-G1-EMBEDDED-DORY-ENTER-V1");fd(b,t);u64(b,k);
    u64(b,k+1);u64(b,h);fe(b,s.D0);fe(b,s.D1);fe(b,s.D2);return sha(b);
}
Digest leave(const Digest&t,std::size_t k){
    Bytes b;frame(b,"REXP-G1-EMBEDDED-DORY-EXIT-V1");fd(b,t);u64(b,k);
    return sha(b);
}
Digest final_digest(const Digest&t,const G1&r){
    Bytes b;frame(b,"REXP-G1-FINAL-OUTPUT-V1");fd(b,t);fe(b,r);return sha(b);
}
RexpRoundMessage msg0(const PreparedStatement&s){
    return{s.E0(),s.F0(),s.TL0(),s.TR0()};
}
DoryCRS make_level_crs(const PreparedPublicParameters&p,std::size_t level){
    if (level > p.d()) throw std::invalid_argument("invalid generator level");
    DoryCRS v;v.d=p.d()-level;v.n=p.n()>>level;
    v.Gamma.assign(p.Gamma().begin(),p.Gamma().begin()+v.n);
    v.Lambda.assign(p.Lambda().begin(),p.Lambda().begin()+v.n);
    v.digest=ComputeDoryCRSDigest(v);return v;
}
DoryPrecomputation make_level_pc(
    const PreparedPublicParameters&p,std::size_t level,const Digest& digest){
    if (level > p.d()) throw std::invalid_argument("invalid precomp level");
    DoryPrecomputation v;
    v.X.assign(p.X().begin()+level,p.X().end());
    v.Delta1R.assign(p.Delta1R().begin()+level,p.Delta1R().end());
    v.Delta2R.assign(p.Delta2R().begin()+level,p.Delta2R().end());
    v.crs_digest=digest;
    v.pairing_product_terms=4*(p.n()>>level)-3;
    return v;
}
void require_bound(
    const PreparedPublicParameters&p,const PreparedStatement&s) {
    if (s.crsDigest()!=p.digest() || s.H().size()!=p.n())
        throw std::invalid_argument("prepared statement belongs to another CRS");
}
bool proof_shape(const PreparedPublicParameters&p,const RexpProof&proof) {
    if (proof.doryProofs.size()!=p.d()
        || proof.dynamicRoundMessages.size()!=(p.d()?p.d()-1:0)
        || !proof.R.isValid() || !proof.R.isValidOrder()) return false;
    for(std::size_t k=0;k<p.d();++k) {
        if(proof.doryProofs[k].rounds.size()!=p.d()-k-1) return false;
    }
    return true;
}
bool validate_proof_gt(
    const RexpProof& proof, RexpProofValidationMetrics* metrics) {
    const auto start=metrics?Clock::now():Clock::time_point{};
    std::size_t checked=0;
    auto check=[&](const GT& x){
        ++checked;
        if(!mcl::bn::isValidGT(x))
            throw std::invalid_argument("proof GT is outside target subgroup");
    };
    for(const auto& message:proof.dynamicRoundMessages){
        check(message.E);check(message.F);check(message.TL);check(message.TR);
    }
    for(const auto& dory:proof.doryProofs){
        for(const auto& round:dory.rounds){
            check(round.D1L);check(round.D1R);check(round.D2L);
            check(round.D2R);check(round.W1);check(round.W2);
        }
    }
    if(metrics){
        metrics->gt_elements_checked=checked;
        metrics->gt_subgroup_validation_ms=
            std::chrono::duration<double,std::milli>(Clock::now()-start).count();
    }
    return true;
}
bool verify_core(
    const PreparedPublicParameters&p,const PreparedStatement&s,
    const RexpProof&proof,bool reference,bool already_validated,
    RexpVerifyMetrics*out) {
    const auto start=out?Clock::now():Clock::time_point{};
    if(out)*out=RexpVerifyMetrics{};
    try {
        require_bound(p,s);
        if(!proof_shape(p,proof))return false;
        if(!already_validated&&!validate_proof_gt(proof,nullptr))return false;
        Digest t=initial(p,s);
        GT d1=s.D1Initial();
#ifdef REXP_ENABLE_PROFILING
        if(out)out->per_dory.reserve(p.d());
#endif
        for(std::size_t k=0;k<p.d();++k){
            const std::size_t m=p.n()>>k,h=m/2;
            const RexpRoundMessage r=k?proof.dynamicRoundMessages[k-1]:msg0(s);
            const Digest tr=round_digest(t,k,m,r);
            const Fr rho=ChallengeNonzeroFr(tr,"REXP-G1-RHO-V1",k);
            const Fr ri=inverse(rho);
            const DoryStatement ds{
                mul(mul(d1,power(r.E,rho)),power(r.F,ri)),
                mul(r.TL,power(r.TR,rho)),
                mul(p.X()[k+1],power(p.Delta2R()[k],ri))};
            const Digest din=enter(tr,k,h,ds);
            Digest dend;
            if(reference) {
                if(!VerifyEmbeddedReference(
                    p.levelCRS(k+1),p.levelPrecomputation(k+1),ds,
                    proof.doryProofs[k],din,&dend))return false;
            } else {
                VerifyMetrics dm;
                const auto& level_crs=p.levelCRS(k+1);
                const auto& level_pc=p.levelPrecomputation(k+1);
                if(!VerifyEmbeddedDeferred(
                    level_crs,level_pc,ds,
                    proof.doryProofs[k],din,&dend,out?&dm:nullptr))return false;
                if(out){
                    ++out->dory_verifications;
                    out->gt_multiexponentiations+=dm.gt_multiexp_calls;
                    out->dory_terminal_pairings+=dm.terminal_pairings;
                    out->actual_gt_bases+=dm.actual_gt_bases;
                    out->transcript_ms+=dm.transcript_ms;
                    out->gt_multiexp_ms+=dm.gt_multiexp_ms;
                    out->dory_pairing_ms+=dm.terminal_pairing_ms;
#ifdef REXP_ENABLE_PROFILING
                    out->per_dory.push_back({
                        k, level_crs.d, dm.symbolic_atom_insertions,
                        dm.actual_gt_bases, dm.coalesced_duplicate_bases,
                        dm.zero_coefficients_removed,
                        dm.identity_bases_removed, dm.gt_multiexp_calls,
                        dm.terminal_pairings});
#endif
                }
            }
            t=leave(dend,k);
            d1=ds.D1;
        }
        (void)final_digest(t,proof.R);
        const auto ps=out?Clock::now():Clock::time_point{};
        if(p.Lambda()[0].isZero())return false;
        GT q;mcl::bn::pairing(q,proof.R,p.Lambda()[0]);
        if(out){out->final_pairings=1;out->final_pairing_ms=
            std::chrono::duration<double,std::milli>(Clock::now()-ps).count();}
        const bool ok=q==d1;
        if(out)out->total_ms=
            std::chrono::duration<double,std::milli>(Clock::now()-start).count();
        return ok;
    } catch (...) {
        return false;
    }
}
}

RawRexpCRS GenerateRawCRS(std::size_t d,std::string_view seed) {


    const DoryCRS c=GenerateDoryCRS(d,seed);
    return{c.d,c.n,c.Gamma,c.Lambda};
}

PreparedPublicParameters PreparePublicParameters(const RawRexpCRS& rawCrs) {
    initialize();
    check_dimension(rawCrs.d,rawCrs.n);
    if(rawCrs.Gamma.size()!=rawCrs.n||rawCrs.Lambda.size()!=rawCrs.n)
        throw std::invalid_argument("raw CRS vector length mismatch");
    for(const auto&x:rawCrs.Gamma)require_point(x,true,"CRS G1 point");
    for(const auto&x:rawCrs.Lambda)require_point(x,true,"CRS G2 point");
    PreparedPublicParameters out;
    out.d_=rawCrs.d;out.n_=rawCrs.n;
    out.Gamma_=rawCrs.Gamma;out.Lambda_=rawCrs.Lambda;
    out.digest_=crs_digest(rawCrs);
    out.X_.reserve(out.d_+1);
    for(std::size_t k=0;k<=out.d_;++k){
        const std::size_t m=out.n_>>k;
        out.X_.push_back(pp(out.Gamma_,0,out.Lambda_,0,m));
        out.pairingProductTerms_+=m;
    }
    out.Delta1R_.reserve(out.d_);out.Delta2R_.reserve(out.d_);
    for(std::size_t k=0;k<out.d_;++k){
        const std::size_t m=out.n_>>k,h=m/2;
        out.Delta1R_.push_back(pp(out.Gamma_,h,out.Lambda_,0,h));
        out.Delta2R_.push_back(pp(out.Gamma_,0,out.Lambda_,h,h));
        out.pairingProductTerms_+=2*h;
    }
    out.levelCRS_.resize(out.d_+1);
    out.levelPrecomputation_.resize(out.d_+1);


    for(std::size_t level=1;level<=out.d_;++level){
        out.levelCRS_[level]=make_level_crs(out,level);
        out.levelPrecomputation_[level]=
            make_level_pc(out,level,out.levelCRS_[level].digest);
    }
    return out;
}

RawRexpStatement GenerateRawStatement(const PreparedPublicParameters&p) {
    initialize();
    RawRexpStatement out;
    out.H.reserve(p.n());
    G1 base;
    static constexpr char domain[]="REXP-G1-PUBLIC-H-BASE-V1";
    mcl::bn::hashAndMapToG1(base,domain,sizeof(domain)-1);
    require_point(base,true,"statement base");
    for(std::size_t i=0;i<p.n();++i){
        Fr scalar;do{scalar.setByCSPRNG();}while(scalar.isZero());
        G1 point;G1::mul(point,base,scalar);out.H.push_back(point);
    }
    return out;
}

PreparedStatement PrepareStatement(
    const PreparedPublicParameters&p,const RawRexpStatement&rawStatement) {
    if(rawStatement.H.size()!=p.n())
        throw std::invalid_argument("raw statement length mismatch");
    for(const auto&x:rawStatement.H)
        require_point(x,false,"statement G1 point");
    PreparedStatement out;
    out.H_=rawStatement.H;out.crsDigest_=p.digest();
    out.D1Initial_=pp(out.H_,0,p.Lambda(),0,p.n());
    if(p.d()>0){
        const std::size_t h=p.n()/2;
        out.E0_=pp(out.H_,h,p.Lambda(),0,h);
        out.F0_=pp(out.H_,0,p.Lambda(),h,h);
        out.TL0_=pp(out.H_,0,p.Lambda(),0,h);
        out.TR0_=out.E0_;
    }else{
        out.E0_.setOne();out.F0_.setOne();
        out.TL0_=out.D1Initial_;out.TR0_.setOne();
    }
    out.digest_=statement_digest(
        p,out.H_,out.D1Initial_,out.E0_,out.F0_,out.TL0_,out.TR0_);
    return out;
}

RexpSetupResult Setup(std::size_t d,std::string_view seed){
    RexpSetupResult out;
    out.rawCRS=GenerateRawCRS(d,seed);
    out.params=PreparePublicParameters(out.rawCRS);
    out.rawStatement=GenerateRawStatement(out.params);
    out.statement=PrepareStatement(out.params,out.rawStatement);
    out.proverInput.H=out.rawStatement.H;
    return out;
}

RexpProof Prove(const PreparedPublicParameters&p,const PreparedStatement&s,
                const RexpProverInput&in){
    require_bound(p,s);
    if(in.H!=s.H())throw std::invalid_argument("prover input differs from H");
    RexpProof proof;Digest t=initial(p,s);std::vector<G1> H=s.H();
    GT d1=s.D1Initial();
    for(std::size_t k=0;k<p.d();++k){
        const std::size_t m=p.n()>>k,h=m/2;RexpRoundMessage r;
        if(k==0)r=msg0(s);
        else{
            r={pp(H,h,p.Lambda(),0,h),pp(H,0,p.Lambda(),h,h),
               pp(H,0,p.Lambda(),0,h),pp(H,h,p.Lambda(),0,h)};
            proof.dynamicRoundMessages.push_back(r);
        }
        const Digest tr=round_digest(t,k,m,r);
        const Fr rho=ChallengeNonzeroFr(tr,"REXP-G1-RHO-V1",k),ri=inverse(rho);
        std::vector<G1> next(h);std::vector<G2> theta(h);
        for(std::size_t i=0;i<h;++i){
            next[i]=add(H[i],smul(H[h+i],rho));
            theta[i]=add(p.Lambda()[i],smul(p.Lambda()[h+i],ri));
        }
        const DoryStatement ds{
            mul(mul(d1,power(r.E,rho)),power(r.F,ri)),
            mul(r.TL,power(r.TR,rho)),
            mul(p.X()[k+1],power(p.Delta2R()[k],ri))};
        const Digest din=enter(tr,k,h,ds);Digest dend;
        proof.doryProofs.push_back(
            ProveEmbedded(p.levelCRS(k+1),{next,theta},din,&dend));
        t=leave(dend,k);d1=ds.D1;H=std::move(next);
    }
    proof.R=H[0];(void)final_digest(t,proof.R);return proof;
}

bool VerifyPrepared(
    const PreparedPublicParameters&p,const PreparedStatement&s,
    const RexpProof&proof,RexpVerifyMetrics*metrics){
    return verify_core(p,s,proof,false,false,metrics);
}
bool VerifyValidatedProof(
    const PreparedPublicParameters&p,const PreparedStatement&s,
    const ValidatedRexpProof&proof,RexpVerifyMetrics*metrics){
    if(proof.d()!=p.d())return false;
    return verify_core(p,s,proof.proof(),false,true,metrics);
}
bool VerifyOptimized(
    const PreparedPublicParameters&p,const PreparedStatement&s,
    const ValidatedRexpProof&proof){
    return VerifyValidatedProof(p,s,proof);
}
bool VerifyReference(
    const PreparedPublicParameters&p,const PreparedStatement&s,
    const RexpProof&proof){
    return verify_core(p,s,proof,true,false,nullptr);
}
bool Verify(const RawRexpCRS&c,const RawRexpStatement&s,
            const RexpProof&proof,RexpVerifyMetrics*metrics){
    try{
        const PreparedPublicParameters p=PreparePublicParameters(c);
        const PreparedStatement prepared=PrepareStatement(p,s);
        return VerifyPrepared(p,prepared,proof,metrics);
    }catch(...){return false;}
}

std::vector<std::uint8_t> SerializeRexpCRS(const RawRexpCRS&c){
    check_dimension(c.d,c.n);
    Bytes b;u64(b,c.d);u64(b,c.n);
    for(const auto&x:c.Gamma)fe(b,x);for(const auto&x:c.Lambda)fe(b,x);
    const Digest d=crs_digest(c);raw(b,d.data(),d.size());return b;
}
std::vector<std::uint8_t> SerializeRexpPrecomputation(
    const PreparedPublicParameters&p){
    Bytes b;for(const auto&x:p.X())fe(b,x);
    for(const auto&x:p.Delta1R())fe(b,x);
    for(const auto&x:p.Delta2R())fe(b,x);return b;
}
std::vector<std::uint8_t> SerializeRexpStatement(const RawRexpStatement&s){
    Bytes b;for(const auto&x:s.H)fe(b,x);return b;
}
std::vector<std::uint8_t> SerializePreparedStatement(const PreparedStatement&s){
    Bytes b;for(const auto&x:s.H())fe(b,x);fe(b,s.D1Initial());fe(b,s.E0());
    fe(b,s.F0());fe(b,s.TL0());fe(b,s.TR0());raw(b,s.digest().data(),32);
    return b;
}
std::vector<std::uint8_t> SerializeRexpProof(
    const RexpProof&p,std::size_t d){
    if(d>=std::numeric_limits<std::size_t>::digits
       ||p.doryProofs.size()!=d
       ||p.dynamicRoundMessages.size()!=(d?d-1:0))
        throw std::invalid_argument("wrong proof shape");
    Bytes b;
    for(std::size_t k=0;k<d;++k){
        if(p.doryProofs[k].rounds.size()!=d-k-1)
            throw std::invalid_argument("wrong embedded Dory depth");
        if(k){const auto&r=p.dynamicRoundMessages[k-1];
            fe(b,r.E);fe(b,r.F);fe(b,r.TL);fe(b,r.TR);}
        for(const auto&r:p.doryProofs[k].rounds){
            fe(b,r.D1L);fe(b,r.D1R);fe(b,r.D2L);
            fe(b,r.D2R);fe(b,r.W1);fe(b,r.W2);}
        fe(b,p.doryProofs[k].PhiFinal);fe(b,p.doryProofs[k].ThetaFinal);
    }
    fe(b,p.R);return b;
}

std::vector<std::uint8_t> SerializeRexpCRSWire(const RawRexpCRS&c){
    check_dimension(c.d,c.n);
    Bytes b;frame(b,"REXP-G1-CRS-WIRE-BN254-V1");u64(b,c.d);u64(b,c.n);
    u64(b,c.Gamma.size());for(const auto&x:c.Gamma)fe(b,x);
    u64(b,c.Lambda.size());for(const auto&x:c.Lambda)fe(b,x);
    return b;
}
std::vector<std::uint8_t> SerializeRexpStatementWire(
    const RawRexpStatement&s,std::size_t n){
    if(s.H.size()!=n)throw std::invalid_argument("statement wire length mismatch");
    Bytes b;frame(b,"REXP-G1-RAW-STATEMENT-WIRE-BN254-V1");u64(b,n);u64(b,s.H.size());
    for(const auto&x:s.H)fe(b,x);return b;
}
std::vector<std::uint8_t> SerializeRexpProofWire(
    const RexpProof&p,std::size_t d){
    const Bytes payload=SerializeRexpProof(p,d);
    Bytes b;frame(b,"REXP-G1-PROOF-WIRE-BN254-V1");u64(b,d);
    u64(b,p.dynamicRoundMessages.size());u64(b,p.doryProofs.size());
    std::size_t dynamicIndex=0;
    for(std::size_t k=0;k<d;++k){
        if(k){const auto&r=p.dynamicRoundMessages[dynamicIndex++];
            fe(b,r.E);fe(b,r.F);fe(b,r.TL);fe(b,r.TR);}
        u64(b,p.doryProofs[k].rounds.size());
        for(const auto&r:p.doryProofs[k].rounds){
            fe(b,r.D1L);fe(b,r.D1R);fe(b,r.D2L);fe(b,r.D2R);fe(b,r.W1);fe(b,r.W2);}
        fe(b,p.doryProofs[k].PhiFinal);fe(b,p.doryProofs[k].ThetaFinal);
    }
    fe(b,p.R);return b;
}
RawRexpCRS DeserializeRexpCRSWire(const Bytes&bytes){
    Reader r(bytes);r.expect("REXP-G1-CRS-WIRE-BN254-V1");
    RawRexpCRS c;c.d=static_cast<std::size_t>(r.readU64());
    c.n=static_cast<std::size_t>(r.readU64());check_dimension(c.d,c.n);
    const auto ng=r.readU64();if(ng!=c.n)throw std::invalid_argument("G1 count mismatch");
    c.Gamma.reserve(c.n);for(std::size_t i=0;i<c.n;++i)c.Gamma.push_back(r.element<G1>());
    const auto nl=r.readU64();if(nl!=c.n)throw std::invalid_argument("G2 count mismatch");
    c.Lambda.reserve(c.n);for(std::size_t i=0;i<c.n;++i)c.Lambda.push_back(r.element<G2>());
    r.finish();(void)PreparePublicParameters(c);return c;
}
RawRexpStatement DeserializeRexpStatementWire(
    const Bytes&bytes,std::size_t expected_n){
    Reader r(bytes);r.expect("REXP-G1-RAW-STATEMENT-WIRE-BN254-V1");
    if(r.readU64()!=expected_n||r.readU64()!=expected_n)
        throw std::invalid_argument("statement count mismatch");
    RawRexpStatement s;s.H.reserve(expected_n);
    for(std::size_t i=0;i<expected_n;++i)s.H.push_back(r.element<G1>());
    r.finish();for(const auto&x:s.H)require_point(x,false,"statement G1 point");
    return s;
}
RexpProof DeserializeRexpProofWire(
    const Bytes&bytes,std::size_t expected_d,
    RexpProofValidationMetrics*metrics){
    const auto decodeStart=Clock::now();
    Reader r(bytes);r.expect("REXP-G1-PROOF-WIRE-BN254-V1");
    if(r.readU64()!=expected_d)throw std::invalid_argument("proof d mismatch");
    const std::size_t dynamic=expected_d?expected_d-1:0;
    if(r.readU64()!=dynamic||r.readU64()!=expected_d)
        throw std::invalid_argument("proof count mismatch");
    RexpProof p;p.dynamicRoundMessages.reserve(dynamic);p.doryProofs.reserve(expected_d);
    for(std::size_t k=0;k<expected_d;++k){
        if(k)p.dynamicRoundMessages.push_back(
            {r.element<GT>(),r.element<GT>(),r.element<GT>(),r.element<GT>()});
        const std::size_t depth=expected_d-k-1;
        if(r.readU64()!=depth)throw std::invalid_argument("Dory depth mismatch");
        DoryProof q;q.rounds.reserve(depth);
        for(std::size_t j=0;j<depth;++j)q.rounds.push_back(
            {r.element<GT>(),r.element<GT>(),r.element<GT>(),
             r.element<GT>(),r.element<GT>(),r.element<GT>()});
        q.PhiFinal=r.element<G1>();q.ThetaFinal=r.element<G2>();
        p.doryProofs.push_back(std::move(q));
    }
    p.R=r.element<G1>();r.finish();
    if(expected_d>=std::numeric_limits<std::size_t>::digits)
        throw std::invalid_argument("proof d too large");
    if(!p.R.isValid()||!p.R.isValidOrder())
        throw std::invalid_argument("invalid proof R G1 point");
    for(const auto& q:p.doryProofs){
        if(!q.PhiFinal.isValid()||!q.PhiFinal.isValidOrder()
           ||!q.ThetaFinal.isValid()||!q.ThetaFinal.isValidOrder())
            throw std::invalid_argument("invalid embedded Dory final point");
    }
    RexpProofValidationMetrics local;
    validate_proof_gt(p,&local);
    if(metrics){
        *metrics=local;
        metrics->canonical_decode_ms=
            std::chrono::duration<double,std::milli>(
                Clock::now()-decodeStart).count()
            - local.gt_subgroup_validation_ms;
    }
    return p;
}
ValidatedRexpProof DeserializeValidatedRexpProofWire(
    const Bytes&bytes,std::size_t d,RexpProofValidationMetrics*metrics){
    RexpProof proof=DeserializeRexpProofWire(bytes,d,metrics);
    ValidatedRexpProof out;out.proof_=std::move(proof);out.d_=d;return out;
}
bool IsValidGTSubgroup(const GT& x) {
    return mcl::bn::isValidGT(x);
}
ValidatedRexpProof ValidateRexpProof(
    const RexpProof&proof,std::size_t d,RexpProofValidationMetrics*metrics){
    if(d>=std::numeric_limits<std::size_t>::digits
       ||proof.doryProofs.size()!=d
       ||proof.dynamicRoundMessages.size()!=(d?d-1:0)
       ||!proof.R.isValid()||!proof.R.isValidOrder())
        throw std::invalid_argument("invalid proof shape");
    for(std::size_t k=0;k<d;++k){
        if(proof.doryProofs[k].rounds.size()!=d-k-1
           ||!proof.doryProofs[k].PhiFinal.isValid()
           ||!proof.doryProofs[k].PhiFinal.isValidOrder()
           ||!proof.doryProofs[k].ThetaFinal.isValid()
           ||!proof.doryProofs[k].ThetaFinal.isValidOrder())
            throw std::invalid_argument("invalid embedded Dory proof");
    }
    validate_proof_gt(proof,metrics);
    ValidatedRexpProof out;out.proof_=proof;out.d_=d;return out;
}
}
