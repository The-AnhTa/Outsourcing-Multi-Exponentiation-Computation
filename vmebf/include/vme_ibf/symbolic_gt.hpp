#pragma once
#include "vme_ibf/gt_multiexp.hpp"
#include <functional>
#include <map>
namespace vme_ibf {
enum class GtAtomKind {PairingX,Delta1R,Delta2R,PairingLLprime,RexpE,RexpF,RexpTL,RexpTR,DoryD1L,DoryD1R,DoryD2L,DoryD2R,DoryW1,DoryW2,BatchU};
struct GtAtomId {GtAtomKind kind{};size_t index{};auto operator<=>(const GtAtomId&)const=default;};
enum class PairingAtomKind {InitialLX,InitialRLprime,TerminalDory,TerminalRexp};
struct PairingAtomId {PairingAtomKind kind{};auto operator<=>(const PairingAtomId&)const=default;};
struct SymbolicGtExpression {std::map<GtAtomId,Fr>gt_terms;std::map<PairingAtomId,Fr>pairing_terms;};
struct PairingInputs {G1 first;G2 second;};
enum class MultiexpMode {Reference,Pippenger};
using GtAtomResolver=std::function<const GT&(GtAtomId)>;
using PairingAtomResolver=std::function<PairingInputs(PairingAtomId)>;
struct PairingProductStats {size_t input_pairing_terms{},coalesced_pairing_terms{},miller_loop_batches{},miller_loop_terms{},final_exponentiations{},ordinary_pairing_calls{};double gt_msm_ms{},multi_pairing_ms{};};
SymbolicGtExpression gt_atom(GtAtomId);SymbolicGtExpression pairing_atom(PairingAtomId);
void multiply_in_place(SymbolicGtExpression&,const SymbolicGtExpression&);SymbolicGtExpression multiplied(const SymbolicGtExpression&,const SymbolicGtExpression&);
void power_in_place(SymbolicGtExpression&,const Fr&);SymbolicGtExpression powered(const SymbolicGtExpression&,const Fr&);void normalize(SymbolicGtExpression&);
GT evaluate_pairing_terms(const SymbolicGtExpression&,const PairingAtomResolver&,size_t* calls=nullptr);
GT evaluate_pairing_terms_reference(const SymbolicGtExpression&,const PairingAtomResolver&,PairingProductStats* stats=nullptr);
GT evaluate_pairing_terms_multi(const SymbolicGtExpression&,const PairingAtomResolver&,PairingProductStats* stats=nullptr);
GT evaluate_symbolic_expression(SymbolicGtExpression,const GtAtomResolver&,const PairingAtomResolver&,MultiexpMode,size_t* pairing_calls=nullptr,PairingProductStats* pairing_stats=nullptr);
}
