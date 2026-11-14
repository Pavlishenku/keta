// Result algebra beyond linear combination:
//  * envelope()            — pointwise min/max over every numeric field with
//    the dominant-source index behind every extremum (bit-copy of the winner,
//    first-index ties, non-finite inputs fail loud);
//  * envelope_factored()   — each source carries [c_min, c_max]; per field,
//    the extremes are max_i max(c_min v_i, c_max v_i) and the matching min:
//    the most unfavourable SINGLE factored source, not a sum of sources (for
//    a sum over admissible factors, see combination_extreme, keta/combination.h);
//  * envelope_concomitant()— pick a privileged scalar, return the ENTIRE
//    companion state of the winning source (design checks need consistent
//    internal-force vectors, never mixed extrema).
#pragma once
#include "keta/statics.h"

namespace keta {

struct StationSrc {
  std::array<unsigned, 6> forces{};  // N Vy Vz T My Mz
  std::array<unsigned, 3> defl{};
  std::array<unsigned, 3> rot{};
};
struct BarSrc {
  std::array<unsigned, 12> end{};
  std::vector<StationSrc> stations;
};
// Dominant-source index behind each enveloped facet-shell scalar. The TWELVE
// enveloped fields, in this order:
//   0 nx  1 ny  2 nxy  3 mx  4 my  5 mxy  6 qx  7 qy  8 m1  9 m2  10 n1  11 n2
// They are enveloped as INDEPENDENT scalars (m1's bound may come from a
// different case than mx's); the per-case principal ANGLES m_theta / n_theta
// are meaningless on a bound and are zeroed there.
using ShellSrc = std::array<unsigned, 12>;

struct EnvelopeBound {
  CaseResult value;
  std::vector<std::array<unsigned, 6>> disp_src;      // per node
  std::vector<std::array<unsigned, 6>> reaction_src;  // per reaction
  std::vector<BarSrc> bar_src;
  std::vector<ShellSrc> shell_src;                    // per facet shell
};
struct Envelope {
  EnvelopeBound lo, hi;
};

Envelope envelope(const std::vector<const CaseResult*>& sources);

struct FactoredSource {
  const CaseResult* result = nullptr;
  double c_min = 1.0, c_max = 1.0;
};
Envelope envelope_factored(const std::vector<FactoredSource>& sources);

struct Privileged {
  enum class What { StationForce, Reaction, Displacement };
  What what = What::StationForce;
  Id bar = 0;
  std::size_t station = 0;
  int comp = 0;
  Id node = 0;
};
struct Concomitant {
  double max = 0, min = 0;
  std::size_t src_max = 0, src_min = 0;
  const CaseResult* state_max = nullptr;  // full companion states
  const CaseResult* state_min = nullptr;
};
Concomitant envelope_concomitant(const std::vector<const CaseResult*>& sources,
                                 const Privileged& p);

}  // namespace keta
