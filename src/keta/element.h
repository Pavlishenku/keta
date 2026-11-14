// Bar element kernels: 3-D Timoshenko frame (12 DOF), truss, with
//   * shear deformation (phi = 12EI/(G Av L^2); Av = 0 -> Euler-Bernoulli),
//   * end releases by static condensation (with recovery re-injection),
//   * rigid-arm end eccentricity T_ecc (composite girder offsets),
//   * Winkler foundation bed (consistent NtN matrices, stiffness only),
//   * consistent + lumped mass, geometric stiffness Kg(N),
//   * fixed-end forces for UDL / linear / point / thermal / self-weight loads
//     (exact Timoshenko shape-function work integrals via Gauss quadrature).
//
// Local DOF order: [uxi uyi uzi rxi ryi rzi uxj uyj uzj rxj ryj rzj].
// Member-force sign convention on recovery: N > 0 tension; diagrams report
// internal forces at stations with the classical convention (M sagging positive
// about the corresponding local axis).
#pragma once
#include "keta/model.h"

namespace keta {

// Rotation: rows of R are the local axes expressed in the global frame.
// For non-vertical bars local z lies in the vertical plane of the member
// (tilted towards +Z); vertical bars take local y = global Y. `roll` then
// rotates y/z about local x.
Mat3 bar_rotation(const Vec3& xi, const Vec3& xj, double roll);

class BarElement {
 public:
  BarElement(const Model& model, const Bar& bar);

  const Bar& bar() const { return *bar_; }
  double L() const { return L_; }
  const Mat3& R() const { return R_; }  // local <- global

  // Nodal-frame (global, at the model nodes, ecc folded) matrices.
  const Mat12& ke() const { return ke_nodal_; }
  Mat12 me(bool lumped = false) const;
  Mat12 kg(double N) const;  // geometric stiffness for axial force N (tension +)

  // Equivalent nodal load vector (global nodal frame) for one span load.
  // Also accumulates the local fixed-end force correction used by recovery.
  Vec12 equivalent_nodal(const Load& load, const Material& mat) const;

  // Recovery: given global nodal displacements (12), return local end
  // displacements (releases re-injected) and local member end forces
  // f_local = k_local u_local - f0_local (f0 = sum of span-load FEFs).
  Vec12 local_displacements(const Vec12& u_nodal, const Vec12& f0_local) const;
  Vec12 end_forces_local(const Vec12& u_local, const Vec12& f0_local) const;

  // Fixed-end (span-load) local force vector for a load — needed by recovery.
  Vec12 f0_local(const Load& load, const Material& mat) const;

  // EXACT beam-column tangent (stability functions) in the nodal frame:
  // T^t (condense(ke_exact_local(N)) + soil) T. N tension-positive.
  Mat12 ke_exact(double N) const;
  // The raw local block behind ke_exact: uncondensed, no Winkler bed, no
  // transform (frame bars; a truss/cable returns its linear local block).
  Mat12 k_exact_local(double N) const;
  // Release handling with a GIVEN raw local matrix (the exact beam-column
  // block of a second-order solve instead of the linear one): span-load
  // condensation f0_k - k_kr k_rr^-1 f0_r, and the re-injection of the
  // released end displacements u_r = k_rr^-1 (f0_r - k_rk u_k).
  Vec12 condense_f0_with(const Mat12& k_raw, const Vec12& f0) const;
  Vec12 local_displacements_with(const Mat12& k_raw, const Vec12& u_nodal, const Vec12& f0) const;
  // Condense an arbitrary local matrix with this bar's release pattern.
  Mat12 condense_matrix(const Mat12& k) const;

  // 12x12 transform from nodal-global to element-end-local displacements
  // (rotation + eccentricity): u_local_ends = T * u_nodal.
  const Mat12& T() const { return T_; }

  const Mat12& k_local() const { return k_local_cond_; }  // condensed local stiffness

  // Condense span-load FEFs consistently with the releases.
  Vec12 condense_f0(const Vec12& f0) const;

 private:
  void build_local_stiffness();
  void build_transform();
  void build_condensation();

  const Model* model_;
  const Bar* bar_;
  const Material* mat_;
  const SectionProps* sec_;
  double L_ = 0;
  Mat3 R_ = Mat3::Identity();
  Mat12 T_ = Mat12::Identity();          // rotation * ecc
  Mat12 k_local_raw_ = Mat12::Zero();    // pure member, before condensation (no soil)
  Mat12 k_local_cond_ = Mat12::Zero();   // after release condensation (no soil)
  Mat12 k_soil_local_ = Mat12::Zero();   // Winkler bed (stiffness only, never in recovery)
  Mat12 ke_nodal_ = Mat12::Zero();       // T^t (k_cond + k_soil) T
  bool has_rel_ = false;
  std::vector<int> rel_, keep_;          // released / kept local DOF indices
  MatX krr_inv_krk_;                     // k_rr^-1 k_rk  (|rel| x 12)
  MatX krr_inv_;                         // k_rr^-1
};

// Shear parameters used across statics/dynamics.
struct ShearPhi {
  double phi_y = 0;  // bending about z (deflection along y)
  double phi_z = 0;  // bending about y (deflection along z)
};
ShearPhi shear_phi(const Material& m, const SectionProps& s, double L);

}  // namespace keta

namespace keta {
// ---- exact beam-column (stability-function) formulation ----
// ST1-style EXACT member stiffness under axial force N (tension-positive):
// the bending blocks use the Livesley/Chen-Lui stability functions —
// trigonometric for compression (u = L sqrt(-N/EI)), hyperbolic for tension —
// which solve the beam-column ODE exactly, so ONE element per member yields
// the exact second-order response and the exact critical load (the first
// tangent singularity of a pinned-pinned single element is exactly kL = pi).
// Below |u| < 1 the functions take their Taylor series in u^2 (exact at
// u = 0, where they reduce to the linear block s = 4, c = 1/2).
// Shear deformation is not folded into the stability functions (documented
// approximation for shear-flexible sections under axial force).
struct StabilityFunctions {
  double s = 4.0, c = 0.5;   // rotational stiffness + carry-over
  double u2 = 0.0;           // signed (kL)^2: >0 compression, <0 tension
};
StabilityFunctions stability_functions(double N, double EI, double L);
// A concentrated dislocation (Load::Kind::PrestrainPoint) given in GLOBAL
// axes (ST1 DEFOR CON … GLO) restated in the bar axes R (rows = local axes):
// the displacement jumps and the rotation jumps turn as vectors.
Load dislocation_local(const Load& l, const Mat3& R);

}  // namespace keta
