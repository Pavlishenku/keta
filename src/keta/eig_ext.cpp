#include "keta/eig_ext.h"

namespace keta {

// ---------------- deflated inverse iteration (ST1 METHODE 1) ----------------

std::vector<PencilPair> inverse_iteration(const FactorCache& Kfac, const SpMat& K, const SpMat& B,
                                          int n_wanted, int max_iters, double tol) {
  require(n_wanted > 0, "inverse iteration: n_wanted must be positive");
  require(max_iters > 0, "inverse iteration: max_iters must be positive");
  require(tol > 0, "inverse iteration: tolerance must be positive");
  Eigen::Index n = K.rows();
  int want = std::min<int>(n_wanted, static_cast<int>(n));

  std::vector<VecX> V, KV;  // converged vectors (K-orthonormal) and their K-images
  std::vector<PencilPair> out;

  auto k_norm = [&](const VecX& v) { return std::sqrt(std::max(0.0, v.dot(VecX(K * v)))); };
  // Two classical Gram-Schmidt passes in the K inner product (same policy as
  // the Lanczos kernel: one pass loses orthogonality on ill-conditioned K).
  auto deflate = [&](VecX& v) {
    for (int pass = 0; pass < 2; ++pass)
      for (std::size_t i = 0; i < V.size(); ++i) v -= v.dot(KV[i]) * V[i];
  };

  for (int mode = 0; mode < want; ++mode) {
    // Deterministic start vector, shifted per mode so a fully deflated start
    // cannot repeat: 1, 1/2, 1/3, ... offset by the mode index.
    VecX x(n);
    for (Eigen::Index i = 0; i < n; ++i) x(i) = 1.0 / static_cast<double>(i + 1 + mode);
    deflate(x);
    double nx = k_norm(x);
    if (!(nx > 0)) break;  // start vector entirely inside the converged span
    x /= nx;

    double mu = 0.0, mu_prev = 0.0;
    bool have_prev = false, converged = false, stationary = false;
    for (int it = 0; it < max_iters; ++it) {
      VecX y = Kfac.solve(VecX(B * x));  // one step of the power iteration on K^-1 B
      deflate(y);
      double ny = k_norm(y);
      if (!(ny > 0)) break;  // B x has no component left outside the converged span
      y /= ny;
      mu = y.dot(VecX(B * y));  // Rayleigh quotient (y is K-normalized)
      x = std::move(y);
      // The Rayleigh quotient settles quadratically, the vector only
      // linearly: next to a close mode (e.g. periods 0.0233 / 0.0222 s) mu is
      // stationary while the vector still fails the residual certificate, so
      // keep iterating until the vector is certified too.
      stationary = have_prev && std::fabs(mu - mu_prev) <= tol * std::max(std::fabs(mu), 1e-300);
      if (stationary && pencil_residual(K, B, mu, x) <= 1e-7) {
        converged = true;
        break;
      }
      mu_prev = mu;
      have_prev = true;
    }
    // A stationary quotient whose vector never reached 1e-7 still goes to the
    // certificate below (1e-6), as before.
    if (!converged && stationary) converged = true;
    if (!converged) break;  // stalled (e.g. two |mu| of equal size, opposite signs)

    PencilPair p;
    p.mu = mu;
    p.x = x;
    p.residual = pencil_residual(K, B, mu, x);
    if (p.residual > 1e-6) break;  // never return an uncertified pair
    out.push_back(std::move(p));
    V.push_back(x);
    KV.push_back(VecX(K * x));
  }
  return out;
}

// ---------------- CTE/VAR buckling ----------------

namespace {

// Freeze the axial + shell membrane fields of one load case (id 0 = "no axial
// state" -> zeros). The shell field rides along so a shell-bearing model's
// CTE/VAR states destabilise the plates too; the per-bar `_fields` entry
// points cannot carry one and leave it zero (shells stiffness + mass only,
// documented in eig_ext.h).
struct FrozenState {
  std::vector<double> N;
  std::vector<Vec3> Sn;
};
FrozenState resolve_state(const Model& model, StaticContext& ctx, Id case_id) {
  FrozenState st;
  st.N.assign(model.bars.size(), 0.0);
  st.Sn.assign(model.shells.size(), Vec3::Zero());
  if (case_id == 0) return st;
  const LoadCase& lc = model.load_case(case_id);
  StaticOptions sopt;
  sopt.stations_per_bar = 2;
  CaseResult r = ctx.solve_case(lc, sopt);
  st.N = axial_field(r);
  st.Sn = shell_membrane_field(r);
  return st;
}

// st1.par tension policy: by default the tension of the applied cases does not
// enter Kg (see the header quote) — clamp the field to compression only.
void apply_tension_policy(std::vector<double>& N, bool tension_stiffening) {
  if (tension_stiffening) return;
  for (double& v : N)
    if (v > 0) v = 0.0;
}

bool all_zero(const std::vector<double>& N) {
  for (double v : N)
    if (v != 0.0) return false;
  return true;
}

bool shells_all_zero(const std::vector<Vec3>& Sn) {
  for (const Vec3& v : Sn)
    if (v.squaredNorm() != 0.0) return false;
  return true;
}

// ST1: « Si on met une charge 0 sur la composante variable, ST1 considère un
// effort normal uniforme de compression des barres » — opt-in here, and then
// lambda is a force in N (the uniform compression is -1 N per bar). A VAR
// case whose destabilising action lives in the SHELLS is a legitimate
// non-zero variable state, so the fallback/error only fires when both fields
// are zero.
void handle_zero_var(std::vector<double>& Nv, const std::vector<Vec3>& Snv,
                     bool uniform_if_zero) {
  if (!all_zero(Nv) || !shells_all_zero(Snv)) return;
  require(uniform_if_zero,
          "buckling cte/var: the VAR axial field is zero — nothing to amplify (a VAR case in pure "
          "tension is zeroed too by the default st1.par tension policy: set tension_stiffening to "
          "keep it, or uniform_var_if_zero for ST1's uniform unit compression fallback)");
  for (double& v : Nv) v = -1.0;
}

void check_field_size(const Model& model, std::vector<double>& N, const char* what) {
  if (N.empty()) {
    N.assign(model.bars.size(), 0.0);
    return;
  }
  require(N.size() == model.bars.size(),
          std::string("buckling cte/var: ") + what + " axial field size mismatch");
  for (double v : N)
    require(finite(v), std::string("buckling cte/var: non-finite ") + what + " axial force");
}

// K' = K + Kg(N_cte), factorized, with the "already beyond critical" gate.
FactorCache factor_cte_state(const Model& model, const StaticContext& ctx,
                             const std::vector<double>& Nc, const std::vector<Vec3>& Snc,
                             SpMat& Kp, Inertia& in) {
  SpMat Kg_c = assemble_geometric(model, ctx.asm_, Nc, &Snc);
  Kp = (ctx.K + Kg_c).eval();
  FactorCache f;
  try {
    f.factorize(Kp);
  } catch (const Error& e) {
    fail(std::string("buckling cte/var: K + Kg(N_cte) is singular — the CTE (fixed) axial state "
                     "sits exactly at a critical state of the structure: ") +
         e.what());
  }
  in = f.inertia();
  require(in.n_neg == 0,
          "buckling cte/var: the CTE (fixed) axial state is already beyond critical — K + Kg(N_cte) "
          "has " + std::to_string(in.n_neg) +
              " negative pivot(s); reduce the fixed case (or check its sign) before amplifying the "
              "variable one");
  return f;
}

BucklingCteVar cte_var_core(const Model& model, StaticContext& ctx, std::vector<double> Nc,
                            std::vector<double> Nv, std::vector<Vec3> Snc, std::vector<Vec3> Snv,
                            int n_modes, const BucklingCteVarOptions& opt) {
  require(n_modes > 0, "buckling cte/var: n_modes must be positive");
  apply_tension_policy(Nc, opt.tension_stiffening);
  apply_tension_policy(Nv, opt.tension_stiffening);
  // The st1.par tension policy is a BAR-AXIAL concept; a shell's membrane
  // TENSOR enters Kg whole (its stabilising and destabilising directions are
  // not separable by a scalar sign).
  if (Snc.empty()) Snc.assign(model.shells.size(), Vec3::Zero());
  if (Snv.empty()) Snv.assign(model.shells.size(), Vec3::Zero());
  handle_zero_var(Nv, Snv, opt.uniform_var_if_zero);

  SpMat Kp;
  Inertia in;
  FactorCache fp = factor_cte_state(model, ctx, Nc, Snc, Kp, in);
  SpMat Kg_v = assemble_geometric(model, ctx.asm_, Nv, &Snv);

  std::vector<PencilPair> pairs;
  if (opt.method == EigMethod::InverseIteration) {
    // Extremal |mu| first = smallest |lambda| first, and the smallest |lambda|
    // may sit on the negative (reversed-loading) side: collect a margin so the
    // positive factors asked for are actually among them.
    int n_collect = std::min(ctx.asm_.dofs.n_free, 2 * n_modes + 4);
    pairs = inverse_iteration(fp, Kp, Kg_v, n_collect, opt.max_iters);
  } else if (ctx.asm_.dofs.n_free <= opt.dense_threshold) {
    pairs = dense_pencil(Kp, Kg_v);
  } else {
    pairs = k_inner_lanczos(fp, Kp, Kg_v, n_modes);
  }

  // K' + lambda Kg_v singular  <=>  Kg_v x = -(1/lambda) K' x: the smallest
  // positive lambda is the most negative mu. Ascending lambda.
  BucklingCteVar out;
  out.N_cte = std::move(Nc);
  out.N_var = std::move(Nv);
  out.cte_inertia = in;
  std::vector<std::pair<double, const PencilPair*>> lam;
  for (const PencilPair& p : pairs)
    if (p.mu < -1e-14 && p.residual < 1e-6) lam.emplace_back(-1.0 / p.mu, &p);
  std::sort(lam.begin(), lam.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
  for (const auto& [l, p] : lam) {
    if (static_cast<int>(out.load_factors.size()) >= n_modes) break;
    out.load_factors.push_back(l);
    out.modes.push_back(expand_unit_crest(ctx.asm_.dofs, p->x));
    out.residuals.push_back(p->residual);
    double best = 0;  // the unit-crest scaling of expand_unit_crest
    for (Eigen::Index i = 0; i < p->x.size(); ++i)
      if (std::fabs(p->x(i)) > std::fabs(best)) best = p->x(i);
    const VecX xs = best != 0 ? VecX(p->x / best) : p->x;
    const VecX kgx = -(Kg_v * xs);
    const double g = xs.dot(kgx);
    out.gen_stiffness.push_back(g);
    const double gamma = g > 0 ? kgx.sum() / std::sqrt(g) : 0.0;  // phi^T KG R, R = 1
    out.participation.push_back(gamma);
  }
  {
    const VecX ones = VecX::Ones(Kg_v.rows());
    out.total_kg = -ones.dot(Kg_v * ones);
  }
  return out;
}

double cte_var_bisection_core(const Model& model, StaticContext& ctx, std::vector<double> Nc,
                              std::vector<double> Nv, std::vector<Vec3> Snc,
                              std::vector<Vec3> Snv, const BucklingCteVarBisectOptions& opt) {
  require(opt.lambda_max > 0 && std::isfinite(opt.lambda_max),
          "buckling cte/var bisection: lambda_max must be a positive finite cap");
  require(opt.tol_rel > 0, "buckling cte/var bisection: tolerance must be positive");
  apply_tension_policy(Nc, opt.tension_stiffening);
  apply_tension_policy(Nv, opt.tension_stiffening);
  if (Snc.empty()) Snc.assign(model.shells.size(), Vec3::Zero());
  if (Snv.empty()) Snv.assign(model.shells.size(), Vec3::Zero());
  handle_zero_var(Nv, Snv, opt.uniform_var_if_zero);
  SpMat Kg_c = assemble_geometric(model, ctx.asm_, Nc, &Snc);
  SpMat Kg_v = assemble_geometric(model, ctx.asm_, Nv, &Snv);

  auto flipped = [&](double lambda) {
    SpMat Kt;
    if (opt.exact_stability) {
      // Exact stability-function tangent at the COMBINED axial state
      // N = N_cte + lambda N_var: the first inertia flip is then the exact
      // critical multiplier even with one element per member.
      std::vector<Triplet> trips;
      for (std::size_t e = 0; e < ctx.asm_.elements.size(); ++e) {
        const Bar& b = model.bars[e];
        if (b.kind == BarKind::Cable) continue;
        Mat12 ke;
        try {
          ke = ctx.asm_.elements[e].ke_exact(Nc[e] + lambda * Nv[e]);
        } catch (const Error&) {
          return true;  // a member driven past its own trig range is beyond critical
        }
        int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
        int g[12];
        for (int c = 0; c < 6; ++c) {
          g[c] = ctx.asm_.dofs.free_ix[static_cast<std::size_t>(ctx.asm_.dofs.full(ni, c))];
          g[6 + c] = ctx.asm_.dofs.free_ix[static_cast<std::size_t>(ctx.asm_.dofs.full(nj, c))];
        }
        for (int a = 0; a < 12; ++a)
          for (int bb = 0; bb < 12; ++bb)
            if (g[a] >= 0 && g[bb] >= 0 && ke(a, bb) != 0)
              trips.emplace_back(g[a], g[bb], ke(a, bb));
      }
      for (std::size_t f = 0; f < ctx.asm_.spring.size(); ++f)
        if (ctx.asm_.spring[f] != 0) {
          int fi = ctx.asm_.dofs.free_ix[f];
          if (fi >= 0) trips.emplace_back(fi, fi, ctx.asm_.spring[f]);
        }
      // Facet shells stay LINEAR under the exact bar tangents: ke plus the
      // membrane Kg at the combined state n_cte + lambda n_var (no
      // stability-function analogue exists for a plate).
      for (const ShellElem& se : ctx.asm_.shells) {
        Vec3 s_n = Snc[se.index] + lambda * Snv[se.index];
        Mat18 kes = se.el.ke();
        if (s_n.squaredNorm() != 0) kes += se.el.kg(s_n(0), s_n(1), s_n(2));
        int gd[18];
        for (int k = 0; k < 3; ++k)
          for (int c = 0; c < 6; ++c)
            gd[6 * k + c] = ctx.asm_.dofs.free_ix[static_cast<std::size_t>(
                ctx.asm_.dofs.full(se.nodes[static_cast<std::size_t>(k)], c))];
        for (int a = 0; a < 18; ++a)
          for (int bb = 0; bb < 18; ++bb)
            if (gd[a] >= 0 && gd[bb] >= 0 && kes(a, bb) != 0)
              trips.emplace_back(gd[a], gd[bb], kes(a, bb));
      }
      Kt = SpMat(ctx.asm_.dofs.n_free, ctx.asm_.dofs.n_free);
      Kt.setFromTriplets(trips.begin(), trips.end());
    } else {
      Kt = ctx.K + Kg_c + lambda * Kg_v;
    }
    FactorCache f;
    try {
      f.factorize(Kt);
    } catch (const Error&) {
      return true;  // exactly singular: counted as flipped
    }
    return f.inertia().n_neg > 0;
  };

  require(!flipped(0.0),
          "buckling cte/var bisection: the CTE (fixed) axial state is already beyond critical");
  double lo = 0.0, hi = 1.0;
  while (!flipped(hi)) {
    lo = hi;
    hi *= 2;
    require(hi <= opt.lambda_max, "buckling cte/var bisection: no critical factor below cap");
  }
  while (hi - lo > opt.tol_rel * hi) {
    double mid = 0.5 * (lo + hi);
    (flipped(mid) ? hi : lo) = mid;
  }
  return 0.5 * (lo + hi);
}

}  // namespace

BucklingCteVar solve_buckling_cte_var(const Model& model, Id case_cte_id, Id case_var_id,
                                      int n_modes, const BucklingCteVarOptions& opt) {
  require(model.continua.empty(), "solve_buckling_cte_var: continuum elements are not supported (linear statics, modal and dynamics take them)");
  model.validate();
  StaticContext ctx = StaticContext::build(model);
  ctx.cache.require_spd();
  FrozenState c = resolve_state(model, ctx, case_cte_id);
  FrozenState v = resolve_state(model, ctx, case_var_id);
  return cte_var_core(model, ctx, std::move(c.N), std::move(v.N), std::move(c.Sn),
                      std::move(v.Sn), n_modes, opt);
}

BucklingCteVar solve_buckling_cte_var_fields(const Model& model, const std::vector<double>& N_cte,
                                             const std::vector<double>& N_var, int n_modes,
                                             const BucklingCteVarOptions& opt) {
  require(model.continua.empty(), "solve_buckling_cte_var_fields: continuum elements are not supported (linear statics, modal and dynamics take them)");
  model.validate();
  std::vector<double> Nc = N_cte, Nv = N_var;
  check_field_size(model, Nc, "CTE");
  check_field_size(model, Nv, "VAR");
  StaticContext ctx = StaticContext::build(model);
  ctx.cache.require_spd();
  // Per-bar signature: no shell prestress state can be carried, so the shells
  // enter with stiffness + mass only (empty fields -> zeros in the core).
  return cte_var_core(model, ctx, std::move(Nc), std::move(Nv), {}, {}, n_modes, opt);
}

double buckling_cte_var_bisection(const Model& model, Id case_cte_id, Id case_var_id,
                                  const BucklingCteVarBisectOptions& opt) {
  require(model.continua.empty(), "buckling_cte_var_bisection: continuum elements are not supported (linear statics, modal and dynamics take them)");
  model.validate();
  StaticContext ctx = StaticContext::build(model);
  FrozenState c = resolve_state(model, ctx, case_cte_id);
  FrozenState v = resolve_state(model, ctx, case_var_id);
  return cte_var_bisection_core(model, ctx, std::move(c.N), std::move(v.N), std::move(c.Sn),
                                std::move(v.Sn), opt);
}

double buckling_cte_var_bisection_fields(const Model& model, const std::vector<double>& N_cte,
                                         const std::vector<double>& N_var,
                                         const BucklingCteVarBisectOptions& opt) {
  require(model.continua.empty(), "buckling_cte_var_bisection_fields: continuum elements are not supported (linear statics, modal and dynamics take them)");
  model.validate();
  std::vector<double> Nc = N_cte, Nv = N_var;
  check_field_size(model, Nc, "CTE");
  check_field_size(model, Nv, "VAR");
  StaticContext ctx = StaticContext::build(model);
  return cte_var_bisection_core(model, ctx, std::move(Nc), std::move(Nv), {}, {}, opt);
}

// ---------------- modal drivers ----------------

ModalResults solve_modal_method(const Model& model, const ModalExtOptions& opt) {
  if (opt.method == EigMethod::Subspace) return solve_modal(model, opt);

  model.validate();
  Assembled asm_ = assemble_elements(model, false, opt.axial_state);
  require(asm_.dofs.n_free > 0, "modal: no free DOFs");
  SpMat K = assemble_stiffness(model, asm_);
  if (opt.axial_state) {
    SpMat Kg = assemble_geometric(model, asm_, *opt.axial_state);
    K = (K + Kg).eval();  // vibration about the loaded state
  }
  SpMat M = assemble_mass(model, asm_, opt.lumped_mass);
  require(M.nonZeros() > 0, "modal: model has no mass (set rho or add point masses)");

  FactorCache f;
  f.factorize(K);
  f.require_spd();
  std::vector<PencilPair> pairs = inverse_iteration(f, K, M, opt.n_modes, opt.max_iters);
  return modal_from_pairs(model, asm_, M, pairs, opt.n_modes);
}

namespace {

// Drop every mode past `keep` from a ModalResults (parallel arrays).
void truncate_modes(ModalResults& r, std::size_t keep) {
  if (keep >= r.omegas.size()) return;
  r.omegas.resize(keep);
  r.frequencies_hz.resize(keep);
  r.modes_free.resize(keep);
  r.modes.resize(keep);
  r.residuals.resize(keep);
  r.participation.resize(keep);
  r.effective_mass.resize(keep);
}

}  // namespace

ModalFrequencyResults solve_modal_frequency(const Model& model, double f_cutoff_hz,
                                            const ModalExtOptions& opt) {
  require(std::isfinite(f_cutoff_hz) && f_cutoff_hz > 0,
          "modal FREQUENCE: the cutoff must be a positive finite frequency [Hz]");
  require(opt.axial_state == nullptr,
          "modal FREQUENCE: a prestress axial state is not supported here — the Sturm count sizes "
          "the problem on (K, M) alone");
  ModalFrequencyResults out;
  out.cutoff_hz = f_cutoff_hz;
  out.expected_count = mode_count_in_band(model, 0.0, f_cutoff_hz, opt.lumped_mass);
  if (out.expected_count <= 0) return out;  // nothing below the cutoff

  ModalExtOptions o = opt;
  o.n_modes = out.expected_count;
  out.modal = solve_modal_method(model, o);

  // Frequencies are ascending: keep the ones at or below the cutoff.
  std::size_t keep = 0;
  while (keep < out.modal.frequencies_hz.size() &&
         out.modal.frequencies_hz[keep] <= f_cutoff_hz * (1 + 1e-9))
    ++keep;
  truncate_modes(out.modal, keep);
  out.complete = static_cast<int>(keep) == out.expected_count;
  if (!out.complete)
    out.warning = "modal FREQUENCE: the Sturm count reports " + std::to_string(out.expected_count) +
                  " mode(s) below " + std::to_string(f_cutoff_hz) + " Hz but only " +
                  std::to_string(keep) +
                  " converged — raise the iteration count or use METHODE 2 (subspace)";
  return out;
}

ModalRatioResults solve_modal_ratio(const Model& model, double target_pct, int nm_increment,
                                    const ModalExtOptions& opt) {
  require(std::isfinite(target_pct) && target_pct > 0 && target_pct <= 100,
          "modal RATIO: the participation target must lie in (0, 100] %");
  require(nm_increment > 0, "modal RATIO: the AVANCE increment must be positive");
  Assembled asm_ = assemble_elements(model);
  int n_free = asm_.dofs.n_free;
  require(n_free > 0, "modal RATIO: no free DOFs");

  ModalRatioResults out;
  int n = std::min(nm_increment, n_free);
  while (true) {
    ModalExtOptions o = opt;
    o.n_modes = n;
    ModalResults r = solve_modal_method(model, o);
    int got = static_cast<int>(r.omegas.size());

    // Cumulated ratio alpha_j = M_j / M_tot per direction, mode by mode (p. 195).
    std::vector<std::array<double, 6>> cum;
    std::array<double, 6> run{};
    for (int i = 0; i < got; ++i) {
      for (int d = 0; d < 6; ++d) {
        double tot = r.total_mass[static_cast<std::size_t>(d)];
        if (tot > 0)
          run[static_cast<std::size_t>(d)] +=
              100.0 * r.effective_mass[static_cast<std::size_t>(i)][static_cast<std::size_t>(d)] /
              tot;
      }
      cum.push_back(run);
    }
    std::array<bool, 6> has_mass{};
    for (int d = 0; d < 6; ++d)
      has_mass[static_cast<std::size_t>(d)] = r.total_mass[static_cast<std::size_t>(d)] > 0;
    bool reached = true;
    for (int d = 0; d < 3; ++d)  // the target applies to the translational directions
      if (has_mass[static_cast<std::size_t>(d)] && run[static_cast<std::size_t>(d)] < target_pct)
        reached = false;

    out.modal = std::move(r);
    out.cumulative_pct = std::move(cum);
    out.reached_pct = run;
    out.has_mass = has_mass;
    out.n_modes = got;
    out.target_reached = reached;
    if (reached) {
      out.warning.clear();
      break;
    }
    if (got < n || n >= n_free) {  // the spectrum (or the solver) is exhausted
      out.warning = "modal RATIO: " + std::to_string(target_pct) +
                    " % not reached with " + std::to_string(got) +
                    " mode(s) — the computed spectrum is exhausted (free DOFs: " +
                    std::to_string(n_free) + ")";
      break;
    }
    n = std::min(n + nm_increment, n_free);
  }
  return out;
}

}  // namespace keta
