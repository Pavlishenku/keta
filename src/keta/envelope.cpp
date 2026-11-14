#include "keta/envelope.h"

namespace keta {

namespace {

void check_shapes(const std::vector<const CaseResult*>& sources) {
  require(!sources.empty(), "envelope: no sources");
  const CaseResult& a = *sources.front();
  for (const CaseResult* s : sources) {
    require(s->node_disp.size() == a.node_disp.size() &&
                s->reactions.size() == a.reactions.size() && s->bars.size() == a.bars.size() &&
                s->shells.size() == a.shells.size(),
            "envelope: source shapes differ");
    for (std::size_t b = 0; b < a.bars.size(); ++b)
      require(s->bars[b].stations.size() == a.bars[b].stations.size(),
              "envelope: station counts differ");
  }
}

// Generic per-scalar envelope: value_of(src) -> scalar; returns (lo, hi, src_lo, src_hi).
template <class F>
void env_scalar(std::size_t n, F&& value_of, double& lo, double& hi, unsigned& slo,
                unsigned& shi) {
  for (std::size_t i = 0; i < n; ++i) {
    double v = value_of(i);
    require(std::isfinite(v), "envelope: non-finite source value");
    if (i == 0 || v < lo) {
      lo = v;
      slo = static_cast<unsigned>(i);
    }
    if (i == 0 || v > hi) {
      hi = v;
      shi = static_cast<unsigned>(i);
    }
  }
}

Envelope envelope_impl(std::size_t n, const CaseResult& shape,
                       const std::function<double(std::size_t, const double*)>& get,
                       const std::vector<const CaseResult*>* plain,
                       const std::vector<FactoredSource>* factored) {
  (void)get;
  Envelope env;
  env.lo.value = shape;
  env.hi.value = shape;
  auto value = [&](std::size_t i, auto&& accessor) -> double {
    if (plain) return accessor(*(*plain)[i]);
    const FactoredSource& fs = (*factored)[i];
    double v = accessor(*fs.result);
    return v;  // factoring handled by caller wrappers below
  };
  (void)value;

  auto for_both = [&](auto&& accessor, auto&& setter, unsigned* slo_out, unsigned* shi_out) {
    double lo = 0, hi = 0;
    unsigned slo = 0, shi = 0;
    env_scalar(
        n,
        [&](std::size_t i) {
          if (plain) return accessor(*(*plain)[i]);
          const FactoredSource& fs = (*factored)[i];
          return accessor(*fs.result);
        },
        lo, hi, slo, shi);
    if (factored) {
      // corner rule: choose c per source term-by-term
      double lo2 = 0, hi2 = 0;
      unsigned l2 = 0, h2 = 0;
      env_scalar(
          n,
          [&](std::size_t i) {
            const FactoredSource& fs = (*factored)[i];
            double v = accessor(*fs.result);
            return std::min(fs.c_min * v, fs.c_max * v);
          },
          lo2, hi2, l2, h2);
      double lo3 = 0, hi3 = 0;
      unsigned l3 = 0, h3 = 0;
      env_scalar(
          n,
          [&](std::size_t i) {
            const FactoredSource& fs = (*factored)[i];
            double v = accessor(*fs.result);
            return std::max(fs.c_min * v, fs.c_max * v);
          },
          lo3, hi3, l3, h3);
      setter(env.lo.value, lo2);
      setter(env.hi.value, hi3);
      if (slo_out) *slo_out = l2;
      if (shi_out) *shi_out = h3;
      return;
    }
    setter(env.lo.value, lo);
    setter(env.hi.value, hi);
    if (slo_out) *slo_out = slo;
    if (shi_out) *shi_out = shi;
  };

  const CaseResult& a = shape;
  env.lo.disp_src.assign(a.node_disp.size(), {});
  env.hi.disp_src.assign(a.node_disp.size(), {});
  for (std::size_t ni = 0; ni < a.node_disp.size(); ++ni)
    for (int c = 0; c < 6; ++c)
      for_both([&](const CaseResult& r) { return r.node_disp[ni](c); },
               [&](CaseResult& r, double v) { r.node_disp[ni](c) = v; },
               &env.lo.disp_src[ni][static_cast<std::size_t>(c)],
               &env.hi.disp_src[ni][static_cast<std::size_t>(c)]);
  env.lo.reaction_src.assign(a.reactions.size(), {});
  env.hi.reaction_src.assign(a.reactions.size(), {});
  for (std::size_t ri = 0; ri < a.reactions.size(); ++ri)
    for (int c = 0; c < 6; ++c)
      for_both([&](const CaseResult& r) { return r.reactions[ri].R(c); },
               [&](CaseResult& r, double v) { r.reactions[ri].R(c) = v; },
               &env.lo.reaction_src[ri][static_cast<std::size_t>(c)],
               &env.hi.reaction_src[ri][static_cast<std::size_t>(c)]);
  env.lo.bar_src.resize(a.bars.size());
  env.hi.bar_src.resize(a.bars.size());
  for (std::size_t b = 0; b < a.bars.size(); ++b) {
    env.lo.bar_src[b].stations.resize(a.bars[b].stations.size());
    env.hi.bar_src[b].stations.resize(a.bars[b].stations.size());
    for (int c = 0; c < 12; ++c)
      for_both([&](const CaseResult& r) { return r.bars[b].end_forces_local(c); },
               [&](CaseResult& r, double v) { r.bars[b].end_forces_local(c) = v; },
               &env.lo.bar_src[b].end[static_cast<std::size_t>(c)],
               &env.hi.bar_src[b].end[static_cast<std::size_t>(c)]);
    for (std::size_t k = 0; k < a.bars[b].stations.size(); ++k) {
      auto comp = [&](const Station& s, int c) -> double {
        switch (c) {
          case 0: return s.N;
          case 1: return s.Vy;
          case 2: return s.Vz;
          case 3: return s.T;
          case 4: return s.My;
          default: return s.Mz;
        }
      };
      auto set_comp = [&](Station& s, int c, double v) {
        switch (c) {
          case 0: s.N = v; break;
          case 1: s.Vy = v; break;
          case 2: s.Vz = v; break;
          case 3: s.T = v; break;
          case 4: s.My = v; break;
          default: s.Mz = v; break;
        }
      };
      for (int c = 0; c < 6; ++c)
        for_both([&](const CaseResult& r) { return comp(r.bars[b].stations[k], c); },
                 [&](CaseResult& r, double v) { set_comp(r.bars[b].stations[k], c, v); },
                 &env.lo.bar_src[b].stations[k].forces[static_cast<std::size_t>(c)],
                 &env.hi.bar_src[b].stations[k].forces[static_cast<std::size_t>(c)]);
      for (int c = 0; c < 3; ++c)
        for_both([&](const CaseResult& r) { return r.bars[b].stations[k].defl(c); },
                 [&](CaseResult& r, double v) { r.bars[b].stations[k].defl(c) = v; },
                 &env.lo.bar_src[b].stations[k].defl[static_cast<std::size_t>(c)],
                 &env.hi.bar_src[b].stations[k].defl[static_cast<std::size_t>(c)]);
      for (int c = 0; c < 3; ++c)
        for_both([&](const CaseResult& r) { return r.bars[b].stations[k].rot(c); },
                 [&](CaseResult& r, double v) { r.bars[b].stations[k].rot(c) = v; },
                 &env.lo.bar_src[b].stations[k].rot[static_cast<std::size_t>(c)],
                 &env.hi.bar_src[b].stations[k].rot[static_cast<std::size_t>(c)]);
    }
  }

  // Facet shells: the 8 tensor components PLUS m1/m2/n1/n2 are enveloped as
  // INDEPENDENT scalars (a design check reads each bound on its own). The
  // per-case principal ANGLES have no meaning on a bound, so they are zeroed.
  env.lo.shell_src.assign(a.shells.size(), ShellSrc{});
  env.hi.shell_src.assign(a.shells.size(), ShellSrc{});
  for (std::size_t s = 0; s < a.shells.size(); ++s) {
    auto comp = [](const ShellResultants& v, int c) -> double {
      switch (c) {
        case 0: return v.nx;
        case 1: return v.ny;
        case 2: return v.nxy;
        case 3: return v.mx;
        case 4: return v.my;
        case 5: return v.mxy;
        case 6: return v.qx;
        case 7: return v.qy;
        case 8: return v.m1;
        case 9: return v.m2;
        case 10: return v.n1;
        default: return v.n2;
      }
    };
    auto set_comp = [](ShellResultants& v, int c, double x) {
      switch (c) {
        case 0: v.nx = x; break;
        case 1: v.ny = x; break;
        case 2: v.nxy = x; break;
        case 3: v.mx = x; break;
        case 4: v.my = x; break;
        case 5: v.mxy = x; break;
        case 6: v.qx = x; break;
        case 7: v.qy = x; break;
        case 8: v.m1 = x; break;
        case 9: v.m2 = x; break;
        case 10: v.n1 = x; break;
        default: v.n2 = x; break;
      }
    };
    for (int c = 0; c < 12; ++c)
      for_both([&](const CaseResult& r) { return comp(r.shells[s].r, c); },
               [&](CaseResult& r, double v) { set_comp(r.shells[s].r, c, v); },
               &env.lo.shell_src[s][static_cast<std::size_t>(c)],
               &env.hi.shell_src[s][static_cast<std::size_t>(c)]);
    env.lo.value.shells[s].r.m_theta = 0;
    env.lo.value.shells[s].r.n_theta = 0;
    env.hi.value.shells[s].r.m_theta = 0;
    env.hi.value.shells[s].r.n_theta = 0;
  }
  return env;
}

}  // namespace

Envelope envelope(const std::vector<const CaseResult*>& sources) {
  check_shapes(sources);
  return envelope_impl(sources.size(), *sources.front(), {}, &sources, nullptr);
}

Envelope envelope_factored(const std::vector<FactoredSource>& sources) {
  std::vector<const CaseResult*> plain;
  for (const auto& s : sources) plain.push_back(s.result);
  check_shapes(plain);
  return envelope_impl(sources.size(), *sources.front().result, {}, nullptr, &sources);
}

Concomitant envelope_concomitant(const std::vector<const CaseResult*>& sources,
                                 const Privileged& p) {
  check_shapes(sources);
  auto pivot = [&](const CaseResult& r) -> double {
    switch (p.what) {
      case Privileged::What::StationForce: {
        for (const auto& b : r.bars)
          if (b.bar == p.bar) {
            const Station& s = b.stations.at(p.station);
            switch (p.comp) {
              case 0: return s.N;
              case 1: return s.Vy;
              case 2: return s.Vz;
              case 3: return s.T;
              case 4: return s.My;
              default: return s.Mz;
            }
          }
        fail("concomitant: unknown bar");
      }
      case Privileged::What::Reaction:
        for (const auto& re : r.reactions)
          if (re.node == p.node) return re.R(p.comp);
        fail("concomitant: unknown reaction node");
      case Privileged::What::Displacement:
        return r.node_disp.at(static_cast<std::size_t>(p.node))(p.comp);
    }
    return 0;
  };
  Concomitant c;
  for (std::size_t i = 0; i < sources.size(); ++i) {
    double v = pivot(*sources[i]);
    require(std::isfinite(v), "concomitant: non-finite pivot");
    if (i == 0 || v > c.max) {
      c.max = v;
      c.src_max = i;
      c.state_max = sources[i];
    }
    if (i == 0 || v < c.min) {
      c.min = v;
      c.src_min = i;
      c.state_min = sources[i];
    }
  }
  return c;
}

}  // namespace keta
