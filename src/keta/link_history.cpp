#include "keta/link_history.h"

#include <algorithm>
#include <cmath>

namespace keta {

namespace {
void check_law(const PlasticLink& l) {
  require(std::isfinite(l.k) && std::isfinite(l.fy) && std::isfinite(l.kp) && l.k > 0 && l.fy > 0 && l.kp >= 0 &&
              l.kp < l.k,
          "link history: a plastic law needs k > 0, fy > 0, 0 <= kp < k");
}
double hardening(const PlasticLink& l) { return l.kp > 0 ? l.k * l.kp / (l.k - l.kp) : 0.0; }
}  // namespace

LinkHistory::LinkHistory(const Model& model, const std::vector<Component>& plastic)
    : model_(std::make_shared<const Model>(model)) {
  std::vector<std::pair<std::size_t, int>> extra;
  for (const Component& c : plastic) {
    check_law(c.law);
    std::size_t li = model.links.size();
    for (std::size_t i = 0; i < model.links.size(); ++i)
      if (model.links[i].id == c.link) li = i;
    require(li < model.links.size(), "link history: no link " + std::to_string(c.link));
    require(c.comp >= 0 && c.comp < 6, "link history: component must be 0..5");
    require(!model.links[li].law[static_cast<std::size_t>(c.comp)],
            "link history: link " + std::to_string(c.link) + " component already carries an elastic law");
    extra.emplace_back(li, c.comp);
  }
  core_ = std::make_shared<const LinkLawCore>(*model_, extra);
  n_elastic_ = core_->comps().size() - plastic.size();
  for (const Component& c : plastic) st_.push_back({c.law, 0.0, 0.0});
  rhs_ = VecX::Zero(core_->n_free());
  v_ = VecX::Zero(core_->n_free());
  d0_.assign(core_->comps().size(), 0.0);
}

std::size_t LinkHistory::index(Id link, int comp) const {
  for (std::size_t j = n_elastic_; j < core_->comps().size(); ++j) {
    const LinkLawCore::Comp& q = core_->comps()[j];
    if (model_->links[q.link].id == link && q.c == comp) return j;
  }
  throw Error("link history: link " + std::to_string(link) + " component " + std::to_string(comp) +
              " is not a plastic component");
}

LinkLaw LinkHistory::step_law(const State& s) const {
  const PlasticLink& l = s.law;
  const double dm = s.up + (s.beta - l.fy) / l.k, dp = s.up + (s.beta + l.fy) / l.k;
  const double ext = std::max(1.0, dp - dm);
  return LinkLaw{{{dm - ext, s.beta - l.fy - l.kp * ext},
                  {dm, s.beta - l.fy},
                  {dp, s.beta + l.fy},
                  {dp + ext, s.beta + l.fy + l.kp * ext}}};
}

void LinkHistory::step(const VecX& rhs, const std::vector<double>& d0) {
  std::vector<LinkLaw> plastic_laws;
  plastic_laws.reserve(st_.size());
  for (const State& s : st_) plastic_laws.push_back(step_law(s));
  std::vector<const LinkLaw*> laws;
  for (std::size_t j = 0; j < n_elastic_; ++j) {
    const LinkLawCore::Comp& q = core_->comps()[j];
    laws.push_back(&*model_->links[q.link].law[static_cast<std::size_t>(q.c)]);
  }
  for (const LinkLaw& l : plastic_laws) laws.push_back(&l);
  core_->solve(rhs, d0, laws, v_);
  // commit: the closed-form return mapping at the converged deformation
  for (std::size_t i = 0; i < st_.size(); ++i) {
    State& s = st_[i];
    const double d = core_->deformation(n_elastic_ + i, v_, d0);
    const double ft = s.law.k * (d - s.up), H = hardening(s.law);
    const double over = std::fabs(ft - s.beta) - s.law.fy;
    if (over > 0) {
      const double sg = ft > s.beta ? 1.0 : -1.0, dg = over / (s.law.k + H);
      s.up += sg * dg;
      s.beta += sg * H * dg;
    }
  }
}

void LinkHistory::apply(Id case_id, int steps) {
  require(steps >= 1, "link history: steps >= 1");
  const LoadCase& lc = model_->load_case(case_id);
  const CaseLoads cl = core_->lower(lc);
  const std::vector<double> dd = core_->d0(cl);
  for (int s = 1; s <= steps; ++s) {
    const double lam = static_cast<double>(s) / steps;
    std::vector<double> d0 = d0_;
    for (std::size_t j = 0; j < d0.size(); ++j) d0[j] += lam * dd[j];
    step(rhs_ + lam * cl.rhs, d0);
  }
  rhs_ += cl.rhs;
  for (std::size_t j = 0; j < d0_.size(); ++j) d0_[j] += dd[j];
  applied_.push_back(case_id);
}

void LinkHistory::set_law(Id link, int comp, const PlasticLink& law) {
  check_law(law);
  const std::size_t j = index(link, comp);
  State& s = st_[j - n_elastic_];
  const double d = core_->deformation(j, v_, d0_);
  const double F = std::clamp(s.law.k * (d - s.up), -law.fy, law.fy);
  s.law = law;
  s.beta = 0;
  s.up = d - F / law.k;
}

double LinkHistory::deformation(Id link, int comp) const { return core_->deformation(index(link, comp), v_, d0_); }
double LinkHistory::slip(Id link, int comp) const { return st_[index(link, comp) - n_elastic_].up; }
double LinkHistory::force(Id link, int comp) const {
  const std::size_t j = index(link, comp);
  const State& s = st_[j - n_elastic_];
  return s.law.k * (core_->deformation(j, v_, d0_) - s.up);
}

CaseResult LinkHistory::result(int stations) const {
  // the cumulative loads: every applied case's loads in one case
  LoadCase all;
  all.id = applied_.empty() ? 0 : applied_.back();
  all.name = "history";
  // (prescribed displacements of one node summed into one load: the stages
  // add their increments)
  for (Id c : applied_)
    for (const Load& l : model_->load_case(c).loads) {
      if (l.kind == Load::Kind::SupportDisp) {
        auto it = std::find_if(all.loads.begin(), all.loads.end(), [&](const Load& q) {
          return q.kind == Load::Kind::SupportDisp && q.node == l.node;
        });
        if (it != all.loads.end()) {
          for (std::size_t k = 0; k < 6; ++k)
            if (l.disp[k]) it->disp[k] = it->disp[k].value_or(0.0) + *l.disp[k];
          continue;
        }
      }
      all.loads.push_back(l);
    }
  const CaseLoads cl = core_->lower(all);
  std::vector<double> F(core_->comps().size());
  for (std::size_t j = 0; j < F.size(); ++j) {
    const double d = core_->deformation(j, v_, d0_);
    if (j < n_elastic_) {
      const LinkLawCore::Comp& q = core_->comps()[j];
      F[j] = model_->links[q.link].law[static_cast<std::size_t>(q.c)]->force(d);
    } else {
      const State& s = st_[j - n_elastic_];
      F[j] = s.law.k * (d - s.up);
    }
  }
  return core_->recover(cl, v_, F, stations, all.id, all.name);
}

}  // namespace keta
