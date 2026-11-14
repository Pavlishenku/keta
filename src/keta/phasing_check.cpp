#include "keta/phasing_check.h"

#include <algorithm>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>

namespace keta {

namespace {
struct TState {
  double applied = 0;  // jack fraction applied so far
  bool live = false, on_bed = false, unbonded = false;
};
std::string id(Id x) { return std::to_string(x); }
}  // namespace

std::vector<PlanIssue> check_plan(const Model& model, const PhasePlan& plan) {
  std::vector<PlanIssue> out;
  std::size_t ip = 0;
  auto err = [&](const char* kind, const std::string& msg) {
    out.push_back({PlanIssue::Severity::Erreur, ip, kind, msg});
  };
  auto warn = [&](const char* kind, const std::string& msg) {
    out.push_back({PlanIssue::Severity::Avertissement, ip, kind, msg});
  };
  // model indexes
  std::unordered_map<Id, const Bar*> bars;
  for (const Bar& b : model.bars) bars.emplace(b.id, &b);
  std::unordered_set<Id> nodes, cases, elinks, rlinks, shell_t;
  for (const Node& n : model.nodes) nodes.insert(n.id);
  for (const LoadCase& c : model.cases) cases.insert(c.id);
  for (const ElasticLink& l : model.links) elinks.insert(l.id);
  for (const RigidLink& l : model.rigid_links) rlinks.insert(l.id);
  for (const ShellTendon& t : model.shell_tendons) shell_t.insert(t.id);
  std::unordered_map<Id, const Tendon*> tendons;
  for (const Tendon& t : model.tendons) tendons.emplace(t.id, &t);
  std::unordered_map<Id, std::vector<Id>> tendons_on_bar;  // bar -> tendons through it
  for (const Tendon& t : model.tendons)
    for (Id b : t.bars) tendons_on_bar[b].push_back(t.id);

  // the initial state, as solve_phased builds it
  std::set<Id> later_bars, later_el, later_rl;
  std::set<std::pair<Id, int>> later_parts;
  for (const Phase& ph : plan.phases) {
    later_bars.insert(ph.activate_bars.begin(), ph.activate_bars.end());
    later_el.insert(ph.activate_links.begin(), ph.activate_links.end());
    later_rl.insert(ph.activate_rigid_links.begin(), ph.activate_rigid_links.end());
    later_parts.insert(ph.activate_parts.begin(), ph.activate_parts.end());
  }
  std::unordered_set<Id> active;
  for (const Bar& b : model.bars)
    if (!later_bars.count(b.id) && b.kind != BarKind::Cable) active.insert(b.id);  // a catenary enters by name only
  std::unordered_set<Id> catenary;  // installed catenary cables (Phase::cable_actions)
  std::unordered_set<Id> el_on, rl_on;
  for (Id l : elinks)
    if (!later_el.count(l)) el_on.insert(l);
  for (Id l : rlinks)
    if (!later_rl.count(l)) rl_on.insert(l);
  std::map<Id, const Support*> supports;  // active, non-implicit, by node
  for (const Support& s : model.supports)
    if (!s.implicit) supports[s.node] = &s;
  std::set<std::pair<Id, int>> cast;  // composite parts cast so far
  for (const auto& [bid, cb] : model.composite)
    for (int i = 0; i < static_cast<int>(cb.parts.size()); ++i)
      if (!later_parts.count({bid, i})) cast.insert({bid, i});
  std::unordered_map<Id, int> degree;  // node -> active bars + supports there
  for (Id b : active) {
    ++degree[bars.at(b)->ni];
    ++degree[bars.at(b)->nj];
  }
  for (const auto& [n, s] : supports) ++degree[n];
  std::unordered_map<Id, TState> ts;
  std::unordered_map<Id, Releases> rel;  // bar releases changed by the plan so far
  std::unordered_set<Id> sections;
  for (const Section& x : model.sections) sections.insert(x.id);
  double date = plan.phases.empty() ? 0.0 : plan.phases.front().time;

  if (plan.phases.empty()) {
    err("plan_vide", "le plan n'a aucune phase");
    return out;
  }
  for (ip = 0; ip < plan.phases.size(); ++ip) {
    const Phase& ph = plan.phases[ip];
    if (ph.time < date) err("date_retrograde", "date " + std::to_string(ph.time) + " antérieure à " + std::to_string(date));
    date = std::max(date, ph.time);
    // 1. activations: supports, bars, links, parts
    for (const Support& s : ph.add_supports) {
      if (!nodes.count(s.node)) {
        err("noeud_inconnu", "appui ajouté au nœud inconnu " + id(s.node));
        continue;
      }
      if (!s.implicit && supports.emplace(s.node, &s).second) ++degree[s.node];
    }
    std::unordered_set<Id> cast_now;
    std::vector<const Bar*> install;  // catenaries: installed after the hinges
    for (Id b : ph.activate_bars) {
      const auto it = bars.find(b);
      if (it == bars.end()) {
        err("barre_inconnue", "activation de la barre inconnue " + id(b));
        continue;
      }
      if (it->second->kind == BarKind::Cable) {
        install.push_back(it->second);
        continue;
      }
      if (!active.insert(b).second) {
        err("barre_deja_active", "barre " + id(b) + " déjà active");
        continue;
      }
      cast_now.insert(b);
      ++degree[it->second->ni];
      ++degree[it->second->nj];
    }
    for (const auto& [b, age] : ph.activate_age)
      if (!cast_now.count(b)) warn("age_sans_effet", "âge donné à la barre " + id(b) + " que la phase n'active pas");
    for (Id b : ph.activate_rel)
      if (!cast_now.count(b)) warn("prolongement_sans_effet", "MODE REL sur la barre " + id(b) + " que la phase n'active pas");
    for (Id l : ph.activate_links) {
      if (!elinks.count(l)) err("liaison_inconnue", "liaison élastique inconnue " + id(l));
      else if (!el_on.insert(l).second) err("liaison_deja_active", "liaison élastique " + id(l) + " déjà active");
    }
    for (Id l : ph.activate_rigid_links) {
      if (!rlinks.count(l)) err("liaison_inconnue", "liaison rigide inconnue " + id(l));
      else if (!rl_on.insert(l).second) err("liaison_deja_active", "liaison rigide " + id(l) + " déjà active");
    }
    for (const auto& [b, pi] : ph.activate_parts) {
      const auto it = model.composite.find(b);
      if (it == model.composite.end()) err("barre_non_composite", "barre " + id(b) + " sans parties de section");
      else if (!active.count(b)) err("partie_sur_barre_inactive", "partie " + std::to_string(pi) + " de la barre inactive " + id(b));
      else if (pi < 0 || static_cast<std::size_t>(pi) >= it->second.parts.size())
        err("partie_inconnue", "partie " + std::to_string(pi) + " inexistante sur la barre " + id(b));
      else if (!cast.insert({b, pi}).second)
        err("partie_deja_coulee", "partie " + std::to_string(pi) + " de la barre " + id(b) + " déjà coulée");
    }
    // 3. deactivations
    for (Id b : ph.deactivate_bars) {
      const auto it = bars.find(b);
      if (it == bars.end()) {
        err("barre_inconnue", "désactivation de la barre inconnue " + id(b));
        continue;
      }
      if (!active.count(b)) {
        err("barre_inactive", "désactivation de la barre inactive " + id(b));
        continue;
      }
      if (it->second->soil.any()) err("barre_sur_sol", "désactivation de la barre " + id(b) + " sur sol élastique");
      const auto tb = tendons_on_bar.find(b);
      if (tb != tendons_on_bar.end())
        for (Id t : tb->second)
          if (ts[t].live || ts[t].on_bed) err("cable_vivant", "barre " + id(b) + " retirée sous le câble tendu " + id(t));
      active.erase(b);
      --degree[it->second->ni];
      --degree[it->second->nj];
    }
    for (Id l : ph.deactivate_links)
      if (!el_on.erase(l)) err("liaison_inactive", "désactivation de la liaison inactive " + id(l));
    for (Id n : ph.remove_support_nodes) {
      if (!nodes.count(n)) err("noeud_inconnu", "retrait d'appui au nœud inconnu " + id(n));
      else if (!supports.erase(n)) err("appui_inactif", "retrait de l'appui inactif au nœud " + id(n));
      else --degree[n];
    }
    // 3b. modifications and jacking
    for (const Support& s : ph.modify_supports) {
      const auto it = supports.find(s.node);
      if (!nodes.count(s.node)) err("noeud_inconnu", "modification d'appui au nœud inconnu " + id(s.node));
      else if (it == supports.end()) err("appui_inactif", "modification de l'appui inactif au nœud " + id(s.node));
      else it->second = &s;
    }
    for (Id n : ph.jack_support_nodes) {
      const auto it = supports.find(n);
      if (!nodes.count(n)) err("noeud_inconnu", "vérinage au nœud inconnu " + id(n));
      else if (it == supports.end()) err("appui_inactif", "vérinage de l'appui inactif au nœud " + id(n));
      else if (it->second->rot != 0.0) err("verinage_rot", "vérinage de l'appui ROT au nœud " + id(n));
    }
    // 3c. section changes and phase hinges (the order of the engine)
    for (const auto& [b, sec] : ph.modify_sections) {
      const auto it = bars.find(b);
      if (it == bars.end()) {
        err("barre_inconnue", "changement de section de la barre inconnue " + id(b));
        continue;
      }
      if (!active.count(b)) err("barre_inactive", "changement de section de la barre inactive " + id(b));
      if (!sections.count(sec)) err("section_inconnue", "section inconnue " + id(sec));
      const Material* mt = nullptr;
      for (const Material& x : model.materials)
        if (x.id == it->second->material) mt = &x;
      if (mt && (mt->concrete || mt->concrete_law))
        err("section_sur_beton_fluant", "changement de section de la barre " + id(b) + " d'un béton qui flue (couler une partie composite)");
      if (model.tapered_of(b) || model.composite.count(b))
        err("section_variable_ou_composite", "changement de section de la barre " + id(b) + " à inertie variable ou composite");
    }
    auto hinge = [&](const Phase::HingeAction& h, bool place) {
      const auto it = bars.find(h.bar);
      if (it == bars.end()) {
        err("barre_inconnue", std::string(place ? "rotule" : "clavage") + " sur la barre inconnue " + id(h.bar));
        return;
      }
      if (!active.count(h.bar)) err("barre_inactive", std::string(place ? "rotule" : "clavage") + " sur la barre inactive " + id(h.bar));
      if (h.end != 0 && h.end != 1) {
        err("extremite_invalide", "extrémité " + std::to_string(h.end) + " (0 ou 1) sur la barre " + id(h.bar));
        return;
      }
      Releases& r = rel.try_emplace(h.bar, it->second->releases).first->second;
      bool any = false;
      for (int c = 0; c < 6; ++c) {
        if (!h.dofs[static_cast<std::size_t>(c)]) continue;
        any = true;
        bool& k = r[static_cast<std::size_t>(6 * h.end + c)];
        if (place && k) err("deja_libere", "rotule sur la barre " + id(h.bar) + " déjà libérée");
        if (!place && !k) err("rien_a_claver", "clavage de la barre " + id(h.bar) + " où rien n'est libéré");
        k = place;
      }
      if (!any) err("rotule_vide", "rotule ou clavage sans composante sur la barre " + id(h.bar));
    };
    for (const Phase::HingeAction& h : ph.place_hinges) hinge(h, true);
    for (const Phase::HingeAction& h : ph.close_hinges) hinge(h, false);
    // 3d. catenary cables: installations, then the actions (engine order)
    for (const Bar* b : install) {
      if (!b->cable) err("catenaire_sans_proprietes", "câble caténaire " + id(b->id) + " sans CableProps");
      if (!catenary.insert(b->id).second) err("catenaire_deja_active", "câble caténaire " + id(b->id) + " déjà installé");
      for (Id n : {b->ni, b->nj})
        if (degree[n] <= 0)
          err("catenaire_ancrage_hors_structure", "câble caténaire " + id(b->id) + " ancré au nœud " + id(n) + " hors de la structure active");
    }
    for (const Phase::CableAction& a : ph.cable_actions) {
      const auto it = bars.find(a.bar);
      if (it == bars.end() || it->second->kind != BarKind::Cable) {
        err("catenaire_inconnue", "action sur la barre " + id(a.bar) + " qui n'est pas un câble caténaire");
        continue;
      }
      if (!catenary.count(a.bar)) {
        err("catenaire_inactive", "action sur le câble caténaire " + id(a.bar) + " non installé");
        continue;
      }
      using K = Phase::CableAction::Kind;
      const bool bad = (a.kind == K::Length || a.kind == K::Horizontal || a.kind == K::EndTension) ? !(a.value > 0)
                       : a.kind == K::Thermal                                                      ? !(a.value > -1)
                                                                                                   : !(a.value >= 1 || a.value == 0);
      if (bad) err("catenaire_valeur_invalide", "action sur le câble caténaire " + id(a.bar) + " : valeur " + std::to_string(a.value));
      if (a.kind == K::Rupture) catenary.erase(a.bar);
    }
    if (!catenary.empty()) {
      // the coupled increment refuses these (typed engine errors)
      bool unilateral = !rl_on.empty();
      for (const auto& [n, sp] : supports) {
        if (sp->decol_axis >= 0) unilateral = true;
        for (int c = 0; c < 6; ++c) unilateral |= sp->oneside[static_cast<std::size_t>(c)] != 0;
      }
      if (unilateral)
        err("catenaire_non_supportee", "câbles caténaires avec appuis unilatéraux ou liaisons rigides actifs");
      for (Id c : catenary) {
        const Bar* b = bars.at(c);
        for (Id n : {b->ni, b->nj})
          if (degree[n] <= 0)
            err("catenaire_ancrage_libre", "câble caténaire " + id(c) + " : son ancrage au nœud " + id(n) + " a quitté la structure");
      }
    }
    // second-order regime: the engine's typed refusals
    for (const auto& [n, v] : ph.imperfection)
      if (!nodes.count(n)) err("noeud_inconnu", "imperfection au nœud inconnu " + id(n));
    if (ph.regime == Phase::Regime::SecondOrder) {
      bool unilateral = !rl_on.empty();
      for (const auto& [n, sp] : supports) {
        if (sp->decol_axis >= 0) unilateral = true;
        for (int c = 0; c < 6; ++c) unilateral |= sp->oneside[static_cast<std::size_t>(c)] != 0;
      }
      if (unilateral) err("second_ordre_non_supporte", "phase du second ordre avec appuis unilatéraux ou liaisons rigides actifs");
      if (!catenary.empty()) err("second_ordre_non_supporte", "phase du second ordre avec câbles caténaires posés");
      if (!model.shells.empty()) err("second_ordre_non_supporte", "phase du second ordre sur un modèle à coques");
      for (Id c : ph.cases)
        if (cases.count(c))
          for (const Load& l : model.load_case(c).loads)
            if (l.kind == Load::Kind::SupportDisp)
              err("second_ordre_non_supporte", "phase du second ordre : déplacement imposé du cas " + id(c));
    }
    for (Id t : ph.inject) {
      if (!tendons.count(t)) err("cable_inconnu", "injection du câble inconnu " + id(t));
      else if (!ts[t].unbonded) warn("injection_sans_effet", "câble " + id(t) + " injecté sans avoir été tendu non injecté");
      else ts[t].unbonded = false;
    }
    for (const Phase::PartialGrout& pg : ph.inject_partial) {
      if (!tendons.count(pg.tendon)) err("cable_inconnu", "injection partielle du câble inconnu " + id(pg.tendon));
      else if (!ts[pg.tendon].unbonded)
        err("injection_partielle_sans_cable_libre", "injection partielle du câble " + id(pg.tendon) + " qui n'est pas tendu non injecté");
      if (pg.s1 < pg.s0) err("intervalle_invalide", "injection partielle du câble " + id(pg.tendon) + " : s1 < s0");
    }
    // 4. loads
    for (Id c : ph.cases) {
      if (!cases.count(c)) {
        err("cas_inconnu", "cas de charge inconnu " + id(c));
        continue;
      }
      for (const Load& l : model.load_case(c).loads) {
        using K = Load::Kind;
        const bool on_bar = l.kind == K::BarUdl || l.kind == K::BarLinear || l.kind == K::BarParabolic || l.kind == K::BarPoint ||
                            l.kind == K::BarThermal || l.kind == K::Prestrain || l.kind == K::PrestrainPoint ||
                            l.kind == K::PrestrainLinear || l.kind == K::PrestrainParabolic;
        if (on_bar && !active.count(l.bar))
          warn("charge_sur_barre_inactive", "cas " + id(c) + " : charge sur la barre inactive " + id(l.bar) + " (sans effet)");
        if ((l.kind == K::NodalForce || l.kind == K::SupportDisp) && degree[l.node] <= 0)
          warn("charge_hors_structure", "cas " + id(c) + " : charge au nœud " + id(l.node) + " hors de la structure active");
      }
    }
    // 5. tensioning
    for (const TensionAction& ta : ph.tension) {
      if (ta.shell) {
        if (!shell_t.count(ta.tendon)) err("cable_inconnu", "câble de coque inconnu " + id(ta.tendon));
        else if (ta.mode != TensionAction::Mode::Jack) err("mode_coque", "câble de coque " + id(ta.tendon) + " : banc non permis");
        continue;
      }
      const auto it = tendons.find(ta.tendon);
      if (it == tendons.end()) {
        err("cable_inconnu", "mise en tension du câble inconnu " + id(ta.tendon));
        continue;
      }
      TState& s = ts[ta.tendon];
      using M = TensionAction::Mode;
      if (ta.mode == M::Bed) {
        if (s.live || s.on_bed) err("banc_deja_tendu", "TENDRE_BANC du câble " + id(ta.tendon) + " déjà tendu");
        s.on_bed = true;
        s.applied = 1.0;
        continue;
      }
      // the engine requires the whole chain active from here on (release,
      // tensioning and detensioning alike)
      for (Id b : it->second->bars)
        if (!active.count(b)) err("cable_sur_barre_inactive", "câble " + id(ta.tendon) + " sur la barre inactive " + id(b));
      if (ta.mode == M::Release) {
        if (!s.on_bed) err("relachement_sans_banc", "RELACHER_BANC du câble " + id(ta.tendon) + " qui n'est pas sur banc");
        s.on_bed = false;
        s.live = true;
      } else if (ta.fraction < s.applied || (ta.fraction < 0 && s.applied <= 0)) {
        if (s.applied <= 0) err("detension_sans_tension", "détension du câble " + id(ta.tendon) + " jamais tendu");
        s.applied = std::max(0.0, ta.fraction);
        s.live = s.applied > 0;
      } else if (ta.fraction == s.applied) {
        warn("tension_sans_effet", "câble " + id(ta.tendon) + " déjà à la fraction " + std::to_string(ta.fraction));
      } else {
        s.applied = ta.fraction;
        s.live = true;
        if (ta.unbonded) s.unbonded = true;
      }
    }
  }
  return out;
}

}  // namespace keta
