#include "keta/draw.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <unordered_map>
#include <unordered_set>

#include "keta/element.h"  // bar_rotation (local axes for the deformed shape)

namespace keta {

namespace {

// ---------------------------------------------------------------------------
// Locale-independent fixed-point formatting.  std::to_string / snprintf follow
// the C locale (a comma decimal separator would silently corrupt every SVG and
// DXF number), so numbers are rendered with integer arithmetic.  Non-finite
// values are a typed error, never a "nan" token in the output.
// ---------------------------------------------------------------------------
std::string fmt_fixed(double v, int nd) {
  require(std::isfinite(v), "draw: non-finite coordinate/value in the display list");
  static const double kPow[] = {1, 10, 100, 1000, 10000, 100000, 1000000};
  double p = kPow[nd];
  bool neg = v < 0;
  long long scaled = std::llround(std::fabs(v) * p);
  if (scaled == 0) neg = false;
  long long ip = scaled / static_cast<long long>(p);
  long long fp = scaled % static_cast<long long>(p);
  std::string s;
  if (neg) s += '-';
  s += std::to_string(ip);
  if (nd > 0) {
    std::string f = std::to_string(fp);
    s += '.';
    s.append(static_cast<std::size_t>(nd) - f.size(), '0');
    s += f;
  }
  return s;
}
std::string fmt3(double v) { return fmt_fixed(v, 3); }

std::string escape_xml(const std::string& in) {
  std::string out;
  out.reserve(in.size());
  for (char c : in) {
    switch (c) {
      case '&': out += "&amp;"; break;
      case '<': out += "&lt;"; break;
      case '>': out += "&gt;"; break;
      case '"': out += "&quot;"; break;
      case '\'': out += "&apos;"; break;
      default: out += c;
    }
  }
  return out;
}

// DXF is a line-oriented (group code, value) stream: a newline inside a value
// would desynchronize every reader.
std::string one_line(const std::string& in) {
  std::string out = in;
  for (char& c : out)
    if (c == '\n' || c == '\r') c = ' ';
  return out;
}

bool contains(const std::vector<Id>& v, Id id) {
  return std::find(v.begin(), v.end(), id) != v.end();
}

int family_index(Family f) { return static_cast<int>(f); }

// The set of nodes a command acts on: with a bar selection, ST1 acts on the
// nodes of the selected bars too (they are what is visible on screen).
std::unordered_set<Id> selected_nodes(const DrawState& st, const Model& m) {
  std::unordered_set<Id> out;
  if (st.sel.all) {
    for (const Node& n : m.nodes) out.insert(n.id);
    return out;
  }
  for (Id n : st.sel.nodes) out.insert(n);
  const std::unordered_set<Id> sel_bars(st.sel.bars.begin(), st.sel.bars.end());
  for (const Bar& b : m.bars)
    if (sel_bars.count(b.id)) {
      out.insert(b.ni);
      out.insert(b.nj);
    }
  return out;
}

std::vector<const Bar*> selected_bars(const DrawState& st, const Model& m, const char* what) {
  std::vector<const Bar*> out;
  // Selection::has_bar over every bar, with the id lists hashed once (a
  // linear membership test per bar was O(bars x selection))
  const Selection& sel = st.sel;
  const std::unordered_set<Id> sel_bars(sel.bars.begin(), sel.bars.end());
  const std::unordered_set<Id> sel_nodes(sel.nodes.begin(), sel.nodes.end());
  auto has_bar = [&](const Bar& b) {
    if (sel.all) return true;
    if (!sel.bars.empty()) return sel_bars.count(b.id) > 0;
    if (!sel.nodes.empty()) return sel_nodes.count(b.ni) > 0 && sel_nodes.count(b.nj) > 0;
    return false;
  };
  for (const Bar& b : m.bars)
    if (has_bar(b)) out.push_back(&b);
  require(!out.empty(), std::string("draw: ") + what + " — the selection contains no bar");
  return out;
}

// Index BarResult by bar id (a CaseResult need not be in model order).
std::unordered_map<Id, const BarResult*> index_bars(const CaseResult& r) {
  std::unordered_map<Id, const BarResult*> ix;
  ix.reserve(r.bars.size());
  for (const BarResult& br : r.bars) ix.emplace(br.bar, &br);
  return ix;
}

void push_polyline(DrawState& st, Family fam, std::vector<Pt> pts) {
  require(pts.size() >= 2, "draw: a polyline needs at least two points");
  for (const Pt& p : pts) require(p.allFinite(), "draw: non-finite device coordinate");
  Prim pr;
  pr.kind = PrimKind::Polyline;
  pr.family = fam;
  pr.color = st.color[static_cast<std::size_t>(family_index(fam))];
  pr.pts = std::move(pts);
  st.items.push_back(std::move(pr));
}

void push_segment(DrawState& st, Family fam, const Pt& a, const Pt& b) {
  require(a.allFinite() && b.allFinite(), "draw: non-finite device coordinate");
  Prim pr;
  pr.kind = PrimKind::Segment;
  pr.family = fam;
  pr.color = st.color[static_cast<std::size_t>(family_index(fam))];
  pr.pts = {a, b};
  st.items.push_back(std::move(pr));
}

void push_text(DrawState& st, Family fam, const Pt& at, const std::string& text, double size) {
  require(at.allFinite(), "draw: non-finite device coordinate");
  Prim pr;
  pr.kind = PrimKind::Text;
  pr.family = fam;
  pr.color = st.color[static_cast<std::size_t>(family_index(fam))];
  pr.pts = {at};
  pr.text = text;
  pr.text_size = size;
  st.items.push_back(std::move(pr));
}

// Fit the projected bounding box of `nodes` into the canvas with a 5 % margin.
constexpr double kMargin = 0.05;

void fit(DrawState& st, const Model& m, const std::unordered_set<Id>& nodes) {
  require(!nodes.empty(), "draw: CENTR — nothing to centre (empty selection)");
  ViewBasis vb = st.basis();
  double x0 = 0, x1 = 0, y0 = 0, y1 = 0;
  bool first = true;
  for (const Node& n : m.nodes) {
    if (!nodes.count(n.id)) continue;
    require(finite(n.x), "draw: non-finite node coordinate");
    Pt p = vb.project(n.x);
    if (first) {
      x0 = x1 = p.x();
      y0 = y1 = p.y();
      first = false;
    } else {
      x0 = std::min(x0, p.x());
      x1 = std::max(x1, p.x());
      y0 = std::min(y0, p.y());
      y1 = std::max(y1, p.y());
    }
  }
  require(!first, "draw: CENTR — nothing to centre (no node in the selection)");
  Canvas c = canvas_of(st.screen);
  double dx = x1 - x0, dy = y1 - y0;
  double usable_w = c.width * (1 - 2 * kMargin), usable_h = c.height * (1 - 2 * kMargin);
  double sx = dx > 0 ? usable_w / dx : std::numeric_limits<double>::infinity();
  double sy = dy > 0 ? usable_h / dy : std::numeric_limits<double>::infinity();
  double s = std::min(sx, sy);
  // A structure that projects to a single point (or a line seen end-on) has no
  // finite fit: fall back to a unit scale so the drawing still exists.
  if (!std::isfinite(s) || s <= 0) s = 1.0;
  st.window.valid = true;
  st.window.center = Pt(0.5 * (x0 + x1), 0.5 * (y0 + y1));
  st.window.scale = s;
  st.window.extent_px = std::max(dx, dy) * s;
}

void ensure_window(DrawState& st, const Model& m) {
  // "Lors d'un premier dessin, le centrage s'effectue par défaut sur la
  // géométrie de la structure" (ST1 p248).
  if (!st.window.valid) fit(st, m, selected_nodes(st, m));
}

// Diagram amplitude in device pixels for a family: the `CENTR ... MAX xmax`
// calibration when given, else 10 % of the fitted structure extent.
double amplitude_px(const DrawState& st, Family fam) {
  double cm = st.max_cm[static_cast<std::size_t>(family_index(fam))];
  if (cm > 0) return cm * kPixelsPerCm;
  double a = 0.10 * st.window.extent_px;
  return a > 0 ? a : 0.10 * std::min(canvas_of(st.screen).width, canvas_of(st.screen).height);
}

// Unit normal to the projected bar, in device pixels: the side positive values
// are plotted on.  A bar seen end-on (zero projected length) has no side; the
// caller skips it.
bool bar_normal(const Pt& pi, const Pt& pj, Pt& n) {
  Pt d = pj - pi;
  double L = d.norm();
  if (!(L > 1e-9)) return false;
  d /= L;
  n = Pt(d.y(), -d.x());
  return true;
}

}  // namespace

// ------------------------------------------------------------------ screen --

Canvas canvas_of(Screen s) {
  switch (s) {
    case Screen::W800_600: return {800, 600};
    case Screen::W1024_768: return {1024, 768};
    case Screen::W1280_1024: return {1280, 1024};
  }
  fail("draw: unknown screen preset");
}

Screen screen_from_name(const std::string& name) {
  std::string k;
  for (char c : name) k += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  if (k == "w800_600") return Screen::W800_600;
  if (k == "w1024_768") return Screen::W1024_768;
  if (k == "w1280_1024") return Screen::W1280_1024;
  fail("draw: unknown ECRAN preset '" + name + "' (w800_600, w1024_768, w1280_1024)");
}

// ----------------------------------------------------------------- plotter --

Plotter plotter_from_name(const std::string& name) {
  std::string k;
  for (char c : name) k += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  if (k == "post") return Plotter::Post;
  if (k == "postp") return Plotter::PostP;
  if (k == "post3") return Plotter::Post3;
  if (k == "post3p") return Plotter::Post3P;
  if (k == "dxf") return Plotter::Dxf;
  if (k == "calc") return Plotter::Calc;
  fail("draw: unknown TRA driver '" + name + "' (post, postp, post3, post3p, dxf, calc)");
}

bool plotter_is_dxf(Plotter p) { return p == Plotter::Dxf; }

std::string plotter_buffer_name(Plotter p) {
  // ST1 writes pos.st1 (PostScript) / dxf.st1 (DXF); keta substitutes SVG for
  // the PostScript and Calcomp drivers and keeps the extension honest.
  switch (p) {
    case Plotter::Post:
    case Plotter::PostP:
    case Plotter::Post3:
    case Plotter::Post3P: return "pos.svg";
    case Plotter::Calc: return "cal.svg";
    case Plotter::Dxf: return "dxf.dxf";
  }
  fail("draw: unknown plotter driver");
}

// ----------------------------------------------------------------- palette --

// ST1 names 16 screen colours (p249) and maps them, per plotter, to a pen table
// held in the parameter file.  keta has no parameter file, so the mapping to RGB
// is fixed here.  The choice: the eight "primary" ST1 names take their obvious
// saturated web value; the six mixed names take the shade their French name
// describes (marron = brown, neutre = mid grey-beige, pistache = pistachio
// green, ciel = sky blue, corail = coral, mauve = mauve).  Blanc is pure white:
// on keta's white sheet a Blanc pen is invisible, exactly as ST1's Noir pen is
// on its black screen.
Rgb color_rgb(Color c) {
  switch (c) {
    case Color::Noir: return {0, 0, 0};
    case Color::Indigo: return {75, 0, 130};
    case Color::Vert: return {0, 128, 0};
    case Color::Cyan: return {0, 176, 176};
    case Color::Rouge: return {204, 0, 0};
    case Color::Magenta: return {200, 0, 200};
    case Color::Marron: return {128, 64, 0};
    case Color::Neutre: return {160, 150, 130};
    case Color::Gris: return {128, 128, 128};
    case Color::Bleu: return {0, 0, 204};
    case Color::Pistache: return {147, 197, 114};
    case Color::Ciel: return {135, 206, 235};
    case Color::Corail: return {255, 127, 80};
    case Color::Mauve: return {180, 120, 220};
    case Color::Jaune: return {230, 200, 0};
    case Color::Blanc: return {255, 255, 255};
  }
  fail("draw: unknown colour");
}

// AutoCAD Colour Index nearest each palette entry (DXF group 62).
int color_aci(Color c) {
  switch (c) {
    case Color::Noir: return 250;      // near-black (ACI 7 flips with the DXF background)
    case Color::Indigo: return 200;
    case Color::Vert: return 3;
    case Color::Cyan: return 4;
    case Color::Rouge: return 1;
    case Color::Magenta: return 6;
    case Color::Marron: return 34;
    case Color::Neutre: return 43;
    case Color::Gris: return 8;
    case Color::Bleu: return 5;
    case Color::Pistache: return 82;
    case Color::Ciel: return 141;
    case Color::Corail: return 22;
    case Color::Mauve: return 211;
    case Color::Jaune: return 2;
    case Color::Blanc: return 255;
  }
  fail("draw: unknown colour");
}

const char* color_name(Color c) {
  switch (c) {
    case Color::Noir: return "NOIR";
    case Color::Indigo: return "INDIGO";
    case Color::Vert: return "VERT";
    case Color::Cyan: return "CYAN";
    case Color::Rouge: return "ROUGE";
    case Color::Magenta: return "MAGENTA";
    case Color::Marron: return "MARRON";
    case Color::Neutre: return "NEUTRE";
    case Color::Gris: return "GRIS";
    case Color::Bleu: return "BLEU";
    case Color::Pistache: return "PISTACHE";
    case Color::Ciel: return "CIEL";
    case Color::Corail: return "CORAIL";
    case Color::Mauve: return "MAUVE";
    case Color::Jaune: return "JAUNE";
    case Color::Blanc: return "BLANC";
  }
  fail("draw: unknown colour");
}

Color color_from_name(const std::string& name) {
  std::string k;
  for (char c : name) k += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  for (int i = 0; i < 16; ++i) {
    Color c = static_cast<Color>(i);
    if (k == color_name(c)) return c;
  }
  fail("draw: unknown colour keyword '" + name +
       "' (NOIR INDIGO VERT CYAN ROUGE MAGENTA MARRON NEUTRE GRIS BLEU PISTACHE CIEL CORAIL "
       "MAUVE JAUNE BLANC)");
}

// ------------------------------------------------------------------ family --

const char* family_name(Family f) {
  switch (f) {
    case Family::Geom: return "GEOM";
    case Family::Num: return "NUM";
    case Family::Effort: return "EFFORT";
    case Family::Depla: return "DEPLA";
    case Family::Contr: return "CONTR";
    case Family::Contrg: return "CONTRG";
    case Family::Press: return "PRESS";
  }
  fail("draw: unknown family");
}

Family family_from_name(const std::string& name) {
  std::string k;
  for (char c : name) k += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  for (int i = 0; i < kFamilyCount; ++i) {
    Family f = static_cast<Family>(i);
    if (k == family_name(f)) return f;
  }
  fail("draw: unknown drawing family '" + name +
       "' (GEOM NUM EFFORT DEPLA CONTR CONTRG PRESS)");
}

// --------------------------------------------------------------- selection --

Selection Selection::of_bars(std::vector<Id> ids) {
  Selection s;
  s.all = false;
  s.bars = std::move(ids);
  return s;
}
Selection Selection::of_nodes(std::vector<Id> ids) {
  Selection s;
  s.all = false;
  s.nodes = std::move(ids);
  return s;
}
bool Selection::has_bar(const Bar& b) const {
  if (all) return true;
  if (!bars.empty()) return contains(bars, b.id);
  if (!nodes.empty()) return contains(nodes, b.ni) && contains(nodes, b.nj);
  return false;
}

// -------------------------------------------------------------- projection --

ViewBasis view_basis(const Vec3& view_dir) {
  require(finite(view_dir), "draw: VUE — non-finite view direction");
  Vec3 v = view_dir;
  // ST1 p246: "si on donne ux=uy=uz=0, ST1 prendra automatiquement ux=uy=uz=1".
  if (v.norm() < 1e-12) v = Vec3(1, 1, 1);
  v.normalize();
  Vec3 z(0, 0, 1);
  Vec3 up = z - z.dot(v) * v;
  if (up.norm() < 1e-9) {
    // Looking straight down (or up) the global Z axis: plan view, +Y up on screen.
    up = Vec3(0, 1, 0);
    up = up - up.dot(v) * v;
  }
  up.normalize();
  ViewBasis vb;
  vb.normal = v;
  vb.up = up;
  vb.right = up.cross(v);
  vb.right.normalize();
  return vb;
}

std::array<Color, kFamilyCount> DrawState::default_colors() {
  std::array<Color, kFamilyCount> c{};
  c[static_cast<std::size_t>(Family::Geom)] = Color::Noir;
  c[static_cast<std::size_t>(Family::Num)] = Color::Gris;
  c[static_cast<std::size_t>(Family::Effort)] = Color::Rouge;
  c[static_cast<std::size_t>(Family::Depla)] = Color::Vert;
  c[static_cast<std::size_t>(Family::Contr)] = Color::Bleu;
  c[static_cast<std::size_t>(Family::Contrg)] = Color::Magenta;
  c[static_cast<std::size_t>(Family::Press)] = Color::Marron;
  return c;
}

Pt DrawState::device_of_screen(const Pt& s) const {
  require(window.valid, "draw: no view window — call centr()/centr_geom() or draw the geometry");
  Canvas c = canvas_of(screen);
  double W = c.width, H = c.height;
  Pt d0(W * 0.5 + window.scale * (s.x() - window.center.x()),
        H * 0.5 - window.scale * (s.y() - window.center.y()));
  // ZOOM is a homothety about a normalized screen point measured from the
  // BOTTOM-left corner (ST1 screen frame); device y grows downward.
  Pt z(zoom_cx * W, (1.0 - zoom_cy) * H);
  return z + zoom * (d0 - z);
}

Pt DrawState::to_device(const Vec3& p) const { return device_of_screen(basis().project(p)); }

// ---------------------------------------------------------------- commands --

void ecran(DrawState& st, Screen s) {
  st.screen = s;
  st.window.valid = false;  // the canvas changed: the fit has to be redone
}

void tra(DrawState& st, Plotter p, bool color) {
  st.plotter = p;
  st.monochrome = !color;
}

void vue(DrawState& st, const Vec3& dir) {
  (void)view_basis(dir);  // validate now, not at the first projection
  st.view = dir;
  st.window.valid = false;
}

void zoom(DrawState& st, double cf) { zoom(st, cf, 0.5, 0.5); }

void zoom(DrawState& st, double cf, double xc, double yc) {
  require(std::isfinite(cf) && cf > 0, "draw: ZOOM — the factor must be finite and positive");
  require(std::isfinite(xc) && std::isfinite(yc), "draw: ZOOM — non-finite centre");
  st.zoom = cf;
  st.zoom_cx = xc;
  st.zoom_cy = yc;
}

void selec(DrawState& st, const Selection& sel) { st.sel = sel; }

void cls(DrawState& st) { st.items.clear(); }

void coul(DrawState& st, Family f, Color c) {
  st.color[static_cast<std::size_t>(family_index(f))] = c;
}

void centr(DrawState& st, const Model& m) {
  st.zoom = 1.0;
  st.zoom_cx = st.zoom_cy = 0.5;
  std::unordered_set<Id> all;
  for (const Node& n : m.nodes) all.insert(n.id);
  fit(st, m, all);
}

void centr_geom(DrawState& st, const Model& m) { fit(st, m, selected_nodes(st, m)); }

void centr_max(DrawState& st, Family f, double xmax_cm) {
  require(std::isfinite(xmax_cm) && xmax_cm > 0,
          "draw: CENTR MAX — xmax must be finite and positive [cm]");
  st.max_cm[static_cast<std::size_t>(family_index(f))] = xmax_cm;
}

// ----------------------------------------------------------------- drawing --

void draw_geom(DrawState& st, const Model& m) {
  m.index();
  ensure_window(st, m);
  std::vector<const Bar*> bars = selected_bars(st, m, "GEOM");
  ViewBasis vb = st.basis();
  for (const Bar* b : bars) {
    Pt pi = st.device_of_screen(vb.project(m.node(b->ni).x));
    Pt pj = st.device_of_screen(vb.project(m.node(b->nj).x));
    push_polyline(st, Family::Geom, {pi, pj});
  }
  // Support glyphs (ST1 draws them with the structure): a triangle whose apex
  // sits on the node marks blocked TRANSLATION, hatch marks under its base mark
  // blocked ROTATION.  A purely elastic support (springs, nothing rigidly
  // blocked) gets the hatched base line without the triangle.
  Canvas c = canvas_of(st.screen);
  double g = st.glyph_size * std::min(c.width, c.height);
  std::unordered_set<Id> nodes = selected_nodes(st, m);
  for (const Support& s : m.supports) {
    if (!s.any() || !nodes.count(s.node)) continue;
    Pt p = st.device_of_screen(vb.project(m.node(s.node).x));
    bool trans = s.fixed[0] || s.fixed[1] || s.fixed[2];
    bool rot = s.fixed[3] || s.fixed[4] || s.fixed[5];
    Pt bl(p.x() - g, p.y() + 1.6 * g), br(p.x() + g, p.y() + 1.6 * g);
    if (trans) {
      push_segment(st, Family::Geom, p, bl);
      push_segment(st, Family::Geom, p, br);
      push_segment(st, Family::Geom, bl, br);
    } else {
      push_segment(st, Family::Geom, bl, br);
    }
    if (rot || !trans) {
      for (int k = 0; k < 3; ++k) {
        double x = bl.x() + (k + 0.5) * (2 * g / 3);
        push_segment(st, Family::Geom, Pt(x, bl.y()), Pt(x - 0.4 * g, bl.y() + 0.5 * g));
      }
    }
  }
}

void draw_num(DrawState& st, const Model& m, const NumOptions& opt) {
  m.index();
  ensure_window(st, m);
  ViewBasis vb = st.basis();
  if (opt.nodes) {
    std::unordered_set<Id> nodes = selected_nodes(st, m);
    require(!nodes.empty(), "draw: NUM — the selection contains no node");
    // Deterministic order: model order, not hash order.
    for (const Node& n : m.nodes) {
      if (!nodes.count(n.id)) continue;
      Pt p = st.device_of_screen(vb.project(n.x));
      push_text(st, Family::Num, Pt(p.x() + 4, p.y() - 4), std::to_string(n.id), st.text_size);
    }
  }
  if (opt.bars) {
    std::vector<const Bar*> bars = selected_bars(st, m, "NUM");
    for (const Bar* b : bars) {
      Pt pi = st.device_of_screen(vb.project(m.node(b->ni).x));
      Pt pj = st.device_of_screen(vb.project(m.node(b->nj).x));
      Pt mid = 0.5 * (pi + pj);
      push_text(st, Family::Num, Pt(mid.x() + 4, mid.y() - 4), std::to_string(b->id),
                st.text_size);
    }
  }
}

// ------------------------------------------------------------- diagrams -----

void draw_bar_diagram(DrawState& st, const Model& m, const std::vector<const CaseResult*>& sources,
                      Family fam, const ScalarField& value) {
  require(!sources.empty(), "draw: a diagram needs at least one result source");
  require(static_cast<bool>(value), "draw: a diagram needs a station value function");
  m.index();
  ensure_window(st, m);
  std::vector<const Bar*> bars = selected_bars(st, m, family_name(fam));
  ViewBasis vb = st.basis();

  std::vector<std::unordered_map<Id, const BarResult*>> ix;
  ix.reserve(sources.size());
  for (const CaseResult* r : sources) {
    require(r != nullptr, "draw: null result source");
    ix.push_back(index_bars(*r));
  }

  // One pass to find the largest ordinate: every source shares one scale, so an
  // envelope's two curves stay comparable.
  double vmax = 0;
  for (std::size_t s = 0; s < sources.size(); ++s)
    for (const Bar* b : bars) {
      auto it = ix[s].find(b->id);
      if (it == ix[s].end()) continue;
      const BarResult& br = *it->second;
      for (std::size_t k = 0; k < br.stations.size(); ++k) {
        double v = value(m, br, k);
        require(std::isfinite(v), "draw: non-finite station value in a diagram");
        vmax = std::max(vmax, std::fabs(v));
      }
    }
  double f = vmax > 0 ? amplitude_px(st, fam) / vmax : 0.0;

  for (std::size_t s = 0; s < sources.size(); ++s)
    for (const Bar* b : bars) {
      auto it = ix[s].find(b->id);
      if (it == ix[s].end()) continue;
      const BarResult& br = *it->second;
      if (br.stations.size() < 2) continue;
      Pt pi = st.device_of_screen(vb.project(m.node(b->ni).x));
      Pt pj = st.device_of_screen(vb.project(m.node(b->nj).x));
      Pt n;
      if (!bar_normal(pi, pj, n)) continue;  // bar seen end-on: no diagram side
      double L = m.bar_chord_length(*b);
      std::vector<Pt> curve;
      curve.reserve(br.stations.size());
      for (std::size_t k = 0; k < br.stations.size(); ++k) {
        double t = L > 0 ? br.stations[k].x / L : 0.0;
        Pt base = pi + t * (pj - pi);
        curve.push_back(base + (f * value(m, br, k)) * n);
      }
      push_polyline(st, fam, curve);
      // Closing baseline: the two end ordinates dropped back onto the bar.
      push_segment(st, fam, pi, curve.front());
      push_segment(st, fam, pj, curve.back());
    }
}

// ------------------------------------------------------------- deformed -----

namespace {

// Global displacement of a bar station: node-line interpolation of the
// undeformed geometry + the local deflection field rotated back to global.
// For an eccentric bar the station line is the element axis, offset from the
// node line; the drawing keeps the node line as the reference (documented
// approximation, the offset is a rigid arm and does not deform).
Vec3 station_disp(const Model& m, const Bar& b, const Station& s) {
  Mat3 R = bar_rotation(m.node(b.ni).x, m.node(b.nj).x, b.roll);
  return R.transpose() * s.defl;
}

void draw_deformed(DrawState& st, const Model& m, const std::vector<const CaseResult*>& sources,
                   DispComp comp) {
  std::vector<const Bar*> bars = selected_bars(st, m, "DEPLA");
  ViewBasis vb = st.basis();
  std::vector<std::unordered_map<Id, const BarResult*>> ix;
  for (const CaseResult* r : sources) {
    require(r != nullptr, "draw: null result source");
    ix.push_back(index_bars(*r));
  }

  auto disp_of = [&](const Bar& b, const Station& s) {
    Vec3 d = station_disp(m, b, s);
    require(finite(d), "draw: non-finite station deflection");
    if (comp == DispComp::All) return d;
    Vec3 o = Vec3::Zero();
    if (comp == DispComp::DX) o.x() = d.x();
    if (comp == DispComp::DY) o.y() = d.y();
    if (comp == DispComp::DZ) o.z() = d.z();
    return o;
  };

  // Amplification: the largest PROJECTED displacement is drawn at the family
  // amplitude (default 10 % of the structure extent, or the CENTR MAX value).
  double dmax = 0;
  for (std::size_t s = 0; s < sources.size(); ++s)
    for (const Bar* b : bars) {
      auto it = ix[s].find(b->id);
      if (it == ix[s].end()) continue;
      for (const Station& stn : it->second->stations)
        dmax = std::max(dmax, vb.project(disp_of(*b, stn)).norm());
    }
  double px_per_m = st.window.scale * st.zoom;
  double amp = amplitude_px(st, Family::Depla);
  double f = dmax > 0 ? amp / (dmax * px_per_m) : 0.0;

  for (std::size_t s = 0; s < sources.size(); ++s)
    for (const Bar* b : bars) {
      auto it = ix[s].find(b->id);
      if (it == ix[s].end()) continue;
      const BarResult& br = *it->second;
      if (br.stations.size() < 2) continue;
      const Vec3& xi = m.node(b->ni).x;
      const Vec3& xj = m.node(b->nj).x;
      double L = m.bar_chord_length(*b);
      std::vector<Pt> curve;
      curve.reserve(br.stations.size());
      for (const Station& stn : br.stations) {
        double t = L > 0 ? stn.x / L : 0.0;
        Vec3 p = xi + t * (xj - xi) + f * disp_of(*b, stn);
        curve.push_back(st.device_of_screen(vb.project(p)));
      }
      push_polyline(st, Family::Depla, curve);
    }
}

// Rotation "diagram": keta's Station has no rotation field, so the curve is
// linear between the two nodal rotations (documented limitation, ST1 plots the
// study sections).
void draw_rotation(DrawState& st, const Model& m, const std::vector<const CaseResult*>& sources,
                   int comp) {
  std::vector<const Bar*> bars = selected_bars(st, m, "DEPLA");
  ViewBasis vb = st.basis();
  double vmax = 0;
  for (const CaseResult* r : sources)
    for (const Bar* b : bars) {
      int i = m.node_index(b->ni), j = m.node_index(b->nj);
      require(i >= 0 && j >= 0, "draw: DEPLA — dangling bar end");
      vmax = std::max({vmax, std::fabs(r->node_disp[static_cast<std::size_t>(i)](3 + comp)),
                       std::fabs(r->node_disp[static_cast<std::size_t>(j)](3 + comp))});
    }
  double f = vmax > 0 ? amplitude_px(st, Family::Depla) / vmax : 0.0;
  for (const CaseResult* r : sources)
    for (const Bar* b : bars) {
      Pt pi = st.device_of_screen(vb.project(m.node(b->ni).x));
      Pt pj = st.device_of_screen(vb.project(m.node(b->nj).x));
      Pt n;
      if (!bar_normal(pi, pj, n)) continue;
      double ri = r->node_disp[static_cast<std::size_t>(m.node_index(b->ni))](3 + comp);
      double rj = r->node_disp[static_cast<std::size_t>(m.node_index(b->nj))](3 + comp);
      std::vector<Pt> curve = {pi + (f * ri) * n, pj + (f * rj) * n};
      push_polyline(st, Family::Depla, curve);
      push_segment(st, Family::Depla, pi, curve.front());
      push_segment(st, Family::Depla, pj, curve.back());
    }
}

void depla_impl(DrawState& st, const Model& m, const std::vector<const CaseResult*>& sources,
                DispComp comp) {
  require(!sources.empty(), "draw: DEPLA needs at least one result source");
  m.index();
  ensure_window(st, m);
  switch (comp) {
    case DispComp::RX: draw_rotation(st, m, sources, 0); return;
    case DispComp::RY: draw_rotation(st, m, sources, 1); return;
    case DispComp::RZ: draw_rotation(st, m, sources, 2); return;
    default: draw_deformed(st, m, sources, comp);
  }
}

}  // namespace

DispComp disp_comp_from_name(const std::string& name) {
  std::string k;
  for (char c : name) k += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  if (k.empty() || k == "ALL" || k == "TOUT") return DispComp::All;
  if (k == "DX") return DispComp::DX;
  if (k == "DY") return DispComp::DY;
  if (k == "DZ") return DispComp::DZ;
  if (k == "RX") return DispComp::RX;
  if (k == "RY") return DispComp::RY;
  if (k == "RZ") return DispComp::RZ;
  fail("draw: unknown DEPLA component '" + name + "' (DX DY DZ RX RY RZ)");
}

void draw_depla(DrawState& st, const Model& m, const CaseResult& r, DispComp comp) {
  depla_impl(st, m, {&r}, comp);
}
void draw_depla(DrawState& st, const Model& m, const Envelope& env, DispComp comp) {
  depla_impl(st, m, {&env.lo.value, &env.hi.value}, comp);
}

// --------------------------------------------------------------- EFFORT -----

ForceComp force_comp_from_name(const std::string& name) {
  std::string k;
  for (char c : name) k += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  if (k == "N") return ForceComp::N;
  if (k == "TY") return ForceComp::TY;
  if (k == "TZ") return ForceComp::TZ;
  if (k == "MX") return ForceComp::MX;
  if (k == "MY") return ForceComp::MY;
  if (k == "MZ") return ForceComp::MZ;
  fail("draw: unknown EFFORT component '" + name + "' (N TY TZ MX MY MZ)");
}

namespace {
ScalarField force_field(ForceComp c) {
  return [c](const Model&, const BarResult& br, std::size_t k) {
    const Station& s = br.stations[k];
    switch (c) {
      case ForceComp::N: return s.N;
      case ForceComp::TY: return s.Vy;
      case ForceComp::TZ: return s.Vz;
      case ForceComp::MX: return s.T;
      case ForceComp::MY: return s.My;
      case ForceComp::MZ: return s.Mz;
    }
    return 0.0;
  };
}
}  // namespace

void draw_effort(DrawState& st, const Model& m, const CaseResult& r, ForceComp comp) {
  draw_bar_diagram(st, m, {&r}, Family::Effort, force_field(comp));
}
void draw_effort(DrawState& st, const Model& m, const Envelope& env, ForceComp comp) {
  draw_bar_diagram(st, m, {&env.lo.value, &env.hi.value}, Family::Effort, force_field(comp));
}

// ---------------------------------------------------------------- CONTR -----

StressFiber fiber_from_name(const std::string& name) {
  std::string k;
  for (char c : name) k += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  if (k == "VY") return StressFiber::VY;
  if (k == "WY") return StressFiber::WY;
  if (k == "VZ") return StressFiber::VZ;
  if (k == "WZ") return StressFiber::WZ;
  if (k == "VYVZ") return StressFiber::VYVZ;
  if (k == "VYWZ") return StressFiber::VYWZ;
  if (k == "WYVZ") return StressFiber::WYVZ;
  if (k == "WYWZ") return StressFiber::WYWZ;
  fail("draw: unknown CONTR fibre '" + name + "' (VY WY VZ WZ VYVZ VYWZ WYVZ WYWZ)");
}

namespace {

// ST1's fibre names are the CARA extreme-fibre distances vy/wy (along local y)
// and vz/wz (along local z).  keta stores no such distances, so they are read
// off the section's generated stress points: vy = max y, wy = min y, vz = max z,
// wz = min z.  A section with no stress point (generic properties given by
// CARA) has no fibre geometry — that is a typed rejection, never a zero stress.
void fiber_ordinates(const Section& sec, StressFiber f, double& y, double& z) {
  require(!sec.stress_points.empty(),
          "draw: CONTR — section '" + sec.name +
              "' has no stress point (generic properties): no fibre ordinate to plot");
  double ymin = sec.stress_points[0].y, ymax = ymin;
  double zmin = sec.stress_points[0].z, zmax = zmin;
  for (const StressPoint& p : sec.stress_points) {
    ymin = std::min(ymin, p.y);
    ymax = std::max(ymax, p.y);
    zmin = std::min(zmin, p.z);
    zmax = std::max(zmax, p.z);
  }
  switch (f) {
    case StressFiber::VY: y = ymax; z = 0; return;
    case StressFiber::WY: y = ymin; z = 0; return;
    case StressFiber::VZ: y = 0; z = zmax; return;
    case StressFiber::WZ: y = 0; z = zmin; return;
    case StressFiber::VYVZ: y = ymax; z = zmax; return;
    case StressFiber::VYWZ: y = ymax; z = zmin; return;
    case StressFiber::WYVZ: y = ymin; z = zmax; return;
    case StressFiber::WYWZ: y = ymin; z = zmin; return;
  }
  fail("draw: unknown CONTR fibre");
}

ScalarField fiber_field(StressFiber f) {
  return [f](const Model& m, const BarResult& br, std::size_t k) {
    const Bar& b = m.bar(br.bar);
    const Section& sec = m.section(b.section);
    double y = 0, z = 0;
    fiber_ordinates(sec, f, y, z);
    const SectionProps& p = sec.props;
    const Station& s = br.stations[k];
    // keta convention: sigma = N/A + My z/Iy - Mz y/Iz, TENSION positive.
    double sigma = s.N / p.A;
    if (p.Iy > 0) sigma += s.My * z / p.Iy;
    if (p.Iz > 0) sigma -= s.Mz * y / p.Iz;
    return sigma;
  };
}

}  // namespace

void draw_contr(DrawState& st, const Model& m, const CaseResult& r, StressFiber fib) {
  draw_bar_diagram(st, m, {&r}, Family::Contr, fiber_field(fib));
}
void draw_contr(DrawState& st, const Model& m, const Envelope& env, StressFiber fib) {
  draw_bar_diagram(st, m, {&env.lo.value, &env.hi.value}, Family::Contr, fiber_field(fib));
}

// ---------------------------------------------------------------- PRESS -----

PressComp press_comp_from_name(const std::string& name) {
  std::string k;
  for (char c : name) k += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  if (k == "PX") return PressComp::PX;
  if (k == "PY") return PressComp::PY;
  if (k == "PZ") return PressComp::PZ;
  if (k == "MX") return PressComp::MX;
  fail("draw: unknown PRESS component '" + name + "' (PX PY PZ MX)");
}

namespace {
ScalarField press_field(PressComp c) {
  return [c](const Model&, const BarResult& br, std::size_t k) {
    // Bars without a Winkler bed carry no pressure vector: flat zero curve.
    if (c == PressComp::MX) return k < br.soil_m.size() ? br.soil_m[k] : 0.0;
    if (k >= br.soil_p.size()) return 0.0;
    const Vec3& p = br.soil_p[k];
    switch (c) {
      case PressComp::PX: return p.x();
      case PressComp::PY: return p.y();
      case PressComp::PZ: return p.z();
      default: return 0.0;
    }
  };
}
}  // namespace

void draw_press(DrawState& st, const Model& m, const CaseResult& r, PressComp comp) {
  draw_bar_diagram(st, m, {&r}, Family::Press, press_field(comp));
}
void draw_press(DrawState& st, const Model& m, const Envelope& env, PressComp comp) {
  draw_bar_diagram(st, m, {&env.lo.value, &env.hi.value}, Family::Press, press_field(comp));
}

// --------------------------------------------------------------- CONTRG -----

void draw_contrg(DrawState& st, const Model& m, const CaseResult& r, const ScalarField& value) {
  draw_bar_diagram(st, m, {&r}, Family::Contrg, value);
}
void draw_contrg(DrawState& st, const Model& m, const Envelope& env, const ScalarField& value) {
  draw_bar_diagram(st, m, {&env.lo.value, &env.hi.value}, Family::Contrg, value);
}

// ------------------------------------------------------------------ SVG -----

namespace {
std::string hex(Rgb c) {
  static const char* d = "0123456789abcdef";
  std::string s = "#";
  for (unsigned char v : {c.r, c.g, c.b}) {
    s += d[v >> 4];
    s += d[v & 0xF];
  }
  return s;
}
std::string stroke_of(const DrawState& st, const Prim& p) {
  return st.monochrome ? std::string("#000000") : hex(color_rgb(p.color));
}
}  // namespace

std::string render_svg(const DrawState& st) {
  Canvas c = canvas_of(st.screen);
  std::string w = std::to_string(c.width), h = std::to_string(c.height);
  std::string out;
  out.reserve(4096 + 64 * st.items.size());
  out += "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" + w + "\" height=\"" + h +
         "\" viewBox=\"0 0 " + w + " " + h + "\">\n";
  out += "<rect x=\"0\" y=\"0\" width=\"" + w + "\" height=\"" + h +
         "\" fill=\"#ffffff\" stroke=\"none\"/>\n";
  out += "<g fill=\"none\" stroke-width=\"1\" stroke-linecap=\"round\" stroke-linejoin=\"round\">\n";
  for (const Prim& p : st.items) {
    std::string col = stroke_of(st, p);
    switch (p.kind) {
      case PrimKind::Segment:
        out += "<line x1=\"" + fmt3(p.pts[0].x()) + "\" y1=\"" + fmt3(p.pts[0].y()) + "\" x2=\"" +
               fmt3(p.pts[1].x()) + "\" y2=\"" + fmt3(p.pts[1].y()) + "\" stroke=\"" + col +
               "\"/>\n";
        break;
      case PrimKind::Polyline: {
        out += "<polyline points=\"";
        for (std::size_t i = 0; i < p.pts.size(); ++i) {
          if (i) out += ' ';
          out += fmt3(p.pts[i].x());
          out += ',';
          out += fmt3(p.pts[i].y());
        }
        out += "\" stroke=\"" + col + "\"/>\n";
        break;
      }
      case PrimKind::Text:
        out += "<text x=\"" + fmt3(p.pts[0].x()) + "\" y=\"" + fmt3(p.pts[0].y()) +
               "\" font-family=\"monospace\" font-size=\"" + fmt3(p.text_size) + "\" fill=\"" +
               col + "\" stroke=\"none\">" + escape_xml(p.text) + "</text>\n";
        break;
    }
  }
  out += "</g>\n</svg>\n";
  return out;
}

// ------------------------------------------------------------------ DXF -----

namespace {
void dxf_pair(std::string& out, int code, const std::string& value) {
  out += std::to_string(code);
  out += '\n';
  out += value;
  out += '\n';
}
void dxf_line(std::string& out, const DrawState& st, const Prim& p, const Pt& a, const Pt& b,
              double H) {
  dxf_pair(out, 0, "LINE");
  dxf_pair(out, 8, family_name(p.family));
  dxf_pair(out, 62, std::to_string(st.monochrome ? 7 : color_aci(p.color)));
  dxf_pair(out, 10, fmt_fixed(a.x(), 4));
  dxf_pair(out, 20, fmt_fixed(H - a.y(), 4));
  dxf_pair(out, 30, "0.0");
  dxf_pair(out, 11, fmt_fixed(b.x(), 4));
  dxf_pair(out, 21, fmt_fixed(H - b.y(), 4));
  dxf_pair(out, 31, "0.0");
}
}  // namespace

std::string render_dxf(const DrawState& st) {
  // Minimal DXF R12 (AC1009): HEADER with the version + drawing extents, then
  // ENTITIES (LINE / TEXT, one layer per drawing family), then EOF.  No TABLES
  // section: R12 readers accept entities on implicitly created layers, and a
  // half-populated LAYER/LTYPE table would be worse than none.
  // The device frame is y-DOWN (SVG); DXF is y-UP, so y_dxf = H - y_device.
  Canvas c = canvas_of(st.screen);
  double H = c.height;
  std::string out;
  out.reserve(4096 + 256 * st.items.size());
  dxf_pair(out, 0, "SECTION");
  dxf_pair(out, 2, "HEADER");
  dxf_pair(out, 9, "$ACADVER");
  dxf_pair(out, 1, "AC1009");
  dxf_pair(out, 9, "$EXTMIN");
  dxf_pair(out, 10, "0.0");
  dxf_pair(out, 20, "0.0");
  dxf_pair(out, 30, "0.0");
  dxf_pair(out, 9, "$EXTMAX");
  dxf_pair(out, 10, fmt_fixed(c.width, 4));
  dxf_pair(out, 20, fmt_fixed(c.height, 4));
  dxf_pair(out, 30, "0.0");
  dxf_pair(out, 0, "ENDSEC");
  dxf_pair(out, 0, "SECTION");
  dxf_pair(out, 2, "ENTITIES");
  for (const Prim& p : st.items) {
    switch (p.kind) {
      case PrimKind::Segment: dxf_line(out, st, p, p.pts[0], p.pts[1], H); break;
      case PrimKind::Polyline:
        for (std::size_t i = 1; i < p.pts.size(); ++i)
          dxf_line(out, st, p, p.pts[i - 1], p.pts[i], H);
        break;
      case PrimKind::Text:
        dxf_pair(out, 0, "TEXT");
        dxf_pair(out, 8, family_name(p.family));
        dxf_pair(out, 62, std::to_string(st.monochrome ? 7 : color_aci(p.color)));
        dxf_pair(out, 10, fmt_fixed(p.pts[0].x(), 4));
        dxf_pair(out, 20, fmt_fixed(H - p.pts[0].y(), 4));
        dxf_pair(out, 30, "0.0");
        dxf_pair(out, 40, fmt_fixed(p.text_size, 4));
        dxf_pair(out, 1, one_line(p.text));
        break;
    }
  }
  dxf_pair(out, 0, "ENDSEC");
  dxf_pair(out, 0, "EOF");
  return out;
}

std::string render(const DrawState& st) {
  return plotter_is_dxf(st.plotter) ? render_dxf(st) : render_svg(st);
}

void save_drawing(const DrawState& st, const std::string& path) {
  std::string text = render(st);
  std::ofstream f(path, std::ios::binary);
  require(f.good(), "draw: cannot open '" + path + "' for writing");
  f.write(text.data(), static_cast<std::streamsize>(text.size()));
  require(f.good(), "draw: write failed on '" + path + "'");
}

std::string sauv(const DrawState& st, const std::string& dir) {
  std::string path = dir;
  if (!path.empty() && path.back() != '/') path += '/';
  path += plotter_buffer_name(st.plotter);
  save_drawing(st, path);
  return path;
}

}  // namespace keta
