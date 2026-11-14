// Headless renderer for the ST1 graphics mode (manual v24, ch. 16-17, p243-258).
//
// ST1's `DESS` mode is an interactive screen plotter: the user picks a screen
// size (`ECRAN`), a plotter driver (`TRA`), a view direction (`VUE`), a zoom
// (`ZOOM`), a selection (`SELEC`), calibrates the diagram scales (`CENTR`) and
// then issues drawing commands (`GEOM`, `NUM`, `DEPLA`, `EFFORT`, `CONTR`,
// `CONTRG`, `PRESS`) whose result can be dumped to a file (`SAUV`).
//
// keta has no window system, so the same command vocabulary is reproduced as a
// pure function of state: every drawing command appends primitives to a display
// list held by `DrawState`, and the list is rendered to a string by
// `render_svg()` / `render_dxf()`.  Nothing here reads a clock, a locale or an
// environment variable: the same state always renders byte-identical output.
//
// --------------------------------------------------------------------------
// Projection (`VUE ux uy uz`, ST1 p246)
// --------------------------------------------------------------------------
// The view vector `v` points FROM the scene TOWARD the observer (it is the
// screen out-of-plane normal), so `VUE 0 0 1` looks straight down (plan view)
// as in the manual's example, and the ST1 default `VUE 1 1 1` is the classic
// isometric.  The orthonormal screen basis is built as
//
//     v̂     = v / |v|                        (out of the screen)
//     up     = ẑ − (ẑ·v̂) v̂,  normalized      (global Z projected on the screen)
//     right  = up × v̂                        (unit, since up ⊥ v̂)
//
// so (right, up, v̂) is right-handed and `right × up = v̂`.  When v̂ is parallel
// to global Z the projection of ẑ vanishes; the fallback is `up = +Y`, which
// gives a plan view with global X to the right and Y up, matching ST1.
// A zero view vector is not an error: ST1 substitutes (1,1,1) (p246), and so do
// we.  A non-finite view vector IS a typed error.
//
// A world point p projects to screen coordinates s = (right·p, up·p) [m], and
// the screen coordinates are mapped to device pixels by the fitted window
// (see `ViewWindow`), y growing DOWNWARD in device space (SVG convention).
//
// --------------------------------------------------------------------------
// Scales
// --------------------------------------------------------------------------
// `CENTR` / `CENTR GEOM` (p248) fit the projected structure into the canvas
// with a 5 % margin on each side.  `ZOOM cf (xc yc)` (p246) then magnifies the
// device image by `cf` about the screen point (xc, yc) given in NORMALIZED
// screen coordinates measured from the bottom-left corner (ST1's screen frame),
// default (0.5, 0.5) = canvas centre.  `ZOOM 2` therefore doubles every drawn
// distance, i.e. halves the world window visible on screen.
//
// `CENTR <EFFORT|DEPLA|CONTR|CONTRG|PRESS> (MAX xmax)` (p248) calibrates a
// diagram family so that its largest ordinate measures `xmax` CENTIMETRES on
// screen.  Centimetres are converted with `kPixelsPerCm` (100 dpi nominal,
// documented below) since a headless canvas has no physical size.  When a
// family has never been calibrated its diagram amplitude defaults to 10 % of
// the fitted structure extent, as in a first ST1 drawing.
//
// --------------------------------------------------------------------------
// Plotter drivers (`TRA nom`, p244) and buffer files (`SAUV`, p245)
// --------------------------------------------------------------------------
// ST1's drivers emit PostScript (`post`, `postp`, `post3`, `post3p` -> buffer
// `pos.st1`), DXF (`dxf` -> buffer `dxf.st1`) or Calcomp (`calc`).  keta
// SUBSTITUTES SVG for every PostScript/Calcomp driver — both are vector page
// languages, SVG is the one every modern viewer reads — and keeps a real DXF
// R12 writer for `dxf`.  The buffer names follow ST1 with the extension of what
// is actually written: `pos.svg` (postscript drivers), `cal.svg` (calc),
// `dxf.dxf` (dxf).  Paper format (A4/A3, landscape/portrait) only changes the
// canvas aspect ST1 previews with; here it is recorded but the canvas stays the
// `ECRAN` preset, so `post` and `post3` render identically.
//
// --------------------------------------------------------------------------
// Sign conventions
// --------------------------------------------------------------------------
// Diagram ordinates are plotted PERPENDICULAR to the projected bar, on the side
// `n = (d_y, −d_x)` of the travel direction d = normalize(P(j) − P(i)) in device
// pixels; for a bar drawn left to right that is "up" on screen for a positive
// value.  Values keep the keta sign convention throughout (N > 0 in tension,
// stresses tension-positive) — ST1 lists N compression-positive, and that flip
// happens at the ST1 import boundary only, never here.
#pragma once
#include <array>
#include <functional>
#include <string>
#include <vector>

#include "keta/envelope.h"
#include "keta/statics.h"

namespace keta {

// ---------------------------------------------------------------- screen ----

// `ECRAN nom` presets (p244).
enum class Screen { W800_600, W1024_768, W1280_1024 };

struct Canvas {
  int width = 1024;
  int height = 768;  // pixels
};

Canvas canvas_of(Screen s);
Screen screen_from_name(const std::string& name);  // "w800_600" ... (typed error)

// Device resolution used to honour `CENTR ... MAX xmax` (xmax in centimetres):
// 100 dpi nominal => 100/2.54 pixels per centimetre.
inline constexpr double kPixelsPerCm = 100.0 / 2.54;

// --------------------------------------------------------------- plotter ----

// `TRA nom` drivers (p244).  See the header note on the PostScript -> SVG
// substitution.
enum class Plotter { Post, PostP, Post3, Post3P, Dxf, Calc };

Plotter plotter_from_name(const std::string& name);
// Buffer file name written by `SAUV` for this driver: "pos.svg" / "cal.svg" /
// "dxf.dxf" (ST1 writes pos.st1 / dxf.st1).
std::string plotter_buffer_name(Plotter p);
bool plotter_is_dxf(Plotter p);

// --------------------------------------------------------------- palette ----

// The 16 ST1 screen colours (p249).  RGB values are keta's choice (ST1 maps
// them to a driver-dependent pen table); see `color_rgb()` in draw.cpp for the
// documented mapping.
enum class Color {
  Noir,
  Indigo,
  Vert,
  Cyan,
  Rouge,
  Magenta,
  Marron,
  Neutre,
  Gris,
  Bleu,
  Pistache,
  Ciel,
  Corail,
  Mauve,
  Jaune,
  Blanc
};

struct Rgb {
  unsigned char r = 0, g = 0, b = 0;
};

Color color_from_name(const std::string& name);  // typed error on an unknown keyword
const char* color_name(Color c);
Rgb color_rgb(Color c);
int color_aci(Color c);  // AutoCAD colour index used by the DXF writer

// ---------------------------------------------------------------- family ----

// A drawing family: the unit of colour (`COUL`, p249) and of diagram scale
// (`CENTR <family> MAX`, p248).  Geom/Num carry a colour but no scale.
enum class Family { Geom, Num, Effort, Depla, Contr, Contrg, Press };
inline constexpr int kFamilyCount = 7;

const char* family_name(Family f);  // also the DXF layer name (upper case)
Family family_from_name(const std::string& name);

// ------------------------------------------------------------- selection ----

// `SELEC BARRE liste` / `SELEC NOEUD liste` / `SELEC TOUT` (p247).
// Default: the whole structure.  When only nodes are selected, a bar is drawn
// iff BOTH its ends are selected (ST1 draws what is fully visible).
struct Selection {
  bool all = true;
  std::vector<Id> bars;
  std::vector<Id> nodes;

  bool has_bar(const Bar& b) const;
  static Selection everything() { return Selection{}; }
  static Selection of_bars(std::vector<Id> ids);
  static Selection of_nodes(std::vector<Id> ids);
};

// ---------------------------------------------------------- display list ----

using Pt = Eigen::Vector2d;  // device pixels, origin top-left, y downward

enum class PrimKind { Segment, Polyline, Text };

struct Prim {
  PrimKind kind = PrimKind::Polyline;
  Family family = Family::Geom;
  Color color = Color::Noir;
  std::vector<Pt> pts;   // Segment: 2 points; Polyline: >= 2; Text: 1 anchor
  std::string text;      // Text only
  double text_size = 0;  // Text only [px]
};

// ------------------------------------------------------------ projection ----

// Orthonormal screen basis built from a view direction (see the header note).
struct ViewBasis {
  Vec3 right = Vec3(1, 0, 0);
  Vec3 up = Vec3(0, 0, 1);
  Vec3 normal = Vec3(0, 1, 0);  // v̂, out of the screen

  Pt project(const Vec3& p) const { return Pt(right.dot(p), up.dot(p)); }
};
ViewBasis view_basis(const Vec3& view_dir);

// Fitted projected-world -> device-pixel mapping (before the zoom).
struct ViewWindow {
  bool valid = false;
  Pt center = Pt::Zero();  // projected world point placed at the canvas centre
  double scale = 1.0;      // pixels per metre
  double extent_px = 0;    // fitted structure extent (max of width/height) [px]
};

// ------------------------------------------------------------ draw state ----

struct DrawState {
  Screen screen = Screen::W1024_768;
  Plotter plotter = Plotter::Post;
  bool monochrome = true;  // ST1: "les pilotes ci-dessus sont monochromes par défaut"

  Vec3 view = Vec3(1, 1, 1);  // `VUE`, default isometric
  double zoom = 1.0;          // `ZOOM cf`
  double zoom_cx = 0.5;       // `ZOOM cf xc yc`, normalized, from the bottom-left
  double zoom_cy = 0.5;

  Selection sel;
  ViewWindow window;

  // Per-family current colour (`COUL`) and diagram calibration (`CENTR MAX`).
  std::array<Color, kFamilyCount> color = default_colors();
  std::array<double, kFamilyCount> max_cm{};  // 0 = never calibrated -> 10 % default

  double text_size = 11.0;  // px, node/bar numbering
  // Support glyph size, as a fraction of min(canvas).  Glyphs and texts are
  // SCREEN-sized: they keep their pixel size under `ZOOM`, as on an ST1 screen
  // (only the projected geometry follows the zoom homothety).
  double glyph_size = 0.018;

  std::vector<Prim> items;  // the display list

  Canvas canvas() const { return canvas_of(screen); }
  ViewBasis basis() const { return view_basis(view); }
  // Device pixel of a world point (window must be valid).
  Pt to_device(const Vec3& p) const;
  Pt device_of_screen(const Pt& s) const;  // projected-world -> device

  static std::array<Color, kFamilyCount> default_colors();
};

// -------------------------------------------------------------- commands ----

void ecran(DrawState& st, Screen s);                 // `ECRAN nom`
void tra(DrawState& st, Plotter p, bool color);      // `TRA nom (COULEUR|MONOCHROME)`
void vue(DrawState& st, const Vec3& dir);            // `VUE ux uy uz`
void zoom(DrawState& st, double cf);                 // `ZOOM cf`
void zoom(DrawState& st, double cf, double xc, double yc);  // `ZOOM cf xc yc`
void selec(DrawState& st, const Selection& sel);     // `SELEC ...`
void cls(DrawState& st);                             // `CLS` — clears the display list
void coul(DrawState& st, Family f, Color c);         // `COUL <family> couleur`

// `CENTR` (whole model, resets the zoom) and `CENTR GEOM` (current selection).
void centr(DrawState& st, const Model& m);
void centr_geom(DrawState& st, const Model& m);
// `CENTR <family> MAX xmax` — largest ordinate of this family = xmax cm on screen.
void centr_max(DrawState& st, Family f, double xmax_cm);

// ---------------------------------------------------------------- drawing ---

// `GEOM` (p249): selected bars as projected polylines + support glyphs.
void draw_geom(DrawState& st, const Model& m);

// `NUM (NOEUD) (BARRE)` (p249): ids at node positions / bar midpoints.
struct NumOptions {
  bool nodes = true;
  bool bars = true;
};
void draw_num(DrawState& st, const Model& m, const NumOptions& opt = {});

// `DEPLA <(DX..RZ)> <CHARG|COMB|ENV|SURCH> liste` (p250).
//  * `All`, `DX`, `DY`, `DZ` draw the deformed shape (all translations, or the
//    single global component), one polyline per bar, from the node displacements
//    and the station deflections rotated back to the global frame.
//  * `RX`, `RY`, `RZ` draw a perpendicular diagram of the nodal rotation; keta's
//    `Station` carries no rotation field, so that diagram is LINEAR between the
//    two end rotations (documented limitation).
enum class DispComp { All, DX, DY, DZ, RX, RY, RZ };
DispComp disp_comp_from_name(const std::string& name);
void draw_depla(DrawState& st, const Model& m, const CaseResult& r, DispComp comp = DispComp::All);
void draw_depla(DrawState& st, const Model& m, const Envelope& env,
                DispComp comp = DispComp::All);  // ENV: two curves (min + max)

// `EFFORT <N,TY,TZ,MX,MY,MZ> ...` (p251).
enum class ForceComp { N, TY, TZ, MX, MY, MZ };
ForceComp force_comp_from_name(const std::string& name);
void draw_effort(DrawState& st, const Model& m, const CaseResult& r, ForceComp comp);
void draw_effort(DrawState& st, const Model& m, const Envelope& env, ForceComp comp);

// `CONTR <VY,WY,VZ,WZ,VYVZ,VYWZ,WYVZ,WYWZ> ...` (p252): normal stress at an
// extreme fibre.  ST1 names fibres by the CARA distances vy/wy/vz/wz; keta's
// sections do not store them, so the fibre ordinates are taken from the extreme
// stress points of the section (vy = max y, wy = min y, vz = max z, wz = min z)
// and sigma = N/A + My z/Iy − Mz y/Iz (tension positive, keta convention).
enum class StressFiber { VY, WY, VZ, WZ, VYVZ, VYWZ, WYVZ, WYWZ };
StressFiber fiber_from_name(const std::string& name);
void draw_contr(DrawState& st, const Model& m, const CaseResult& r, StressFiber fib);
void draw_contr(DrawState& st, const Model& m, const Envelope& env, StressFiber fib);

// `PRESS <PX,PY,PZ,MX> ...` (p254): soil pressure of bars on an elastic bed.
// Bars without a Winkler bed carry no pressure and are skipped.
enum class PressComp { PX, PY, PZ, MX };
PressComp press_comp_from_name(const std::string& name);
void draw_press(DrawState& st, const Model& m, const CaseResult& r, PressComp comp);
void draw_press(DrawState& st, const Model& m, const Envelope& env, PressComp comp);

// `CONTRG i ...` (p253): generalized stress i.  ST1 evaluates a user-defined
// generalized stress; keta has no such vocabulary in the kernel, so the caller
// (the ST1 layer, a script) supplies the evaluated per-station scalar.
using ScalarField = std::function<double(const Model&, const BarResult&, std::size_t station)>;
void draw_contrg(DrawState& st, const Model& m, const CaseResult& r, const ScalarField& value);
void draw_contrg(DrawState& st, const Model& m, const Envelope& env, const ScalarField& value);

// Generic diagram driver behind EFFORT/CONTR/CONTRG/PRESS: one polyline per bar
// per source, all sources sharing one scale.  Exposed for callers with their own
// per-station quantity.
void draw_bar_diagram(DrawState& st, const Model& m, const std::vector<const CaseResult*>& sources,
                      Family fam, const ScalarField& value);

// ---------------------------------------------------------------- output ----

std::string render_svg(const DrawState& st);
std::string render_dxf(const DrawState& st);
// Dispatch on the current plotter driver (dxf -> DXF, everything else -> SVG).
std::string render(const DrawState& st);

// `SAUV`: write the current display list to `plotter_buffer_name(st.plotter)`
// inside `dir`; returns the path written.
std::string sauv(const DrawState& st, const std::string& dir = ".");
// `TRA 'nomfichier.ext'` + `SAUV`: write to an explicit path with the current driver.
void save_drawing(const DrawState& st, const std::string& path);

}  // namespace keta
