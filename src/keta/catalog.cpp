#include "keta/catalog.h"

#include <cctype>

namespace keta {

namespace {

// One row: designation number and nominal h, b, tw, tf, r in MILLIMETRES (as
// tabulated by EN 10365; converted to metres on lookup).
struct Row {
  int size;
  double h, b, tw, tf, r;
};

const Row kIPE[] = {
    {80, 80, 46, 3.8, 5.2, 5},      {100, 100, 55, 4.1, 5.7, 7},   {120, 120, 64, 4.4, 6.3, 7},
    {140, 140, 73, 4.7, 6.9, 7},    {160, 160, 82, 5.0, 7.4, 9},   {180, 180, 91, 5.3, 8.0, 9},
    {200, 200, 100, 5.6, 8.5, 12},  {220, 220, 110, 5.9, 9.2, 12}, {240, 240, 120, 6.2, 9.8, 15},
    {270, 270, 135, 6.6, 10.2, 15}, {300, 300, 150, 7.1, 10.7, 15}, {330, 330, 160, 7.5, 11.5, 18},
    {360, 360, 170, 8.0, 12.7, 18}, {400, 400, 180, 8.6, 13.5, 21}, {450, 450, 190, 9.4, 14.6, 21},
    {500, 500, 200, 10.2, 16.0, 21}, {550, 550, 210, 11.1, 17.2, 24}, {600, 600, 220, 12.0, 19.0, 24},
};

const Row kHEA[] = {
    {100, 96, 100, 5.0, 8.0, 12},    {120, 114, 120, 5.0, 8.0, 12},   {140, 133, 140, 5.5, 8.5, 12},
    {160, 152, 160, 6.0, 9.0, 15},   {180, 171, 180, 6.0, 9.5, 15},   {200, 190, 200, 6.5, 10.0, 18},
    {220, 210, 220, 7.0, 11.0, 18},  {240, 230, 240, 7.5, 12.0, 21},  {260, 250, 260, 7.5, 12.5, 24},
    {280, 270, 280, 8.0, 13.0, 24},  {300, 290, 300, 8.5, 14.0, 27},  {320, 310, 300, 9.0, 15.5, 27},
    {340, 330, 300, 9.5, 16.5, 27},  {360, 350, 300, 10.0, 17.5, 27}, {400, 390, 300, 11.0, 19.0, 27},
    {450, 440, 300, 11.5, 21.0, 27}, {500, 490, 300, 12.0, 23.0, 27}, {550, 540, 300, 12.5, 24.0, 27},
    {600, 590, 300, 13.0, 25.0, 27}, {650, 640, 300, 13.5, 26.0, 27}, {700, 690, 300, 14.5, 27.0, 27},
    {800, 790, 300, 15.0, 28.0, 30}, {900, 890, 300, 16.0, 30.0, 30}, {1000, 990, 300, 16.5, 31.0, 30},
};

const Row kHEB[] = {
    {100, 100, 100, 6.0, 10.0, 12},    {120, 120, 120, 6.5, 11.0, 12},   {140, 140, 140, 7.0, 12.0, 12},
    {160, 160, 160, 8.0, 13.0, 15},    {180, 180, 180, 8.5, 14.0, 15},   {200, 200, 200, 9.0, 15.0, 18},
    {220, 220, 220, 9.5, 16.0, 18},    {240, 240, 240, 10.0, 17.0, 21},  {260, 260, 260, 10.0, 17.5, 24},
    {280, 280, 280, 10.5, 18.0, 24},   {300, 300, 300, 11.0, 19.0, 27},  {320, 320, 300, 11.5, 20.5, 27},
    {340, 340, 300, 12.0, 21.5, 27},   {360, 360, 300, 12.5, 22.5, 27},  {400, 400, 300, 13.5, 24.0, 27},
    {450, 450, 300, 14.0, 26.0, 27},   {500, 500, 300, 14.5, 28.0, 27},  {550, 550, 300, 15.0, 29.0, 27},
    {600, 600, 300, 15.5, 30.0, 27},   {650, 650, 300, 16.0, 31.0, 27},  {700, 700, 300, 17.0, 32.0, 27},
    {800, 800, 300, 17.5, 33.0, 30},   {900, 900, 300, 18.5, 35.0, 30},  {1000, 1000, 300, 19.0, 36.0, 30},
};

const Row kHEM[] = {
    {100, 120, 106, 12.0, 20.0, 12},   {120, 140, 126, 12.5, 21.0, 12},  {140, 160, 146, 13.0, 22.0, 12},
    {160, 180, 166, 14.0, 23.0, 15},   {180, 200, 186, 14.5, 24.0, 15},  {200, 220, 206, 15.0, 25.0, 18},
    {220, 240, 226, 15.5, 26.0, 18},   {240, 270, 248, 18.0, 32.0, 21},  {260, 290, 268, 18.0, 32.5, 24},
    {280, 310, 288, 18.5, 33.0, 24},   {300, 340, 310, 21.0, 39.0, 27},  {320, 359, 309, 21.0, 40.0, 27},
    {340, 377, 309, 21.0, 40.0, 27},   {360, 395, 308, 21.0, 40.0, 27},  {400, 432, 307, 21.0, 40.0, 27},
    {450, 478, 307, 21.0, 40.0, 27},   {500, 524, 306, 21.0, 40.0, 27},  {550, 572, 306, 21.0, 40.0, 27},
    {600, 620, 305, 21.0, 40.0, 27},   {650, 668, 305, 21.0, 40.0, 27},  {700, 716, 304, 21.0, 40.0, 27},
    {800, 814, 303, 21.0, 40.0, 30},   {900, 910, 302, 21.0, 40.0, 30},  {1000, 1008, 302, 21.0, 40.0, 30},
};

struct Family {
  const char* prefix;
  const Row* rows;
  std::size_t n;
};

const Family kFamilies[] = {
    {"IPE", kIPE, std::size(kIPE)},
    {"HEA", kHEA, std::size(kHEA)},
    {"HEB", kHEB, std::size(kHEB)},
    {"HEM", kHEM, std::size(kHEM)},
};

}  // namespace

std::string normalize_designation(const std::string& name) {
  std::string out;
  out.reserve(name.size());
  for (char c : name) {
    if (std::isspace(static_cast<unsigned char>(c)) || c == '-' || c == '_') continue;
    out.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
  }
  return out;
}

std::optional<Section> catalog_section(const std::string& designation) {
  const std::string norm = normalize_designation(designation);
  for (const Family& f : kFamilies) {
    const std::string prefix = f.prefix;
    if (norm.compare(0, prefix.size(), prefix) != 0) continue;
    const std::string num = norm.substr(prefix.size());
    if (num.empty() || num.size() > 5) return std::nullopt;
    for (char c : num)
      if (!std::isdigit(static_cast<unsigned char>(c))) return std::nullopt;
    const int size = std::stoi(num);
    for (std::size_t i = 0; i < f.n; ++i) {
      const Row& r = f.rows[i];
      if (r.size != size) continue;
      Section s;
      s.name = prefix + std::to_string(r.size);
      s.shape = ShapeI{r.h * 1e-3, r.b * 1e-3, r.tw * 1e-3, r.tf * 1e-3, r.r * 1e-3};
      s.finalize();
      return s;
    }
    return std::nullopt;
  }
  return std::nullopt;
}

std::vector<std::string> catalog_designations() {
  std::vector<std::string> out;
  for (const Family& f : kFamilies)
    for (std::size_t i = 0; i < f.n; ++i) out.push_back(std::string(f.prefix) + std::to_string(f.rows[i].size));
  return out;
}

}  // namespace keta
