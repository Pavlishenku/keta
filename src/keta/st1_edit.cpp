#include "keta/st1_edit.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace keta::st1 {

namespace {

std::string fit(std::string s, int w) {
  if (static_cast<int>(s.size()) > w) return std::string(static_cast<std::size_t>(w), '*');
  return std::string(static_cast<std::size_t>(w) - s.size(), ' ') + s;
}

// Drop the optional leading zero of "0.xxx" / "-0.xxx" when s exceeds w.
std::string squeeze(std::string s, int w) {
  if (static_cast<int>(s.size()) <= w) return s;
  std::size_t p = s[0] == '-' ? 1 : 0;
  if (s.compare(p, 2, "0.") == 0) s.erase(p, 1);
  return s;
}

}  // namespace

std::string fortran_i(long long v, int w) { return fit(std::to_string(v), w); }

std::string fortran_f(double x, int w, int d) {
  if (!std::isfinite(x)) return std::string(static_cast<std::size_t>(w), '*');
  d = std::clamp(d, 0, 30);
  char buf[128];
  std::snprintf(buf, sizeof buf, "%#.*f", d, x == 0 ? 0.0 : x);
  return fit(squeeze(buf, w), w);
}

std::string fortran_g(double x, int w, int d) {
  if (!std::isfinite(x)) return std::string(static_cast<std::size_t>(w), '*');
  const int fw = w - 4;
  d = std::clamp(d, 1, 30);
  char buf[128];
  if (x == 0) {
    std::snprintf(buf, sizeof buf, "%#.*f", d - 1, 0.0);
    return fit(buf, fw) + "    ";
  }
  // d significant digits: the decade k of the ROUNDED value decides F or E.
  std::snprintf(buf, sizeof buf, "%.*e", d - 1, std::fabs(x));
  const char* e = std::strchr(buf, 'e');
  const int e10 = std::atoi(e + 1);
  const int k = e10 + 1;  // 10^(k-1) <= |x| < 10^k
  if (k >= 0 && k <= d) {
    char f[128];
    std::snprintf(f, sizeof f, "%#.*f", d - k, x);
    return fit(squeeze(f, fw), fw) + "    ";
  }
  // E edit: sign 0.ddddd E+ee
  std::string digits(1, buf[0]);
  digits.append(buf + 2, static_cast<std::size_t>(d - 1));
  char ex[16];
  const int ee = k;
  if (std::abs(ee) <= 99) std::snprintf(ex, sizeof ex, "E%c%02d", ee < 0 ? '-' : '+', std::abs(ee));
  else std::snprintf(ex, sizeof ex, "%c%03d", ee < 0 ? '-' : '+', std::abs(ee));
  std::string s = std::string(x < 0 ? "-" : "") + "0." + digits + ex;
  return fit(squeeze(s, w), w);
}

}  // namespace keta::st1
