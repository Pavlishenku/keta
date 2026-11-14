// European hot-rolled I/H profile catalog: IPE 80-600, HEA / HEB / HEM 100-1000.
//
// Nominal dimensions (h, b, tw, tf, root radius r) per EN 10365:2017 (formerly
// Euronorm 19-57 for IPE and 53-62 for HE), stored in metres. The profiles are
// returned as finalized `ShapeI` sections WITH their root fillets, so A, Iy,
// Iz, Wel and Wpl reproduce the tabulated values (fillet spandrels exact), J
// follows El Darwish & Johnston (1-3 % of the tables) and Avz the EN 1993-1-1
// rolled-section rule. A fillet-free outline would read A ~3-4 % low and J
// 10-30 % low.
//
// Lookup is tolerant of case, spaces, dashes and underscores: "IPE300",
// "ipe 300" and "IPE-300" are the same profile. `normalize_designation` is the
// canonical key and is shared with every consumer that compares section names
// (the Select by Section component).
#pragma once
#include <optional>
#include <string>
#include <vector>

#include "keta/sections.h"

namespace keta {

// Uppercase and strip whitespace, '-' and '_': " ipe 300 " -> "IPE300".
std::string normalize_designation(const std::string& name);

// The finalized section of a catalog profile (name = canonical designation,
// e.g. "IPE300"), or nullopt for an unknown family / size.
std::optional<Section> catalog_section(const std::string& designation);

// Every catalog designation, family by family in ascending size.
std::vector<std::string> catalog_designations();

}  // namespace keta
