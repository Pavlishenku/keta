// keta::param — a human-readable description of the change between two
// document snapshots. Used by an undo history panel ("+2 Slider",
// "Slider: 6 → 8", "moved 3 nodes").
//
// APPROXIMATE BY DESIGN: entries come from a JSON diff, not from the edit
// sites. Pure and exception-free: unreadable input yields {"edit"};
// byte-identical input yields {}.
//
//   1. node-type multisets: "+N <Name>" / "−N <Name>" per type, adds before
//      removes, each sorted by name. <Name> is the type id's last segment,
//      '_' -> ' ', first letter upper ("keta.select_by_section" ->
//      "Select by section");
//   2. canvas groups and notes, always: "+1 group 'Deck'", "−1 group",
//      "group renamed" / "group edited", "+1 note", "note edited",
//      "moved N notes";
//   3. only when the type multisets match: per key, "<Name>: old → new" when
//      exactly one primitive param field changed (numbers trimmed, strings
//      quoted), else "<Name>: params changed"; position-only changes count
//      toward "moved N node(s)" (listed last);
//   4. only when the whole key -> type map matches: "+N wire(s)" /
//      "−N wire(s)";
//   5. more than 6 entries: the first 5 plus "… and N more".
//
// keta's JSON numbers are doubles, so 1 and 1.0 compare equal: a param edit
// between them is not reported.
#pragma once
#include <string>
#include <vector>

namespace keta::param {

std::vector<std::string> describe_change(const std::string& old_json, const std::string& new_json);

}  // namespace keta::param
