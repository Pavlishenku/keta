// ST1 edition: the RESU tables of an ST1 session written as ST1 itself prints
// them (the resu_*_non_pagine.txt layout), so a keta result file
// reads, diffs and compares line for line against an ST1 one.
//
// ST1 is a Fortran program and its tables are Fortran edit descriptors:
//
//   block header   ' --- RESULTATS DU CHARGEMENT :' I9 ' ---'
//   node table     I11, 4X, n(1X, G11.5)            (Dx Dy Dz Rx Ry Rz)
//   support table  I11 (support), I11 (node), [1X,A3 Min/Max | 4X], n(1X, G11.5)
//   bar tables     I11, F11.3, 1X, A3 (Min/Max or blank), n G13.5
//                  (the Max line of a pair starts with 22 blanks)
//   cable table    I7, n G14.5                       (X Y Z Tension Var.Tension)
//
// with gfortran's G editing: |x| in [0.1, 10^d) after rounding to d digits ->
// F(w-4).(d-k) followed by 4 blanks, otherwise E(w).d with the leading zero
// dropped when the sign does not fit ('-.16903E-02'); zero -> F(w-4).(d-1).
//
// The SIGN CONVENTIONS are ST1's (established on the twelve reference decks):
// N is compression-positive, the shears and moments Ty, Tz, Mx, My, Mz have the
// sign opposite to keta's station convention, displacements and reactions are
// keta's own (global axes).
#pragma once
#include <string>

namespace keta::st1 {

// Fortran Gw.d edit of x (gfortran rules, see the header note).
std::string fortran_g(double x, int w, int d);
// Fortran Fw.d / Iw edits (right-justified; asterisks when it does not fit).
std::string fortran_f(double x, int w, int d);
std::string fortran_i(long long v, int w);

}  // namespace keta::st1
