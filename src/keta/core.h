// keta core — shared types, errors, small math helpers.
// SI base units everywhere: m, N, Pa, kg, rad, s. Global frame: right-handed XYZ, Z up.
// Element local frame: x along the bar (node i -> j), y/z principal axes of the section.
#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace keta {

using Id = std::int64_t;
using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;
using Vec6 = Eigen::Matrix<double, 6, 1>;
using Mat6 = Eigen::Matrix<double, 6, 6>;  // one node's 6 DOFs (support stiffness, nodal mass)
using Vec12 = Eigen::Matrix<double, 12, 1>;
using Mat12 = Eigen::Matrix<double, 12, 12>;
using VecX = Eigen::VectorXd;
using MatX = Eigen::MatrixXd;
using SpMat = Eigen::SparseMatrix<double>;  // CSC
using Triplet = Eigen::Triplet<double>;

inline constexpr double kPi = 3.14159265358979323846;
inline constexpr double kGravity = 9.80665;  // m/s^2

// Every model/analysis failure is a typed error, never a silent NaN.
class Error : public std::runtime_error {
 public:
  explicit Error(std::string msg) : std::runtime_error(std::move(msg)) {}
};

[[noreturn]] inline void fail(const std::string& msg) { throw Error(msg); }
inline void require(bool cond, const std::string& msg) {
  if (!cond) fail(msg);
}

inline bool finite(double v) { return std::isfinite(v); }
inline bool finite(const Vec3& v) { return v.allFinite(); }

// Per-node DOFs, fixed order: ux uy uz rx ry rz.
inline constexpr int kDofPerNode = 6;

}  // namespace keta
