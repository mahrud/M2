#include "interface/toric.h"

#include <M2/config.h>
#include <M2/math-include.h>

#include "interrupted.hpp"
#include "interface/gmp-util.h"
#include "m2tbb.hpp"
#include "matrices/matrix-con.hpp"
#include "matrices/matrix.hpp"

#include <gmp.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <libnormaliz/cone.h>

#include <algorithm>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

using Integer = mpz_class;
using Ray = std::vector<Integer>;
using ConeKeys = std::vector<int>;
using FanCones = std::vector<ConeKeys>;

struct ThreadLimitGuard {
  int old;
  bool active;
  explicit ThreadLimitGuard(int n) : old(0), active(n > 0) {
    if (active) old = libnormaliz::set_thread_limit(n);
  }
  ~ThreadLimitGuard() { if (active) libnormaliz::set_thread_limit(old); }
};

void clearNormalizInterrupt() { libnormaliz::nmz_interrupted = 0; }

bool decodeCones(M2_arrayint a, FanCones& result)
{
  if (a == nullptr || a->len < 1) return false;
  size_t p = 0;
  const size_t n = static_cast<size_t>(a->array[p++]);
  result.clear();
  result.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    if (p >= static_cast<size_t>(a->len)) return false;
    const int m = a->array[p++];
    if (m < 0 || p + static_cast<size_t>(m) > static_cast<size_t>(a->len))
      return false;
    ConeKeys c;
    c.reserve(static_cast<size_t>(m));
    for (int j = 0; j < m; ++j) c.push_back(a->array[p++]);
    result.push_back(std::move(c));
  }
  return p == static_cast<size_t>(a->len);
}

bool validFan(const Matrix* input, const FanCones& cones)
{
  if (input == nullptr || input->get_ring() == nullptr) return false;
  const int nrays = static_cast<int>(input->n_rows());
  for (const ConeKeys& c : cones) {
    std::set<int> seen;
    for (int i : c)
      if (i < 0 || i >= nrays || !seen.insert(i).second) return false;
  }
  return true;
}

std::vector<Ray> readRays(const Matrix* input)
{
  std::vector<Ray> rays(input->n_rows(), Ray(input->n_cols()));
  for (size_t i = 0; i < input->n_rows(); ++i)
    for (size_t j = 0; j < input->n_cols(); ++j)
      rays[i][j] = Integer(input->elem(i, j).get_mpz());
  return rays;
}

libnormaliz::Matrix<Integer> normalizMatrix(const std::vector<Ray>& rays,
                                             const ConeKeys& cone)
{
  libnormaliz::Matrix<Integer> result(cone.size(), rays.front().size());
  for (size_t i = 0; i < cone.size(); ++i)
    for (size_t j = 0; j < rays.front().size(); ++j)
      result[i][j] = rays[cone[i]][j];
  return result;
}

std::vector<int> sortedCone(std::vector<int> c)
{
  std::sort(c.begin(), c.end());
  c.erase(std::unique(c.begin(), c.end()), c.end());
  return c;
}

bool contains(const ConeKeys& c, const ConeKeys& face)
{
  return std::includes(c.begin(), c.end(), face.begin(), face.end());
}

struct SmithStatus {
  bool fullRank;
  bool unimodular;
};

SmithStatus smithStatus(const std::vector<Ray>& rays, const ConeKeys& c)
{
  if (c.empty()) return {true, true};
  if (c.size() > rays.front().size()) return {false, false};
  fmpz_mat_t A, S;
  fmpz_mat_init(A, rays.front().size(), c.size());
  fmpz_mat_init(S, rays.front().size(), c.size());
  for (size_t j = 0; j < c.size(); ++j)
    for (size_t i = 0; i < rays.front().size(); ++i)
      fmpz_set_mpz(fmpz_mat_entry(A, i, j), rays[c[j]][i].get_mpz_t());
  fmpz_mat_snf(S, A);
  bool full = true, unit = true;
  for (size_t i = 0; i < c.size(); ++i) {
    if (fmpz_is_zero(fmpz_mat_entry(S, i, i))) full = false;
    if (!fmpz_is_one(fmpz_mat_entry(S, i, i))) unit = false;
  }
  fmpz_mat_clear(A);
  fmpz_mat_clear(S);
  return {full, full && unit};
}

bool matrixRankFull(const std::vector<Ray>& rays, const ConeKeys& c)
{ return smithStatus(rays, c).fullRank; }

bool isSmooth(const std::vector<Ray>& rays, const ConeKeys& c)
{ return smithStatus(rays, c).unimodular; }

bool findIndependentRows(const std::vector<Ray>& rays, const ConeKeys& face,
                         size_t at, std::vector<size_t>& rows)
{
  if (rows.size() == face.size()) {
    fmpz_mat_t A;
    fmpz_mat_init(A, face.size(), face.size());
    for (size_t j = 0; j < face.size(); ++j)
      for (size_t i = 0; i < face.size(); ++i)
        fmpz_set_mpz(fmpz_mat_entry(A, i, j),
                     rays[face[j]][rows[i]].get_mpz_t());
    fmpz_t det;
    fmpz_init(det);
    fmpz_mat_det(det, A);
    const bool independent = !fmpz_is_zero(det);
    fmpz_clear(det);
    fmpz_mat_clear(A);
    return independent;
  }
  const size_t dimension = rays.front().size();
  for (size_t i = at; i < dimension; ++i) {
    rows.push_back(i);
    if (findIndependentRows(rays, face, i + 1, rows)) return true;
    rows.pop_back();
  }
  return false;
}

bool isRelativeInterior(const std::vector<Ray>& rays, const ConeKeys& face,
                        const Ray& candidate)
{
  if (face.empty() || candidate.size() != rays.front().size()) return false;
  std::vector<size_t> rows;
  if (!findIndependentRows(rays, face, 0, rows)) return false;
  fmpz_mat_t A;
  fmpz_mat_init(A, face.size(), face.size());
  for (size_t j = 0; j < face.size(); ++j)
    for (size_t i = 0; i < face.size(); ++i)
      fmpz_set_mpz(fmpz_mat_entry(A, i, j),
                   rays[face[j]][rows[i]].get_mpz_t());
  fmpz_t det, numerator;
  fmpz_init(det);
  fmpz_init(numerator);
  fmpz_mat_det(det, A);
  const int sign = fmpz_sgn(det);
  bool interior = sign != 0;
  for (size_t j = 0; interior && j < face.size(); ++j) {
    fmpz_mat_t replaced;
    fmpz_mat_init_set(replaced, A);
    for (size_t i = 0; i < face.size(); ++i)
      fmpz_set_mpz(fmpz_mat_entry(replaced, i, j),
                   candidate[rows[i]].get_mpz_t());
    fmpz_mat_det(numerator, replaced);
    interior = fmpz_sgn(numerator) == sign;
    fmpz_mat_clear(replaced);
  }
  fmpz_clear(numerator);
  fmpz_clear(det);
  fmpz_mat_clear(A);
  return interior;
}

void combinations(const ConeKeys& c, size_t k, size_t at, ConeKeys& current,
                   std::vector<ConeKeys>& result)
{
  if (current.size() == k) {
    result.push_back(current);
    return;
  }
  for (size_t i = at; i + (k - current.size()) <= c.size(); ++i) {
    current.push_back(c[i]);
    combinations(c, k, i + 1, current, result);
    current.pop_back();
  }
}

bool insertUnique(FanCones& cones, const ConeKeys& c,
                  std::unordered_set<std::string>& seen)
{
  ConeKeys s = sortedCone(c);
  std::string key;
  for (int x : s) { key += std::to_string(x); key += ','; }
  if (!seen.insert(key).second) return false;
  cones.push_back(std::move(s));
  return true;
}

FanCones pullingSubdivision(const std::vector<Ray>& rays, const FanCones& input,
                            int strategy, int seed, int threads, bool verbose)
{
  FanCones output;
  std::unordered_set<std::string> seen;
  std::vector<int> order(rays.size());
  std::iota(order.begin(), order.end(), 0);
  if (seed != 0) {
    std::mt19937 generator(static_cast<unsigned int>(seed));
    std::shuffle(order.begin(), order.end(), generator);
  }
  std::vector<size_t> position(rays.size());
  for (size_t i = 0; i < order.size(); ++i) position[order[i]] = i;
  for (const ConeKeys& sigma : input) {
    if (sigma.size() <= rays.front().size() &&
        matrixRankFull(rays, sigma)) {
      insertUnique(output, sigma, seen);
      continue;
    }
    ConeKeys localOrder = sigma;
    std::sort(localOrder.begin(), localOrder.end(),
              [&](int a, int b) { return position[a] < position[b]; });
    if ((strategy & 1) != 0) std::reverse(localOrder.begin(), localOrder.end());
    auto local = normalizMatrix(rays, localOrder);
    const auto& localGenerators = local.get_elements();
    libnormaliz::Cone<Integer> cone(libnormaliz::Type::cone, local);
    cone.setVerbose(verbose);
    const auto& triangulation = cone.getTriangulation(
        libnormaliz::ConeProperty::PullingTriangulation);
    const auto& generators = triangulation.second.get_elements();
    for (const auto& simplex : triangulation.first) {
      ConeKeys result;
      for (auto key : simplex.key) {
        if (key >= generators.size())
          throw std::runtime_error("invalid Normaliz simplex key");
        auto it = std::find(localGenerators.begin(), localGenerators.end(), generators[key]);
        if (it == localGenerators.end())
          throw std::runtime_error("Normaliz triangulation generator mismatch");
        size_t generatorIndex = static_cast<size_t>(it - localGenerators.begin());
        if (generatorIndex >= localOrder.size())
          throw std::runtime_error("Normaliz triangulation generator index mismatch");
        result.push_back(localOrder[generatorIndex]);
      }
      insertUnique(output, result, seen);
    }
  }
  return output;
}

bool chooseHilbertRay(const std::vector<Ray>& rays, const ConeKeys& face,
                      int strategy, int threads, bool verbose, Ray& answer)
{
  clearNormalizInterrupt();
  try {
    auto local = normalizMatrix(rays, face);
    libnormaliz::Cone<Integer> cone(libnormaliz::Type::cone, local);
    cone.setVerbose(verbose);
    libnormaliz::ConeProperties props(libnormaliz::ConeProperty::HilbertBasis);
    if ((strategy & 1) != 0)
      props.set(libnormaliz::ConeProperty::PrimalMode);
    else
      props.set(libnormaliz::ConeProperty::DualMode);
    cone.compute(props);
    bool found = false;
    Integer bestScore;
    Ray best;
    for (const auto& candidate : cone.getHilbertBasis()) {
      bool old = false;
      for (const Ray& r : local.get_elements())
        if (candidate == r) old = true;
      if (!old && isRelativeInterior(rays, face, candidate)) {
        Integer score = 0;
        for (const Integer& x : candidate) score += abs(x);
        if (!found || score < bestScore) {
          found = true;
          bestScore = score;
          best = candidate;
        }
      }
    }
    if (found) { answer = std::move(best); return true; }
  } catch (const libnormaliz::InterruptException&) {
    ERROR("rawSmoothFan interrupted");
    return false;
  } catch (const std::exception& e) {
    ERROR(e.what());
    return false;
  }
  ERROR("Normaliz found no Hilbert basis ray in the relative interior of the singular face");
  return false;
}

void appendRay(std::vector<Ray>& rays, Ray ray)
{
  Integer g = 0;
  for (const Integer& x : ray) g = gcd(g, abs(x));
  if (g > 1) for (Integer& x : ray) x /= g;
  rays.push_back(std::move(ray));
}

bool smoothRefinement(std::vector<Ray>& rays, FanCones& cones,
                      int strategy, int limit, int threads, bool verbose)
{
  const size_t d = rays.front().size();
  bool changed;
  do {
    changed = false;
    for (size_t k = 2; k <= d; ++k) {
    while (true) {
      if (system_interrupted()) {
        ERROR("rawSmoothFan interrupted");
        return false;
      }
      std::vector<std::vector<ConeKeys>> local(cones.size());
      mtbb::parallel_for(mtbb::blocked_range<int>{0, static_cast<int>(cones.size())},
        [&](const mtbb::blocked_range<int>& range) {
          for (int q = range.begin(); q != range.end(); ++q) {
            if (isSmooth(rays, cones[q])) continue;
            std::vector<ConeKeys> faces;
            ConeKeys current;
            combinations(cones[q], k, 0, current, faces);
            for (const ConeKeys& face : faces)
              if (!isSmooth(rays, face)) local[q].push_back(face);
          }
        });
      if (system_interrupted()) {
        ERROR("rawSmoothFan interrupted");
        return false;
      }

      FanCones batch;
      std::set<ConeKeys> unique;
      for (const auto& faces : local)
        for (const ConeKeys& face : faces) unique.insert(face);
      for (const ConeKeys& face : unique) batch.push_back(face);
      if (batch.empty()) break;

      for (const ConeKeys& face : batch) {
        if (limit > 0 && static_cast<int>(cones.size()) >= limit) return true;
        Ray ray;
        if (!chooseHilbertRay(rays, face, strategy, threads, verbose, ray))
          return false;
        const int newIndex = static_cast<int>(rays.size());
        appendRay(rays, std::move(ray));
        FanCones next;
        std::unordered_set<std::string> seen;
        for (const ConeKeys& sigma : cones) {
          if (!contains(sigma, face)) { insertUnique(next, sigma, seen); continue; }
          for (int removed : face) {
            ConeKeys child;
            for (int x : sigma) if (x != removed) child.push_back(x);
            child.push_back(newIndex);
            insertUnique(next, child, seen);
          }
        }
        cones.swap(next);
        changed = true;
      }
    }
    }
  } while (changed);
  return true;
}

const Matrix* encodeFan(const Matrix* input, const std::vector<Ray>& rays,
                        const FanCones& cones)
{
  const Ring* R = input->get_ring();
  size_t width = 2 + rays.front().size();
  for (const ConeKeys& c : cones) width = std::max(width, 2 + c.size());
  MatrixConstructor mat(R->make_FreeModule(rays.size() + cones.size()), width);
  size_t row = 0;
  auto put = [&](size_t i, int value) {
    mpz_ptr z = newitem(__mpz_struct);
    mpz_init_set_si(z, value);
    mpz_reallocate_limbs(z);
    mat.set_entry(row, i, ring_elem(z));
  };
  auto putBig = [&](size_t i, const Integer& value) {
    mpz_ptr z = newitem(__mpz_struct);
    mpz_init_set(z, value.get_mpz_t());
    mpz_reallocate_limbs(z);
    mat.set_entry(row, i, ring_elem(z));
  };
  for (size_t i = 0; i < rays.size(); ++i) {
    put(0, 0);
    put(1, static_cast<int>(rays[i].size()));
    for (size_t j = 0; j < rays[i].size(); ++j) putBig(2 + j, rays[i][j]);
    ++row;
  }
  for (size_t i = 0; i < cones.size(); ++i) {
    put(0, 1);
    put(1, static_cast<int>(cones[i].size()));
    for (size_t j = 0; j < cones[i].size(); ++j) put(2 + j, cones[i][j]);
    ++row;
  }
  return mat.to_matrix();
}

const Matrix* run(const Matrix* raysInput, M2_arrayint encoded, bool smooth,
                  int strategy, int seed, int limit, int threads, bool verbose)
{
  FanCones cones;
  if (raysInput == nullptr || raysInput->get_ring() != globalZZ) {
    ERROR("fan rays must be a matrix over ZZ");
    return nullptr;
  }
  if (!decodeCones(encoded, cones) || !validFan(raysInput, cones)) {
    ERROR("invalid ray/cone fan encoding");
    return nullptr;
  }
  std::vector<Ray> rays = readRays(raysInput);
  if (rays.empty() || rays.front().empty()) {
    ERROR("fan must have nonempty rays");
    return nullptr;
  }
  try {
    ThreadLimitGuard guard(threads);
    clearNormalizInterrupt();
    if (smooth) {
      cones = pullingSubdivision(rays, cones, strategy, seed, threads, verbose);
      if (!smoothRefinement(rays, cones, strategy, limit, threads, verbose))
        return nullptr;
    } else {
      cones = pullingSubdivision(rays, cones, strategy, seed, threads, verbose);
    }
    return encodeFan(raysInput, rays, cones);
  } catch (const libnormaliz::InterruptException&) {
    ERROR("rawSmoothFan interrupted");
    return nullptr;
  } catch (const std::exception& e) {
    ERROR(e.what());
    return nullptr;
  }
}

}

const Matrix* rawSimplicialFan(const Matrix* rays, M2_arrayint cones,
                               int strategy, int seed, int limit,
                               int threads, bool verbose)
{ return run(rays, cones, false, strategy, seed, limit, threads, verbose); }

const Matrix* rawSmoothFan(const Matrix* rays, M2_arrayint cones,
                           int strategy, int seed, int limit,
                           int threads, bool verbose)
{ return run(rays, cones, true, strategy, seed, limit, threads, verbose); }
