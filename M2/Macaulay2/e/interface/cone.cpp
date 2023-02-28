// Written in 2021 by Mahrud Sayrafi

#include "interface/cone.h"

#include <M2/math-include.h>

#include "debug.hpp"
#include "interface/gmp-util.h"
#include "interface/matrix.h"
#include "matrix-con.hpp"
#include "matrix.hpp"
#include "relem.hpp"

#include <libnormaliz/cone.h>
#include <vector>

typedef mpz_class Integer;

/**
 * \ingroup cones
 */

const Matrix /* or null */ *rawFourierMotzkin(const Matrix *C)
{
  try
    {
      // TODO: generalize the input type, in particular to allow lineality space
      const Ring *R = C->get_ring();
      const size_t c = C->n_cols();  // rank of ambient lattice
      const size_t r = C->n_rows();  // number of cone inequalities

      auto ineqs = libnormaliz::Matrix<Integer>(r, c);
      for (size_t i = 0; i < r; i++)
        for (size_t j = 0; j < c; j++)
	  // libnormaliz uses A*x >= 0, Macaulay2 uses A*x <= 0
          ineqs[i][j] = (-1) * static_cast<Integer>(C->elem(i, j).get_mpz());

      auto cone = libnormaliz::Cone<Integer>(libnormaliz::Type::inequalities, ineqs);
      auto rays = cone.getExtremeRays();
      size_t n = rays.size();  // number of extremal rays

      MatrixConstructor mat(R->make_FreeModule(n), c);
      for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < c; j++)
          {
            mpz_ptr z = newitem(__mpz_struct);
            mpz_init_set(z, rays[i][j].get_mpz_t());
            mpz_reallocate_limbs(z);
            mat.set_entry(i, j, ring_elem(z));
          }

      return mat.to_matrix();
  } catch (const exc::engine_error &e)
    {
      ERROR(e.what());
      return nullptr;
  }
}

const Matrix /* or null */ *rawFourierMotzkinEqs(const Matrix *A, const Matrix *B)
{
  try
    {
      // TODO: generalize the input type, in particular to allow lineality space
      const Ring *R = A->get_ring();
      const size_t c = A->n_cols();  // rank of ambient lattice
      const size_t r = A->n_rows();  // number of cone inequalities

      auto ineqs = libnormaliz::Matrix<Integer>(r, c);
      for (size_t i = 0; i < r; i++)
        for (size_t j = 0; j < c; j++)
	  // libnormaliz uses A*x >= 0, Macaulay2 uses A*x <= 0
          ineqs[i][j] = (-1) * static_cast<Integer>(A->elem(i, j).get_mpz());

      const size_t c2 = B->n_cols();  // rank of ambient lattice
      const size_t r2 = B->n_rows();  // number of cone equations

      auto eqs = libnormaliz::Matrix<Integer>(r2, c2);
      for (size_t i = 0; i < r2; i++)
        for (size_t j = 0; j < c2; j++)
          eqs[i][j] = static_cast<Integer>(B->elem(i, j).get_mpz());

      auto cone = libnormaliz::Cone<Integer>(
        libnormaliz::Type::inequalities, ineqs,
        libnormaliz::Type::equations, eqs);
      auto rays = cone.getExtremeRays();
      auto span = cone.getMaximalSubspace();
      size_t m = rays.size();  // number of extremal rays
      size_t n = span.size();  // number of maximal subspaces

      // TODO: how can we return two separate matrices?
      // currently last row is wasted to return the number m
      MatrixConstructor mat(R->make_FreeModule(m + n + 1), c);
      mpz_ptr z = newitem(__mpz_struct);
      mpz_init_set_ui(z, (unsigned int) m);
      mat.set_entry(m + n, 0, ring_elem(z));
      for (size_t j = 0; j < c; j++)
	{
	  for (size_t i = 0; i < m; i++)
	    {
	      mpz_ptr z = newitem(__mpz_struct);
	      mpz_init_set(z, rays[i][j].get_mpz_t());
	      mpz_reallocate_limbs(z);
	      mat.set_entry(i, j, ring_elem(z));
	    }
	  for (size_t i = 0; i < n; i++)
	    {
	      mpz_ptr z = newitem(__mpz_struct);
	      mpz_init_set(z, span[i][j].get_mpz_t());
	      mpz_reallocate_limbs(z);
	      mat.set_entry(i + m, j, ring_elem(z)); // shifted by m
	    }
	}

      return mat.to_matrix();
  } catch (const exc::engine_error &e)
    {
      ERROR(e.what());
      return nullptr;
  }
}

const Matrix /* or null */ *rawHilbertBasis(const Matrix *C)
{
  try
    {
      // TODO: Check that C is over ZZ
      // TODO: for cones over QQ, lift to ZZ first
      // TODO: Normaliz also supports algebraic cones defined over
      // algebraic number fields embedded in RR
      const Ring *R = C->get_ring();
      const size_t c = C->n_cols();  // rank of ambient lattice
      const size_t r = C->n_rows();  // number of cone rays

      auto rays = libnormaliz::Matrix<Integer>(r, c);
      for (size_t i = 0; i < r; i++)
        for (size_t j = 0; j < c; j++)
          rays[i][j] = static_cast<Integer>(C->elem(i, j).get_mpz());

      auto cone = libnormaliz::Cone<Integer>(libnormaliz::Type::cone, rays);
      // cone.compute(libnormaliz::ConeProperty::HilbertBasis,
      //              libnormaliz::ConeProperty::DefaultMode);
      auto HB = cone.getHilbertBasis();
      size_t n = HB.size();  // number of basis elements

      MatrixConstructor mat(R->make_FreeModule(n), c);
      for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < c; j++)
          {
            mpz_ptr z = newitem(__mpz_struct);
            mpz_init_set(z, HB[i][j].get_mpz_t());
            mpz_reallocate_limbs(z);
            mat.set_entry(i, j, ring_elem(z));
          }

      return mat.to_matrix();
  } catch (const exc::engine_error &e)
    {
      ERROR(e.what());
      return nullptr;
  }
}

const Matrix /* or null */ *rawInteriorVector(const Matrix *C)
{
  try
    {
      // WIP: this doesn't work in all cases yet
      // TODO: generalize the input type
      const Ring *R = C->get_ring();
      const size_t c = C->n_cols();  // rank of ambient lattice
      const size_t r = C->n_rows();  // number of cone inequalities

      auto rays = libnormaliz::Matrix<Integer>(r, c);
      for (size_t i = 0; i < r; i++)
        for (size_t j = 0; j < c; j++)
          rays[i][j] = static_cast<Integer>(C->elem(i, j).get_mpz());

      auto cone = libnormaliz::Cone<Integer>(libnormaliz::Type::cone, rays);
      cone.compute(libnormaliz::ConeProperty::ExtremeRays);
      auto vect = cone.getGrading();
      // auto gcd = cone.getGradingDenom();

      MatrixConstructor mat(R->make_FreeModule(1), c);
      for (size_t j = 0; j < c; j++)
        {
          mpz_ptr z = newitem(__mpz_struct);
          mpz_init_set(z, vect[j].get_mpz_t());
          mpz_reallocate_limbs(z);
          mat.set_entry(0, j, ring_elem(z));
        }

      return mat.to_matrix();
  } catch (const exc::engine_error &e)
    {
      ERROR(e.what());
      return nullptr;
  }
}
