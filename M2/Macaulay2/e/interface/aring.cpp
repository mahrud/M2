#include "interface/aring.h"

#include <algorithm>
#include <utility>
#include <vector>
#include <memory>

#include "aring-gf-flint-big.hpp"
#include "aring-gf-flint.hpp"
#include "aring-gf-m2.hpp"
#include "aring-glue.hpp"
#include "aring-qq.hpp"
#include "aring-tower.hpp"
#include "aring-zz-flint.hpp"
#include "aring-zzp-flint.hpp"
#include "aring-zzp.hpp"
#include "exceptions.hpp"
#include "polyring.hpp"
#include "relem.hpp"

const RingQQ *globalQQ;

void initializeRationalRing()
{
  globalQQ = RingQQ::create();
}

const RingQQ *rawARingQQ() { return globalQQ; }
const Ring * /* or null */ rawARingZZFlint()
{
  return M2::ConcreteRing<M2::ARingZZ>::create();
}

const Ring * /* or null */ rawARingQQFlint()
{
  return M2::ConcreteRing<M2::ARingQQFlint>::create();
}

const Ring /* or null */ *rawARingZZp(unsigned long p)
{
  if (p <= 1 || p >= 32750)
    {
      ERROR("ZZP: expected a prime number p in range 2 <= p <= 32749");
      return 0;
    }
  return M2::ConcreteRing<M2::ARingZZp>::create(p);
}

// TODO: p can be an fmpz
const Ring /* or null */ *rawARingZZpFlint(unsigned long p)
{
  return M2::ConcreteRing<M2::ARingZZpFlint>::create(p);
}

static const PolynomialRing * /* or null */ checkGaloisFieldInput(
    const RingElement *f)
{
  // Check that the ring R of f is a polynomial ring in one var over a ZZ/p
  // If any of these fail, then return 0.
  const PolynomialRing *R = f->get_ring()->cast_to_PolynomialRing();
  if (R == 0)
    {
      ERROR("expected poly ring of the form ZZ/p[x]/(f)");
      return 0;
    }
  if (R->n_vars() != 1)
    {
      ERROR("expected poly ring of the form ZZ/p[x]/(f)");
      return 0;
    }
  if (R->n_quotients() != 1)
    {
      ERROR("expected poly ring of the form ZZ/p[x]/(f)");
      return 0;
    }
  if (R->characteristic() == 0)
    {
      ERROR("expected poly ring of the form ZZ/p[x]/(f)");
      return 0;
    }
  return R;
}

const Ring /* or null */ *rawARingGaloisFieldM2(const RingElement *f)
{
  const PolynomialRing *R = checkGaloisFieldInput(f);
  if (R == 0) return 0;  // error message has already been logged
  try
    {
      return M2::ConcreteRing<M2::ARingGFM2>::create(*R, f->get_value());
  } catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return NULL;
  }
}

const Ring /* or null */ *rawARingGaloisFieldFlintBig(const RingElement *f)
{
  const PolynomialRing *R = checkGaloisFieldInput(f);
  if (R == 0) return 0;  // error message has already been logged
  try
    {
      return M2::ConcreteRing<M2::ARingGFFlintBig>::create(*R, f->get_value());
  } catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return NULL;
  }
}

const Ring /* or null */ *rawARingGaloisFieldFlintZech(const RingElement *f)
{
  const PolynomialRing *R = checkGaloisFieldInput(f);
  if (R == 0) return 0;  // error message has already been logged
  try
    {
      return M2::ConcreteRing<M2::ARingGFFlint>::create(*R, f->get_value());
  } catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return NULL;
  }
}

// TODO: prime can be an fmpz
const Ring /* or null */ *rawARingGaloisField(long prime, long dimension)
{
  if (prime < 2)
    {
      ERROR("GF: expected a prime number p");
      return 0;
    }
  if (dimension < 0)
    {
      ERROR("GF: expected non-negative extension degree");
      return 0;
    }
  try
    {
      return M2::ConcreteRing<M2::ARingGFFlint>::create(prime, dimension);
  } catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return NULL;
  }
}

M2_arrayintOrNull rawARingGFPolynomial(const Ring *R)
{
  auto *RGF = dynamic_cast<const M2::ConcreteRing<M2::ARingGFFlint> *>(R);
  if (RGF == 0)
    {
      ERROR("expected a GaloisField");
      return 0;
    }
  const M2::ARingGFFlint &A = RGF->ring();
  return A.getModulusPolynomialCoeffs();
}

M2_arrayintOrNull rawARingGFCoefficients(const RingElement *f)
{
  auto *RGF = dynamic_cast<const M2::ConcreteRing<M2::ARingGFFlint> *>(f->get_ring());
  if (RGF == 0)
    {
      ERROR("expected a GaloisField");
      return 0;
    }
  const M2::ARingGFFlint &A = RGF->ring();
  M2::ARingGFFlint::ElementType a;
  A.from_ring_elem(a, f->get_value());
  std::vector<long> poly;
  A.getSmallIntegerCoefficients(a, poly);
  return stdvector_to_M2_arrayint(poly);
}

const Ring /* or null */ *rawARingTower1(const Ring *K, M2_ArrayString names)
{
  try
    {
      const M2::ConcreteRing<M2::ARingZZpFlint> *Kp =
          dynamic_cast<const M2::ConcreteRing<M2::ARingZZpFlint> *>(K);
      if (Kp == 0)
        {
          ERROR("expected a base ring ZZ/p");
          return NULL;
        }
      const M2::ARingZZpFlint &A = Kp->ring();

      // Get the names into the correct form:
      auto varnames = M2_ArrayString_to_stdvector(names);
      M2::ARingTower *T = M2::ARingTower::create(A, varnames);
      return M2::ConcreteRing<M2::ARingTower>::create(
          std::unique_ptr<M2::ARingTower>(T));
  } catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return NULL;
  }
}

const Ring /* or null */ *rawARingTower2(const Ring *R1,
                                         M2_ArrayString new_names)
{
  try
    {
      const M2::ConcreteRing<M2::ARingTower> *K =
          dynamic_cast<const M2::ConcreteRing<M2::ARingTower> *>(R1);
      if (K == 0)
        {
          ERROR("expected a tower ring");
          return NULL;
        }
      const M2::ARingTower &A = K->ring();

      auto new_varnames = M2_ArrayString_to_stdvector(new_names);
      M2::ARingTower *T = M2::ARingTower::create(A, new_varnames);
      return M2::ConcreteRing<M2::ARingTower>::create(
          std::unique_ptr<M2::ARingTower>(T));
  } catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return NULL;
  }
}

const Ring /* or null */ *rawARingTower3(const Ring *R1,
                                         engine_RawRingElementArray eqns)
{
  try
    {
      const M2::ConcreteRing<M2::ARingTower> *K =
          dynamic_cast<const M2::ConcreteRing<M2::ARingTower> *>(R1);
      if (K == 0)
        {
          ERROR("expected a tower ring");
          return NULL;
        }
      const M2::ARingTower &A = K->ring();

      std::vector<M2::ARingTower::ElementType> extensions;

      for (int i = 0; i < eqns->len; i++)
        {
          const RingElement *f = eqns->array[i];
          M2::ARingTower::ElementType f1;
          if (f->get_ring() != R1)
            {
              ERROR("extension element has incorrect base ring");
              return 0;
            }
          A.from_ring_elem(f1, f->get_value());
          extensions.push_back(f1);
        }
      M2::ARingTower *T = M2::ARingTower::create(A, extensions);
      return M2::ConcreteRing<M2::ARingTower>::create(
          std::unique_ptr<M2::ARingTower>(T));
  } catch (const exc::engine_error& e)
    {
      ERROR(e.what());
      return 0;
  }
}

// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e "
// indent-tabs-mode: nil
// End:
