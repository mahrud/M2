#pragma once

#include <string>
#include <vector>
#include <iostream>

#include "interface/monomial-ordering.h"

#include "ExponentVector.hpp"

typedef struct mon_part_rec_
{
  enum MonomialOrdering_type type;
  int nvars;
  int *wts;
} * mon_part;

struct MonomialOrdering
{
  unsigned int _hash;
  unsigned int len;
  mon_part array[1];
};

class MonomialOrderings {
public:
  static std::string toString(const MonomialOrdering *mo);
  static MonomialOrdering* join(const std::vector<MonomialOrdering*>& M);
  static MonomialOrdering* product(const std::vector<MonomialOrdering*>& M);
  static MonomialOrdering* Lex(int nvars);
  static MonomialOrdering* Lex2(int nvars);
  static MonomialOrdering* Lex4(int nvars);
  static MonomialOrdering* GRevLex(int nvars);
  static MonomialOrdering* GRevLex2(int nvars);
  static MonomialOrdering* GRevLex4(int nvars);
  static MonomialOrdering* GRevLex(const std::vector<int>& wts);
  static MonomialOrdering* GRevLex2(const std::vector<int>& wts);
  static MonomialOrdering* GRevLex4(const std::vector<int>& wts);
  static MonomialOrdering* RevLex(int nvars);
  static MonomialOrdering* Weights(const std::vector<int>& wts);
  static MonomialOrdering* GroupLex(int nvars);
  static MonomialOrdering* GroupRevLex(int nvars);
  static MonomialOrdering* PositionUp();
  static MonomialOrdering* PositionDown();

  static MonomialOrdering* GRevLex(const std::vector<int>& wts, int packing);
};

typedef int *monomial;

typedef int32_t
    deg_t;  // this is the integer type to use for degrees and weights

typedef const int *const_monomial;

struct mo_block
{
  enum MonomialOrdering_type typ;
  int nvars;
  int nslots;
  int first_exp;
  int first_slot;
  int nweights;
  deg_t *weights;
};

struct MonomialOrder
{
  int nvars;
  int nslots;
  int nblocks;
  int nblocks_before_component;
  int nslots_before_component;
  int component_up; /* bool */
  deg_t
      *degs; /* 0..nvars: heuristic degree of each variable.  degs[nvars] = 1.
                 Assumption: degs[i] >= 1, for all i, and should be an integer.
                 Any graded rev lex block assumes graded wrt these degrees. */
  struct mo_block *blocks; /* 0..nblocks-1 with each entry a struct mo_block */
  int *is_laurent; /* 0..nvars-1: 0 or 1: 1 means negative exponents allowed */

public:
  enum overflow_type { OVER, OVER1, OVER2, OVER4 } * overflow;
  std::vector<int> overflow_flags() const;
};

MonomialOrder *monomialOrderMake(const MonomialOrdering *mo);
void monomialOrderFree(MonomialOrder *mo);
void monomialOrderEncodeFromActualExponents(const MonomialOrder *mo,
                                            const_exponents a,
                                            monomial b);
void monomialOrderDecodeToActualExponents(const MonomialOrder *mo,
                                          const_monomial a,
                                          exponents_t b);

std::vector<bool> laurentVariables(const MonomialOrder* mo);

bool monomialOrderingToMatrix(
    const struct MonomialOrdering& mo,
    std::vector<int>& mat,
    bool& base_is_revlex,
    int& component_direction,     // -1 is Down, +1 is Up, 0 is not present
    int& component_is_before_row  // -1 means: at the end. 0 means before the
                                  // order.
    // and r means considered before row 'r' of the matrix.
    );

/*
// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e "
// indent-tabs-mode: nil
// End:
*/
