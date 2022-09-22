#pragma once

#include <string>
#include <vector>

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

struct MonomialOrder_rec
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
};

typedef struct MonomialOrder_rec MonomialOrder;
MonomialOrder *monomialOrderMake(const MonomialOrdering *mo);
void monomialOrderFree(MonomialOrder *mo);
void monomialOrderEncodeFromActualExponents(const MonomialOrder *mo,
                                            const_exponents a,
                                            monomial b);
void monomialOrderDecodeToActualExponents(const MonomialOrder *mo,
                                          const_monomial a,
                                          exponents_t b);

std::vector<bool> laurentVariables(const MonomialOrder* mo);

namespace MO { /** MO = monomial order building blocks */
  enum order_type {
    Lex = 1,            /**< lexicographic order on a specified
                             number of variables */
    GRevLex = 4,        /**< graded reverse lexicographic order
                           on a specified number of variables,
                           with (positive) weights for the grading */
    RevLex = 10,        /**< reverse lexicographic order
                           on a specified number of variables.
                           Not a global order! */
    Weights = 11,
    GroupLex = 12,    /* Lex order here */
    GroupRevLex = 13, /* Rev lex order here */
    Position = 16       /**< 1 or -1 for Position => Up, Position => Down.
                           This determines the order on free modules */
  };

  class Block
  {
    order_type mType;
  
    /** number of bits per exponent/value for this type (currently 8, 16, or 32)
     *
     */
    int mPacking;

    /** number of variables for this block
     *  for Lex, GRevLex, RevLex, GroupLex, GroupRevLex: number of variables
     *  for GRevLexWeights: number of variables
     *  PositionUp, PositionDown: 0
     *  Weights: 0.
     */
    int mNumVars; 

    /// index of the first variable for this block
    int mStartVariable;

    /// the starting 
    int mStartEncoded;
    int mNumEncoded;

    std::vector<int> mWeights; // Only used for MO::Weights, MO::GRevLex

  public:
    /** constructor */
    Block(order_type t,
          const std::vector<int>& data,
          int& next_var,
          int& next_slot
          );
  };
};  // end namespace MO

class NewAgeMonomialOrder
{
private:
  // The sequence of blocks for the monomial order
  std::vector<MO::Block> mParts; 

  /// Number of variables
  int mNumVars;

  /// Size of encoded monomial
  int mEncodedSize;

  /// which block contains the component
  int mComponentLocation;

  /// if true, components are encoded via negation, so comparison will be faster
  bool mPositionUp;
  
  ///
  std::vector<int> mIsLocalVariable;
  std::vector<int> mIsGroupVariable;
  
public:
  using MOInput = std::vector<std::pair<MO::order_type, std::vector<int>>>;
  NewAgeMonomialOrder(int nvars, MOInput input);

  MOInput toMonomialType() const;

public:
  // informational functions
public:
  // encoding/decoding/comparison of monomials.
};

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
