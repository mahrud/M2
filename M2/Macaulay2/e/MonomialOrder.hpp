#pragma once

#include <string>
#include <vector>
#include <iostream>

#include "interface/monomial-ordering.h" // Just for enum MonomialOrdering_type

/// TODO: this is in some sense a todo list for my current modifications
/// In Macaulay2, a MonomialOrder is a data structure that allows teh following options:
///   - Create a monomial order from its parts
///   - Encode a monomial into a list of ints for which comparison and multipliciation is easy.
///   - Decode a monomial back to an exponent vector, also a "varpower" product (sparse monomial).
///   - Comparison function
///   - Multiplication of two monomials, with overflow checking.
///   - to/from a "matrix" description of the monomial order.  This is lossy, in that
///     the matrix representation doesn't know about the packing information in the order.
/// Additionally, there are methods for querying aspects of the monomial order.
///   - How many slots (i.e. ints in the encoded monomial) correspond to the first d parts of the monomial order?)
///   - Is this a global ordering?
///   - Which variables are local?
///   - How many (and which) variables allow negative exponents?
///   - Is this graded reverse lex order?  Is it Lex order?
///   - What else is needed?
/// Questions remaining:
///   - where best to put torsion in monoids?  In the mult function?

/// MO = monomial order types, internal functions and types for monomial orders */
namespace MO {
  enum order_type {
    Lex = 1, /**< lexicographic order on a specified number of variables */
    GRevLex = 4, /**< graded reverse lexicographic order
                           on a specified number of variables,
                           with weights 1 for each variable */
    GRevLexWeights = 7, /**< graded reverse lexicographic order
                       with specified positive weights, on a number of variables */
    RevLex = 10, /**< reverse lexicographic order
                           on a specified number of variables.
                           Not a global order! */
    Weights = 11,
    GroupLex = 12, // previously MO_LAURENT
    GroupRevLex = 13, // previously MO_LAURENT_REVLEX
    PositionUp = 15,
    PositionDown = 16
  };
  //Previous names and values  
  // MO_LEX = 1, 
  // MO_GREVLEX = 4,
  // MO_GREVLEX_WTS = 7,
  // MO_REVLEX = 10,
  // MO_WEIGHTS = 11,
  // MO_LAURENT = 12,        /* Lex order here */
  // MO_LAURENT_REVLEX = 13, /* Rev lex order here */
  // MO_POSITION_UP = 15,
  // MO_POSITION_DOWN = 16

  class Block
  {
  private:
    enum order_type mType;
  
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
    Block(enum order_type t,
          const std::vector<int>& data,
          int& next_var,
          int& next_slot
          );

    enum order_type type() const { return mType; }
    int numVars() const { return mNumVars; }
    auto weights() const -> const std::vector<int>& { return mWeights; }
  };
}; // end namespace MO

///////////////////////////////////////////////////////////////////////

class NewAgeMonomialOrder
{
  friend std::ostream& operator<<(std::ostream&, const NewAgeMonomialOrder& mo);
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

  using MOInput = std::vector<std::pair<enum MO::order_type, std::vector<int>>>;
  NewAgeMonomialOrder(int nvars, MOInput input);

  MOInput toMonomialType() const;

  /// TODO
  bool isWellDefined(int verbose_level) const;

  /// TODO
  std::ostream debugOutput(std::ostream& o) const;
public:
  // informational functions
public:
  // encoding/decoding/comparison of monomials.
};

std::ostream& operator<<(std::ostream& o, const NewAgeMonomialOrder& mo);

/*
// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e "
// indent-tabs-mode: nil
// End:
*/
