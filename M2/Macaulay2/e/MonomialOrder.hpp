#pragma once

#include <string>
#include <vector>
#include <iostream>

using ExponentVector = int*;
using EncodedMonomial = int*;

///TODO: maybe we need 2 classes here.
/// 1. MonomialOrder
/// 2. MonomialEncoder (and decoder).
///   this is like the current intermal monomial ordering code.
///   it has a list of parts (monomial order parts, and (firstslot, numslot, packing, overflow) info for each.
///   operations
///   - encodeExponentVector
///   - decodeToExponentVector
///   - mult (with overflow checks) (this needs the encoder for the over flow fields)
///   - comparison (this might not need the Encoder at all).

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
    Packing, /**< # exponents per (currently 32 bit) word.  Only affects Lex, GRevLex, GRevLexWeights
                     (but not the weight part).  Does not affect the monomial order. */
    Lex, /**< lexicographic order on a specified number of variables */
    GRevLex, /**< graded reverse lexicographic order
                           on a specified number of variables,
                           with weights 1 for each variable */
    GRevLexWeights, /**< graded reverse lexicographic order
                       with specified positive weights, on a number of variables */
    RevLex, /**< reverse lexicographic order
                           on a specified number of variables.
                           Not a global order! */
    Weights,
    GroupLex, // previously MO_LAURENT
    GroupRevLex, // previously MO_LAURENT_REVLEX
    PositionUp,
    PositionDown
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

  std::string toString(MO::order_type t);

  class MonomialBlock
  {
  private:
    enum order_type mType;
  
    /** This is the number of exponents to pack in to a word.  This is
        the one piece of data that isn't part of the monomial order.
        But we keep it, because it is easier to collect the
        information once from the front end
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
    std::vector<int> mWeights; // Only used for MO::Weights, MO::GRevLex
  public:
    /** constructor */
    MonomialBlock(enum order_type t,
          const std::vector<int>& data,
          int packing,
          int& next_var
          );

    enum order_type type() const { return mType; }
    int numVars() const { return mNumVars; }
    int firstVar() const { return mStartVariable; }
    auto weights() const -> const std::vector<int>& { return mWeights; }
    int packing() const { return mPacking; }
    void debugDisplay(std::ostream& o) const;
    std::pair<enum MO::order_type, std::vector<int>> toMonomialType() const;
  };

  class EncodingBlock
  {
  private:
    MonomialBlock mBlock;

    /// encoded locations (start, length).
    int mStartEncoded;
    int mNumEncoded;
  public:
    /** constructor */
    EncodingBlock(MonomialBlock b,
                  int& next_slot
                  );
    
    //TODO: reinstate the ones needed here.
    // enum order_type type() const { return mType; }
    // int numVars() const { return mNumVars; }
    // auto weights() const -> const std::vector<int>& { return mWeights; }
    void debugDisplay(std::ostream& o) const;
  };
}; // end namespace MO


///////////////////////////////////////////////////////////////////////

class NewAgeMonomialOrder
{
  friend std::ostream& operator<<(std::ostream&, const NewAgeMonomialOrder& mo);
private:
  /// Number of variables
  int mNumVars;

  // The sequence of blocks for the monomial order
  std::vector<MO::MonomialBlock> mParts; 

  // Informational about the order
  std::vector<int> mIsLocalVariable;
  std::vector<int> mIsGroupVariable;
  
public:
  using MOInput = std::vector<std::pair<enum MO::order_type, std::vector<int>>>;
  NewAgeMonomialOrder(int nvars, const MOInput& input);

  MOInput toMonomialType() const;

  /// TODO
  bool isWellDefined(int verbose_level) const;

  void debugDisplay(std::ostream& o) const;

  int numVars() const { return mNumVars; }
public:
  ////// informational functions ///////

  // Is the monomial order Lex or GRevLex?
  
  
  // The list of indices of variables which are < 1.
  std::vector<int> localVariables() const;
  std::vector<int> invertibleVariables() const;

  std::vector<bool> invertibleVariablesPredicate() const; // aka laurentVariables
  std::vector<bool> localVariablesPredicate() const; // aka non term order variables.

  bool hasLocalVariables() const { return localVariables().size() > 0; }
  bool hasInvertibleVariables() const { return invertibleVariables().size() > 0; }
public:
  //TODO  unencoded comparison.
  int compareExponentVectors(const ExponentVector &a, const ExponentVector* b) const;

  //TODO: is an unencoded monomial in the subring where the first n parts are 0?

  //TODO: to and from std::vector<int> data for mgb and the front end.

  bool monomialOrderingToMatrix(
                           std::vector<int>& mat,
                           bool& base_is_revlex,
                           int& component_direction,     // -1 is Down, +1 is Up, 0 is not present
                           int& component_is_before_row  // -1 means: at the end. 0 means before the
                           // order.
                           // and r means considered before row 'r' of the matrix.
                           ) const;
};

std::ostream& operator<<(std::ostream& o, const std::vector<int>& a);
std::ostream& operator<<(std::ostream& o, MO::order_type t);
std::ostream& operator<<(std::ostream& o, const std::pair<enum MO::order_type, std::vector<int>>& a);
std::ostream& operator<<(std::ostream& o, const NewAgeMonomialOrder& mo);
std::ostream& operator<<(std::ostream& o, const NewAgeMonomialOrder::MOInput& a);


/////////////// Encoder code /////////////////////////////
class MonomialEncoder
{
private:
  std::vector<MO::EncodingBlock> mParts; 

  /// Number of variables
  int mNumVars;

  /// Size of encoded monomial (includes component, if present).
  int mEncodedSize;

  /// the component comes after the following block #.
  int mComponentOccursBeforeThisBlock;

  /// the component value in the encoded monomial is at this location
  int mEncodedComponentLocation;
  
  /// if true, components are encoded via negation, so comparison will be faster
  bool mPositionUp;
  
public:
  MonomialEncoder(const NewAgeMonomialOrder& mo);

  //TODO.
  int compareEncodedMonomials();
  void multEncodedMonomials(); // can overflow
  void encode(const ExponentVector& a, EncodedMonomial& result) const;
  void decode(const EncodedMonomial& a, ExponentVector& b) const;

  void debugDisplay(std::ostream& o) const;
};


/*
// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e "
// indent-tabs-mode: nil
// End:
*/
