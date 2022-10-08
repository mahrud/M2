#include <iostream>
#include <sstream>

#include "MonomialOrder.hpp"
#include "exceptions.hpp"

namespace MO {
  /// creates one block of the monomial order.
  /// data contains the information for the type
  ///  Lex, RevLex, GRevLex, 
  ///    data[0]: #vars in block
  ///    data[1]: packing factor: 1, 2, 4.
  ///    notes: GRevLex has a default grading vector: 1,1,...,1.
  ///  GroupLex, GroupRevLex (packing factor is 1: only one number per int/word).
  ///    data[0]: #vars in block
  ///    notes: no packing factor is allowed. (?)
  ///  Position
  ///    data[0] = 1 or -1.  1 is UP, -1 is DOWN.
  ///  Weights
  ///    data[0] = #elements in weight vector
  ///    data[1] = begin loc for weight vector
  ///    data[2 .. 2 - 1 + data[0]]: weight vector itself.

  ///  data[0]: packing size (8, 16, 32, or should it be 1,2,4?)
  ///  data[1]: #variables (meaning depends on the type).
  ///  GRevLex: data[1] = #vars in this block.
  ///  GRevLexWeights: data[1] = #vars in this block
  ///   data[2+0..2+(nvars-1)] the positive weight vector
  ///  Lex: data[1] = #vars in block
  ///  RevLex: data[1] = #vars in block
  ///  Weights: data[1] = start variable
  ///           data[2] = #weights
  ///           data[3]..data[3-1+#weights]: weight vector,
  ///              possible with negative entries.
  ///  Position: data[1] = 1 "up"
  //               or:    -1 "down"
  ///  GroupLex, GroupRevLex: data[1] = #vars in block
  Block::Block(order_type t,
               const std::vector<int>& data,
               int& next_var,
               int& next_slot
               )
    : mType(t)
      // mPacking(data[1]),
      // mNumVars(data[0]),
      // mStartVariable(next_var),
      // mStartEncoded(next_slot)
  {
    switch (t) {
    case MO::Lex:
    case MO::RevLex:
    case MO::GRevLex:
    case MO::GroupLex: // Lex
    case MO::GroupRevLex:
      mNumVars = data[0];
      mPacking = data[1];
      break;
    case MO::GRevLexWeights:
      mNumVars = data[0];
      mPacking = data[1];
      std::copy(data.cbegin()+2, data.cend(), std::back_inserter(mWeights));
      //TODO: check that these are all positive.
      break;
    case MO::Weights:
      mStartVariable = data[0];
      std::copy(data.cbegin()+1, data.cend(), std::back_inserter(mWeights));
      break;
    case MO::PositionUp:
    case MO::PositionDown:
      // this is handled higher up the food chain
      break;
    };
    next_var += mNumVars;
  }
}; // end of MO namespace

NewAgeMonomialOrder::NewAgeMonomialOrder(int numvars, MOInput input)
  : mNumVars(numvars)
{
  int next_var = 0;
  int next_slot = 0;
  
  // Create each block in the monomial order in turn.
  for (const auto &i : input)
    {
      if (i.first == MO::PositionUp)
        {
          mComponentLocation = next_slot;
          mPositionUp = true;
        }
      else if (i.first == MO::PositionDown)
        {
          mComponentLocation = next_slot;
          mPositionUp = false;
        }
      mParts.push_back(MO::Block(i.first, i.second, next_var, next_slot));
    }

  mEncodedSize = next_slot;
}

std::ostream& operator<<(std::ostream& o, const std::vector<int>& a)
{
  o << "{";
  for (int i = 0; i < a.size(); i++)
    {
      if (i > 0) o << ",";
      o << a[i];
    }
  o << "}";
  return o;
}

  
std::ostream& operator<<(std::ostream& o, const NewAgeMonomialOrder& mo) {
  o << "MonomialOrder => {";
  for (int i = 0; i < mo.mParts.size(); ++i)
    {
      const auto& p = mo.mParts[i];
      if (i == 0)
        o << "\n    ";
      else
        o << ",\n    ";
      switch (p.type())
        {
        case MO::Lex:
          o << "Lex => " << p.numVars();
          break;
        case MO::GRevLex:
          o << "GRevLex => " << p.numVars();
          break;
        case MO::GRevLexWeights:
          o << "GRevLex => " << p.weights();
          break;
        case MO::RevLex:
          o << "RevLex => " << p.numVars();
          break;
        case MO::Weights:
          o << "Weights => ";
          //TODO print weights
          //o << p.mWeights;
          break;
        case MO::GroupLex:
          o << "GroupLex => " << p.numVars();
          break;
        case MO::GroupRevLex:
          o << "GroupRevLex => " << p.numVars();
          break;
        case MO::PositionUp:
          o << "Position => Up";
          break;
        case MO::PositionDown:
          o << "Position => Down";
          break;
        default:
          o << "UNKNOWN";
          break;
        }
    }
  o << "\n    }";
  return o;
}

// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e "
// indent-tabs-mode: nil
// End:

