#include <iostream>
#include <sstream>

#include "MonomialOrder.hpp"
#include "exceptions.hpp"

namespace MO {
  MonomialBlock::MonomialBlock(order_type t,
                               const std::vector<int>& data,
                               int packing,
                               int& next_var
               )
    : mType(t),
      mPacking(packing),
      mNumVars(0),
      mStartVariable(next_var) // will change for MO::Weights.
  {
    if (t == MO::PositionUp or t == MO::PositionDown)
      {
        mPacking = 1;
        mStartVariable = 0;
        return;
      }
    
    if (data.size() == 0)
      {
        throw exc::engine_error("expected each monomial block to have positive number of variables");
      }

    if (t == MO::Weights)
      {
        mNumVars = 0;
        mStartVariable = data[0];
        std::copy(data.cbegin()+1, data.cend(), std::back_inserter(mWeights));
        return;
      }

    mNumVars = (data.size() >= 2 ? data.size() : data[0]);
    if (mNumVars <= 0)
      {
        throw exc::engine_error("expected each monomial block to have positive number of variables");
      }

    switch (t) {
    case MO::GroupLex: // Lex
    case MO::GroupRevLex:
      mPacking = 1;
      break;
    case MO::GRevLexWeights:
      //TODO: check that these are all positive.
      std::copy(data.cbegin(), data.cend(), std::back_inserter(mWeights));
      break;
    case MO::Lex:
    case MO::RevLex:
    case MO::GRevLex:
    // The following should not occur here.
    case MO::Weights:
    case MO::PositionUp:
    case MO::PositionDown:
      break;
    case MO::Packing:
      throw exc::engine_error("internal error: received incorrect monomial block");
    };

    next_var += mNumVars;
  }

  std::vector<int> prepend(int a, const std::vector<int>& b)
  {
    std::vector<int> result;
    result.push_back(a);
    for (auto c : b) result.push_back(c);
    return result;
  }
  
  std::pair<enum MO::order_type, std::vector<int>> MonomialBlock::toMonomialType() const
  {
    std::pair<enum MO::order_type, std::vector<int>> result;
    std::vector<int> newwts; // I don't know how to easily get around this?
    switch (mType)
      {
      case MO::Lex:            return {MO::Lex, {mNumVars}};
      case MO::RevLex:         return {MO::RevLex, {mNumVars}};
      case MO::GroupLex:       return {MO::GroupLex, {mNumVars}};
      case MO::GroupRevLex:    return {MO::GroupRevLex, {mNumVars}};
      case MO::GRevLex:        return {MO::GRevLex, {mNumVars}};
      case MO::GRevLexWeights: return {MO::GRevLexWeights, mWeights};
      case MO::Weights:        return {MO::Weights, prepend(mStartVariable, mWeights)};
      case MO::PositionUp:     return {MO::PositionUp, {}};
      case MO::PositionDown:   return {MO::PositionUp, {}};
      case MO::Packing:        return {MO::Packing, {100}}; // This should never occur.  TODO: error here
      }
  }

  std::string toString(MO::order_type t)
  {
    switch (t) {
    case MO::Lex:            return "Lex";
    case MO::RevLex:         return "RevLex";
    case MO::GroupLex:       return "GroupLex";
    case MO::GroupRevLex:    return "GroupRevLex";
    case MO::GRevLex:        return "GRevLex";
    case MO::GRevLexWeights: return "GRevLexWeights";
    case MO::Weights:        return "Weights";
    case MO::PositionUp:     return "PositionUp";
    case MO::PositionDown:   return "PositionDown";
    case MO::Packing:        return "Packing";
    }
  }
  
  void MonomialBlock::debugDisplay(std::ostream& o) const
  {
    o << "  " << mType;
    o << "[#vars=" << mNumVars << " start=" << mStartVariable;
    if (mWeights.size() > 0)
      o << "  weights=" << mWeights;
    //    if (mPacking != 1)
      o << " packing " << mPacking;
    o << "]";
  }

  EncodingBlock::EncodingBlock(MonomialBlock b,
                  int& next_slot
                  )
    : mBlock(b),
      mStartEncoded(next_slot)
  {
    switch (b.type()) {
    case MO::Lex:
    case MO::RevLex:
      mNumEncoded = b.numVars() / b.packing();
      break;
    case MO::GroupLex: // Lex
    case MO::GroupRevLex:
      mNumEncoded = b.numVars();
      break;
    case MO::GRevLex:
      mNumEncoded = 1 + b.numVars() / b.packing();
      break;
    case MO::GRevLexWeights:
      mNumEncoded = 1 + b.numVars() / b.packing();
      break;
    case MO::Weights:
      mNumEncoded = 1;
      break;
    case MO::PositionUp:
    case MO::PositionDown:
    case MO::Packing:
      // this is handled higher up the food chain
      break;
    };
    next_slot += mNumEncoded;
  }

  void EncodingBlock::debugDisplay(std::ostream& o) const
  {
    mBlock.debugDisplay(o);
    o << "  encoded [start,#] = [" << mStartEncoded << ", " << mNumEncoded << "]";
  }
  
}; // end of MO namespace

std::ostream& operator<<(std::ostream& o, MO::order_type t)
{
  return (o << MO::toString(t));
}

std::ostream& operator<<(std::ostream& o, MO::MonomialBlock b)
{
  o << b.type() << " => ";
  switch (b.type()) {
  case MO::Lex:
  case MO::RevLex:
  case MO::GroupLex:
  case MO::GroupRevLex:
  case MO::GRevLex:
    o << b.numVars();
    break;
  case MO::GRevLexWeights:
  case MO::Weights:
    o << b.weights();
    break;
  case MO::PositionUp:
  case MO::PositionDown:
    o << "{}";
    break;
  case MO::Packing:
    o << b.packing();
    break;
  }
  return o;    
  //  return (o << MO::toString(t));
}

std::ostream& operator<<(std::ostream& o, const std::pair<enum MO::order_type, std::vector<int>>& a)
{
  o << "{" << a.first; //<< " => " << a.second << "}";
  return o;
}

std::ostream& operator<<(std::ostream& o, const NewAgeMonomialOrder::MOInput& a)
  {
    for (const auto& a1 : a)
      {
        o << "{" << a1.first << ", " << a1.second << "}" << std::endl;
      }
    return o;
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

NewAgeMonomialOrder::NewAgeMonomialOrder(int numvars, const MOInput& input)
  : mNumVars(numvars)
{
  int next_var = 0;
  int packing = 1;

  // Create each block in the monomial order in turn.
  for (const auto &i : input)
    {
      // We do not store directly the Packing elements
      if (i.first == MO::Packing)
        {
          if (i.second.size() != 1)
            throw exc::engine_error("packing size not given correctly");
          packing = i.second[0];
          if (not (packing == 1 or packing == 2 or packing == 4))
            throw exc::engine_error("packing size must be 1, 2 or 4");
        }
      else
        mParts.push_back(MO::MonomialBlock(i.first, i.second, packing, next_var));
    }

  if (numvars != next_var)
    throw exc::engine_error("wrong number of variables");
}

auto NewAgeMonomialOrder::toMonomialType() const -> MOInput
{
  MOInput result;
  int packing = 1;
  for (const auto &m : mParts)
    {
      if (m.packing() != packing
          and (m.type() == MO::Lex or m.type() == MO::GRevLex or m.type() == MO::GRevLexWeights))
        {
          result.push_back({MO::Packing, {m.packing()}});
          packing = m.packing();
        }
      result.push_back(m.toMonomialType());
    }
  return result;
}

std::ostream& operator<<(std::ostream& o, const NewAgeMonomialOrder& mo) {
  int packing = 1;
  o << "MonomialOrder => {";
  for (int i = 0; i < mo.mParts.size(); ++i)
    {
      const auto& b = mo.mParts[i];
      if (i == 0)
        o << "\n    ";
      else
        o << ",\n    ";
      if (b.packing() != packing
          and (b.type() == MO::Lex or b.type() == MO::GRevLex or b.type() == MO::GRevLexWeights))
        {
          o << "Packing => " << b.packing() << ",\n    ";
        }
      o << b;
    }
  o << "\n    }";    
  return o;
}

void NewAgeMonomialOrder::debugDisplay(std::ostream& o) const
{
  o << "MonomialOrder:" << std::endl;
  o << "  #vars    = " << mNumVars << std::endl;
  for (const auto& m : mParts)
    {
      m.debugDisplay(o);
      o << std::endl;
    }
}

//////////////////////////// Encoder code /////////////////////////////////////
MonomialEncoder::MonomialEncoder(const NewAgeMonomialOrder& mo)
{
  int next_var = 0;
  int next_slot = 0;
  int packing = 1;

  mComponentOccursBeforeThisBlock = -1;
#if 0
  // TODO: reinstate the code in this block.
  // Create each block in the monomial order in turn.
  for (const auto &i : input)
    {
      if (i.first == MO::PositionUp)
        {
          mComponentOccursBeforeThisBlock = mParts.size();
          mEncodedComponentLocation = next_slot;
          mPositionUp = true;
          next_slot++;
          continue;
        }
      else if (i.first == MO::PositionDown)
        {
          mComponentOccursBeforeThisBlock = mParts.size();
          mEncodedComponentLocation = next_slot;
          mPositionUp = false;
          next_slot++;
          continue;
        }
      mParts.push_back(MO::Block(i.first, i.second, packing, next_var, next_slot));
    }

  if (mComponentOccursBeforeThisBlock == -1)
    {
      // In this case, Position wasn't given.  We add in Position => Up here.
      mComponentOccursBeforeThisBlock = mParts.size();
      mEncodedComponentLocation = next_slot;
      mPositionUp = true;
      next_slot++;
    }
  if (numvars != next_var)
    throw exc::engine_error("wrong number of variables");
#endif  
  mEncodedSize = next_slot;
}

void MonomialEncoder::debugDisplay(std::ostream& o) const
{
  o << "MonomialEncoder:" << std::endl;
  o << "  #vars    = " << mNumVars << std::endl;
  o << "  #encoded = " << mEncodedSize << std::endl;
  o << "  #compslot = " << mEncodedComponentLocation << std::endl;
  o << "  #compblock = " << mComponentOccursBeforeThisBlock << std::endl;
  o << "  #comp up = " << (mPositionUp ? "true" : "false") << std::endl;
  for (const auto& m : mParts)
    {
      m.debugDisplay(o);
      o << std::endl;
    }
}

// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e "
// indent-tabs-mode: nil
// End:

