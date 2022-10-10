#include <gtest/gtest.h>
#include <iostream>

#include "MonomialOrder.hpp"

TEST(MonomialOrder, createFromArray)
{
  std::vector<std::pair<MO::order_type, std::vector<int>>> array1 {
    {MO::Lex, {3}},
    {MO::Weights, {0, 1, 2, 3, 4}},
    {MO::GRevLex, {4}},
    {MO::GRevLexWeights, {1, 1, 1, 1, 2}},
    {MO::PositionUp, {}}
  };
  NewAgeMonomialOrder mo1 {12, array1};
  auto array2 = mo1.toMonomialType();
  EXPECT_EQ(array1, array2);
  std::cout << "--- display monomialType---" << std::endl;
  std::cout << mo1.toMonomialType() << std::endl;
  std::cout << "--- debug display ---" << std::endl;
  mo1.debugDisplay(std::cout);
  
  NewAgeMonomialOrder mo2 {12, {
    {MO::Lex, {3}},
    {MO::Weights, {0, 1, 2, 3, 4}},
    {MO::GRevLex, {4}},
    {MO::GRevLexWeights, {1, 1, 1, 1, 2}},
    {MO::PositionUp, {}}
    }};

  std::cout << "mo2 = " << mo2 << std::endl;
  mo2.debugDisplay(std::cout);


  NewAgeMonomialOrder mo3 {9, {
    {MO::Lex, {3}},
    {MO::RevLex, {2}},
    {MO::GroupRevLex, {4}},
    {MO::PositionUp, {}}
    }};

  std::cout << "mo3 = " << mo3 << std::endl;
  mo3.debugDisplay(std::cout);

}

TEST(MonomialOrder, packing)
{
  std::vector<std::pair<MO::order_type, std::vector<int>>> array1 {
    {MO::Packing, {4}},
    {MO::Lex, {3}},
    {MO::Weights, {0, 1, 2, 3, 4}},
    {MO::Packing, {2}},
    {MO::GRevLex, {4}},
    {MO::GRevLexWeights, {1, 1, 1, 1, 2}},
    {MO::PositionUp, {}}
  };
  NewAgeMonomialOrder mo1 {12, array1};
  auto array2 = mo1.toMonomialType();
  EXPECT_EQ(array1, array2);
  std::cout << mo1.toMonomialType() << std::endl;
  std::cout << "--- debug display ---" << std::endl;
  mo1.debugDisplay(std::cout);

  std::cout << "mo1 = " << mo1 << std::endl;
  
}

// TEST(MonomialOrder, arrayFromMonomialOrder)
// {
// }

// TEST(MonomialOrder, matrixFromMonomialOrder)
// {
// }

