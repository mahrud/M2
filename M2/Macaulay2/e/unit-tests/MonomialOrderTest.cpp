#include <gtest/gtest.h>
#include <iostream>

#include "MonomialOrder.hpp"

TEST(MonomialOrder, createFromArray)
{
  std::vector<std::pair<MO::order_type, std::vector<int>>> array1 {
    {MO::Lex, {4, 32}},
    {MO::Weights, {0, 1, 2, 3, 4}},
    {MO::GRevLexWeights, {4, 32, 1, 1, 1, 1}},
    {MO::PositionUp, {}}
  };
  NewAgeMonomialOrder mo1 {8, array1};

  NewAgeMonomialOrder mo2 {8, {
    {MO::Lex, {4, 32}},
    {MO::Weights, {0, 1, 2, 3, 4}},
    {MO::GRevLexWeights, {4, 32, 1, 1, 1, 1}},
    {MO::PositionUp, {}}
    }};

  std::cout << "mo2 = " << mo2 << std::endl;
}

// TEST(MonomialOrder, createFromMatrix)
// {
// }

// TEST(MonomialOrder, arrayFromMonomialOrder)
// {
// }

// TEST(MonomialOrder, matrixFromMonomialOrder)
// {
// }

