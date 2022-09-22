#include <gtest/gtest.h>

#include "MonomialOrder.hpp"

TEST(MonomialOrder, createFromArray)
{
  std::vector<std::pair<MO::order_type, std::vector<int>>> array1 {
    {MO::Lex, {4, 32}},
    {MO::Weights, {0, 1, 2, 3, 4}},
    {MO::GRevLex, {4, 32, 1, 1, 1, 1}},
    {MO::Position, {1}}
  };
  NewAgeMonomialOrder mo1 {8, array1};
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

