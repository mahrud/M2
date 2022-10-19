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

TEST(MonomialOrder, localVariables)
{
  NewAgeMonomialOrder mo1 {10, {
      {MO::Packing, {4}},
      {MO::Lex, {3}},
      {MO::RevLex, {2}},
      {MO::GroupRevLex, {2}},
      {MO::PositionUp, {}},
      {MO::Weights, {0, 0,0,0,0,0,0,0,1, -1, 1}},
      {MO::Lex, {3}}
    }};

  std::cout << mo1 << std::endl;

  auto locals = mo1.localVariablesPredicate();
  std::vector<bool> locals_ans {
    false, false, false,
    true, true,
    true, true,
    false, true, false
  };
  EXPECT_EQ(locals, locals_ans);

  auto inverts = mo1.invertibleVariablesPredicate();
  std::vector<bool> inverts_ans {
    false, false, false,
    false, false,
    true, true,
    false, false, false
  };
  EXPECT_EQ(inverts, inverts_ans);
}

TEST(MonomialOrder, matrixFromMonomialOrder)
{
  NewAgeMonomialOrder mo1 {6, {
      {MO::Packing, {4}},
      {MO::Lex, {3}},
      {MO::PositionUp, {}},
      {MO::Weights, {0, 0,0,0,1, -1, 1}},
      {MO::Lex, {3}}
    }};

  std::vector<int> mat;
  bool base_is_revlex;
  int component_dir;
  int comp_is_before_row;
  EXPECT_TRUE(mo1.monomialOrderingToMatrix(mat, base_is_revlex, component_dir, comp_is_before_row));
  std::cout << mat << std::endl;
  EXPECT_FALSE(base_is_revlex);
  EXPECT_EQ(component_dir, 1);
  EXPECT_EQ(comp_is_before_row, 3);
}

TEST(MonomialOrder, matrixFromMonomialOrder2)
{
  NewAgeMonomialOrder mo1 {6, {
      {MO::Packing, {4}},
      {MO::Lex, {3}},
      {MO::PositionUp, {}},
      {MO::GRevLex, {3}}
    }};

  std::vector<int> mat;
  bool base_is_revlex;
  int component_dir;
  int comp_is_before_row;
  EXPECT_TRUE(mo1.monomialOrderingToMatrix(mat, base_is_revlex, component_dir, comp_is_before_row));
  std::cout << mat << std::endl;
  EXPECT_TRUE(base_is_revlex);
  EXPECT_EQ(component_dir, 1);
  EXPECT_EQ(comp_is_before_row, 3);
}

