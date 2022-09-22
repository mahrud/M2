#include "interface/monomial-ordering.h"

#include "MonomialOrder.hpp"

M2_arrayint rawMonomialOrderingToMatrix(const struct MonomialOrdering *mo)
{
  bool base;
  std::vector<int> mat;
  M2_arrayint result = nullptr;
  int component_is_before_row = 0;
  int component_direction = 0;
  if (monomialOrderingToMatrix(
          *mo, mat, base, component_direction, component_is_before_row))
    {
      int top = static_cast<int>(mat.size());
      result = M2_makearrayint(top + 3);
      for (int i = 0; i < top; i++) result->array[i] = mat[i];
      result->array[top] = (base ? 1 : 0);
      result->array[top + 1] = component_direction;
      result->array[top + 2] = component_is_before_row;
    }
  return result;
}

// Local Variables:
// indent-tabs-mode: nil
// End:
