#include "MonomialOrderOld.hpp"
#include "exceptions.hpp"
#include "error.h"
#include <iostream>
#include <sstream>

static struct mon_part_rec_ *mo_make(enum MonomialOrdering_type type,
                                     int nvars,
                                     const int *wts)
{
  mon_part result;
  result = getmemstructtype(mon_part);
  result->type = type;
  result->nvars = nvars;
  if (wts != nullptr)
    {
      int i;
      result->wts = getmematomicvectortype(int, nvars);
      for (i = 0; i < nvars; i++) result->wts[i] = wts[i];
    }
  else
    result->wts = nullptr;
  return result;
}

static MonomialOrdering *make_mon_order(int n)
{
  static unsigned int next_hash = 23023421;
  MonomialOrdering *z = getmemarraytype(MonomialOrdering *, n);
  z->len = n;
  z->_hash = next_hash++;
  int i;
  for (i = 0; i < n; i++) z->array[i] = nullptr;
  return z;
}

static MonomialOrdering *M2_mo_offset(const MonomialOrdering *mo, int offset)
{
  int i, j;
  MonomialOrdering *result = make_mon_order(mo->len);
  for (i = 0; i < mo->len; i++)
    {
      mon_part p = mo->array[i];
      if (p->type != MO_WEIGHTS)
        result->array[i] = mo_make(p->type, p->nvars, p->wts);
      else
        {
          mon_part q = mo_make(MO_WEIGHTS, offset + p->nvars, nullptr);
          q->wts = getmemvectortype(int, q->nvars);
          for (j = 0; j < offset; j++) q->wts[j] = 0;
          for (; j < q->nvars; j++) q->wts[j] = p->wts[j - offset];
        }
    }
  return result;
}

static bool is_good(mon_part p)
{
  switch (p->type)
    {
      case MO_LEX:
      case MO_LEX2:
      case MO_LEX4:
      case MO_GREVLEX:
      case MO_GREVLEX2:
      case MO_GREVLEX4:
      case MO_GREVLEX_WTS:
      case MO_GREVLEX2_WTS:
      case MO_GREVLEX4_WTS:
      case MO_LAURENT:
      case MO_NC_LEX:
      case MO_LAURENT_REVLEX:
      case MO_REVLEX:
      case MO_WEIGHTS:
        return (p->nvars > 0);
      case MO_POSITION_UP:
      case MO_POSITION_DOWN:
        return true;
    }
  return false;
}


std::ostringstream &toString(std::ostringstream &o, int len, int *p)
{
  o << "{";
  for (int i = 0; i < len; i++)
    {
      if (i > 0) o << ",";
      o << p[i];
    }
  o << "}";
  return o;
}

std::ostringstream &ones(std::ostringstream &o, int len)
{
  o << "{";
  for (int i = 0; i < len; i++)
    {
      if (i > 0) o << ",";
      o << 1;
    }
  o << "}";
  return o;
}

std::string MonomialOrderings::toString(const MonomialOrdering *mo)
{
  std::ostringstream o;
  o << "MonomialOrder => {";
  for (int i = 0; i < mo->len; i++)
    {
      mon_part p = mo->array[i];
      bool p_ones = false;
      if (i == 0)
        o << "\n    ";
      else
        o << ",\n    ";
      switch (p->type)
        {
          case MO_LEX:
            o << "Lex => " << p->nvars;
            break;
          case MO_LEX2:
            o << "LexSmall => " << p->nvars;
            break;
          case MO_LEX4:
            o << "LexTiny => " << p->nvars;
            break;
          case MO_GREVLEX:
            o << "GRevLex => ";
            p_ones = true;
            break;
          case MO_GREVLEX2:
            o << "GRevLexSmall => ";
            p_ones = true;
            break;
          case MO_GREVLEX4:
            o << "GRevLexTiny => ";
            p_ones = true;
            break;
          case MO_GREVLEX_WTS:
            o << "GRevLex => ";
            break;
          case MO_GREVLEX2_WTS:
            o << "GRevLexSmall => ";
            break;
          case MO_GREVLEX4_WTS:
            o << "GRevLexTiny => ";
            break;
          case MO_REVLEX:
            o << "RevLex => " << p->nvars;
            break;
          case MO_WEIGHTS:
            o << "Weights => ";
            break;
          case MO_LAURENT:
            o << "GroupLex => " << p->nvars;
            break;
          case MO_LAURENT_REVLEX:
            o << "GroupRevLex => " << p->nvars;
            break;
          case MO_NC_LEX:
            o << "NCLex => " << p->nvars;
            break;
          case MO_POSITION_UP:
            o << "Position => Up";
            break;
          case MO_POSITION_DOWN:
            o << "Position => Down";
            break;
          default:
            o << "UNKNOWN";
            break;
        }
      if (p->wts != nullptr) { ::toString(o, p->nvars, p->wts); }
      else if (p_ones)
        {
          ::ones(o, p->nvars);
        }
    }
  o << "\n    }";
  return o.str();
}

// Many monomial ordering routines are in monordering.c
// Here are some that use c++ features, so cannot be there.
// @todo Make monordering.{h,c} into a c++ class.

static void write_row(std::vector<int> &grading,
                      int nvars,
                      int which,
                      int value)
{
  for (int i = 0; i < nvars; i++)
    if (i == which)
      grading.push_back(value);
    else
      grading.push_back(0);
}
static void write_weights(std::vector<int> &grading,
                          int nvars,
                          int firstvar,
                          int *wts,
                          int nwts)
// place nvars ints into grading:  0 ... 0 wts[0] wts[1] ... wts[nwts-1] 0 ....
// 0
// where wts[0] is in the 'firstvar' location.  If wts is NULL, treat it as the
// vector with nwts '1's.
{
  for (int i = 0; i < firstvar; i++) grading.push_back(0);
  if (wts == nullptr)
    for (int i = 0; i < nwts; i++) grading.push_back(1);
  else
    for (int i = 0; i < nwts; i++) grading.push_back(wts[i]);
  for (int i = firstvar + nwts; i < nvars; i++) grading.push_back(0);
}

bool monomialOrderingToMatrix(
    const struct MonomialOrdering& mo,
    std::vector<int>& mat,
    bool& base_is_revlex,
    int& component_direction,     // -1 is Down, +1 is Up, 0 is not present
    int& component_is_before_row  // -1 means: at the end. 0 means before the
                                  // order.
    // and r means considered before row 'r' of the matrix.
    )
{
  // a false return value means an error has occurred.
  int nvars = rawNumberOfVariables(&mo);
  base_is_revlex = true;
  enum LastBlock { LEX, REVLEX, WEIGHTS, NONE };
  LastBlock last = NONE;
  int nwts = 0;  // local var used in MO_WEIGHTS section
  int nrows = 0;
  int firstvar = 0;
  component_direction = 0;
  component_is_before_row =
      -2;                   // what should the default value be?  Probably: -1.
  size_t last_element = 0;  // The vector 'mat' will be resized back to this
                            // value if the last part of the order is lex or
                            // revlex.
  for (int i = 0; i < mo.len; i++)
    {
      mon_part p = mo.array[i];
      switch (p->type)
        {
          case MO_LEX:
          case MO_LEX2:
          case MO_LEX4:
            // printf("lex %d\n", p->nvars);
            last_element = mat.size();
            for (int j = 0; j < p->nvars; j++)
              {
                write_row(mat, nvars, firstvar + j, 1);
              }
            last = LEX;
            firstvar += p->nvars;
            nrows += p->nvars;
            break;
          case MO_GREVLEX:
          case MO_GREVLEX2:
          case MO_GREVLEX4:
            // printf("grevlex %d %ld\n", p->nvars, p->wts);
            write_weights(mat, nvars, firstvar, p->wts, p->nvars);
            last_element = mat.size();
            for (int j = p->nvars - 1; j >= 1; --j)
              {
                write_row(mat, nvars, firstvar + j, -1);
              }
            last = REVLEX;
            firstvar += p->nvars;
            nrows += p->nvars;
            break;
          case MO_GREVLEX_WTS:
          case MO_GREVLEX2_WTS:
          case MO_GREVLEX4_WTS:
            // printf("grevlex_wts %d %ld\n", p->nvars, p->wts);
            write_weights(mat, nvars, firstvar, p->wts, p->nvars);
            last_element = mat.size();
            for (int j = p->nvars - 1; j >= 1; --j)
              {
                write_row(mat, nvars, firstvar + j, -1);
              }
            last = REVLEX;
            firstvar += p->nvars;
            nrows += p->nvars;
            break;
          case MO_REVLEX:
            // printf("revlex %d\n", p->nvars);
            last_element = mat.size();
            for (int j = p->nvars - 1; j >= 0; --j)
              {
                write_row(mat, nvars, firstvar + j, -1);
              }
            last = REVLEX;
            firstvar += p->nvars;
            nrows += p->nvars;
            break;
          case MO_WEIGHTS:
            // printf("matsize= %d weights %d p->wts=%lu\n", mat.size(),
            // p->nvars, p->wts);
            nwts = (p->nvars > nvars ? nvars : p->nvars);
            write_weights(mat, nvars, 0, p->wts, nwts);
            nrows++;
            last_element = mat.size();
            last = WEIGHTS;
            break;
          case MO_LAURENT:
          case MO_LAURENT_REVLEX:
          case MO_NC_LEX:
            return false;
            break;
          case MO_POSITION_UP:
            component_direction = 1;
            component_is_before_row = nrows;
            break;
          case MO_POSITION_DOWN:
            component_direction = -1;
            component_is_before_row = nrows;
            break;
          default:
            // DO nothing
            break;
        }
    }
  if (last == LEX)
    {
      // last block was lex, so use lex tie-breaker
      mat.resize(last_element);
      if (nrows == component_is_before_row) component_is_before_row = -1;
      base_is_revlex = false;
    }
  else if (last == REVLEX)
    {
      // last block was revlex, so use revlex tie-breaker
      if (nrows == component_is_before_row) component_is_before_row = -1;
      mat.resize(last_element);
    }
  else
    {
      // last block is a weight vector, so use revlex as the tie-breaker.
      // nothing to change here.
    }
  return true;
}

int moIsLex(const MonomialOrdering *mo)
{
  // The monomial order is lex if what?
  // one lex block, no grevlex blocks, no weightvector blocks.
  // only: lex block and position blocks are allowed.
  int nlex = 0;
  int i;
  for (i = 0; i < mo->len; i++)
    {
      enum MonomialOrdering_type typ = mo->array[i]->type;
      switch (typ)
        {
          case MO_LEX:
          case MO_LEX2:
          case MO_LEX4:
            nlex++;
            break;
          case MO_POSITION_UP:
          case MO_POSITION_DOWN:
            break;
          default:
            return 0;
        }
    }
  return (nlex == 1);
}

int moIsGRevLex(const MonomialOrdering *mo)
{
  // The monomial order is lex if what?
  // one lex block, no grevlex blocks, no weightvector blocks.
  // only: lex block and position blocks are allowed.
  int ngrevlex = 0;
  int i;
  for (i = 0; i < mo->len; i++)
    {
      enum MonomialOrdering_type typ = mo->array[i]->type;
      switch (typ)
        {
          case MO_GREVLEX:
          case MO_GREVLEX2:
          case MO_GREVLEX4:
          case MO_GREVLEX_WTS:
          case MO_GREVLEX2_WTS:
          case MO_GREVLEX4_WTS:
            ngrevlex++;
            break;
          case MO_POSITION_UP:
          case MO_POSITION_DOWN:
            break;
          default:
            return 0;
        }
    }
  return (ngrevlex == 1);
}

int rawNumberOfVariables(const MonomialOrdering *mo)
{
  int i, sum = 0;
  for (i = 0; i < mo->len; i++)
    if (mo->array[i]->type != MO_WEIGHTS) sum += mo->array[i]->nvars;
  return sum;
}

M2_arrayint moGetWeightValues(const MonomialOrdering *mo)
{
  int nvars = rawNumberOfVariables(mo);
  // grab the first weight vector
  if (mo->len == 0) return nullptr;
  if (mo->array[0]->type == MO_WEIGHTS)
    {
      int i;
      M2_arrayint result = M2_makearrayint(nvars);
      int *wts = mo->array[0]->wts;
      for (i = 0; i < mo->array[0]->nvars; i++) result->array[i] = wts[i];
      for (; i < nvars; i++) result->array[i] = 0;
      return result;
    }
  return nullptr;
}

int rawNumberOfInvertibleVariables(const MonomialOrdering *mo)
{
  int i, sum = 0;
  for (i = 0; i < mo->len; i++)
    if (mo->array[i]->type == MO_LAURENT ||
        mo->array[i]->type == MO_LAURENT_REVLEX)
      sum += mo->array[i]->nvars;
  return sum;
}

M2_arrayint rawNonTermOrderVariables(const MonomialOrdering *mo)
// returns a list of the indices of those variables which are less than 1 in
// the given monomial order.
{
  int i, j, sum, nextvar;
  int nvars = rawNumberOfVariables(mo);
  int *gt = getmematomicvectortype(int, nvars);
  for (i = 0; i < nvars; i++)
    gt[i] =
        0;  // 0 means undecided, -1 means non term order, 1 means term order
  // Now we loop through the parts of the monomial order
  nextvar = 0;
  for (i = 0; i < mo->len; i++)
    {
      mon_part p = mo->array[i];
      switch (p->type)
        {
          case MO_LEX:
          case MO_LEX2:
          case MO_LEX4:
          case MO_GREVLEX:
          case MO_GREVLEX2:
          case MO_GREVLEX4:
          case MO_GREVLEX_WTS:
          case MO_GREVLEX2_WTS:
          case MO_GREVLEX4_WTS:
          case MO_LAURENT:
          case MO_NC_LEX:
            for (j = 0; j < p->nvars; j++, nextvar++)
              if (gt[nextvar] == 0) gt[nextvar] = 1;
            break;
          case MO_LAURENT_REVLEX:
          case MO_REVLEX:
            for (j = 0; j < p->nvars; j++, nextvar++)
              if (gt[nextvar] == 0) gt[nextvar] = -1;
            break;
          case MO_WEIGHTS:
            for (j = nextvar; j < p->nvars; j++)
              if (gt[j] == 0)
                {
                  if (p->wts[j] > 0)
                    gt[j] = 1;
                  else if (p->wts[j] < 0)
                    gt[j] = -1;
                }
            break;
          case MO_POSITION_UP:
          case MO_POSITION_DOWN:
            break;
        }
    }
  // At this point every variables' gt should be 1 or -1.
  sum = 0;
  for (i = 0; i < nvars; i++)
    {
      if (gt[i] == 0) INTERNAL_ERROR("gt[i] should not be 0");
      if (gt[i] < 0) sum++;
    }

  // Make an array of this length.
  M2_arrayint result = M2_makearrayint(sum);
  nextvar = 0;
  for (i = 0; i < nvars; i++)
    if (gt[i] < 0) result->array[nextvar++] = i;
  return result;
}

MonomialOrdering *rawLexMonomialOrdering(int nvars, int packing)
{
  MonomialOrdering *result;
  mon_part p;
  enum MonomialOrdering_type typ;

  if (packing == 2)
    typ = MO_LEX2;
  else if (packing == 4)
    typ = MO_LEX4;
  else
    typ = MO_LEX;

  p = mo_make(typ, nvars, nullptr);
  result = make_mon_order(1);
  result->array[0] = p;
  return result;
}

MonomialOrdering /* or null */ *rawGRevLexMonomialOrdering(M2_arrayint degs,
                                                           int packing)
{
  MonomialOrdering *result;
  mon_part p;
  enum MonomialOrdering_type typ;
  int *wts;
  int all_one = 1;
  unsigned int i;
  for (i = 0; i < degs->len; i++)
    if (degs->array[i] <= 0)
      {
        ERROR("grevlex: expected all degrees to be positive");
        return nullptr;
      }
    else if (degs->array[i] > 1)
      all_one = 0;

  if (all_one)
    {
      if (packing == 2)
        typ = MO_GREVLEX2;
      else if (packing == 4)
        typ = MO_GREVLEX4;
      else
        typ = MO_GREVLEX;
      wts = nullptr;
    }
  else
    {
      if (packing == 2)
        typ = MO_GREVLEX2_WTS;
      else if (packing == 4)
        typ = MO_GREVLEX4_WTS;
      else
        typ = MO_GREVLEX_WTS;
      wts = degs->array;
    }

  p = mo_make(typ, degs->len, wts);
  result = make_mon_order(1);
  result->array[0] = p;
  return result;
}

MonomialOrdering *rawRevLexMonomialOrdering(int nvars)
{
  mon_part p = mo_make(MO_REVLEX, nvars, nullptr);
  MonomialOrdering *result = make_mon_order(1);
  result->array[0] = p;
  return result;
}

MonomialOrdering *rawWeightsMonomialOrdering(M2_arrayint wts)
{
  mon_part p = mo_make(MO_WEIGHTS, wts->len, wts->array);
  MonomialOrdering *result = make_mon_order(1);
  result->array[0] = p;
  return result;
}
MonomialOrdering *rawGroupLexMonomialOrdering(int nvars)
{
  mon_part p = mo_make(MO_LAURENT, nvars, nullptr);
  MonomialOrdering *result = make_mon_order(1);
  result->array[0] = p;
  return result;
}
MonomialOrdering *rawGroupRevLexMonomialOrdering(int nvars)
{
  mon_part p = mo_make(MO_LAURENT_REVLEX, nvars, nullptr);
  MonomialOrdering *result = make_mon_order(1);
  result->array[0] = p;
  return result;
}
MonomialOrdering *rawNClexMonomialOrdering(int nvars)
{
  mon_part p = mo_make(MO_NC_LEX, nvars, nullptr);
  MonomialOrdering *result = make_mon_order(1);
  result->array[0] = p;
  return result;
}
MonomialOrdering *rawPositionMonomialOrdering(M2_bool up_or_down)
{
  mon_part p =
      mo_make((up_or_down ? MO_POSITION_UP : MO_POSITION_DOWN), 0, nullptr);
  MonomialOrdering *result = make_mon_order(1);
  result->array[0] = p;
  return result;
}

MonomialOrdering *rawJoinMonomialOrdering(engine_RawMonomialOrderingArray M)
{
  MonomialOrdering *result;
  const MonomialOrdering *mo;
  int i, j, sum, next;
  int nvars_so_far = 0;

  //  sum = 0;
  //  for (i=0; i<M->len; i++)
  //    sum += M->array[i]->len;

  sum = 0;
  for (i = 0; i < M->len; i++)
    {
      mo = M->array[i];
      for (j = 0; j < mo->len; j++)
        if (is_good(mo->array[j])) sum++;
    }

  result = make_mon_order(sum);
  next = 0;
  for (i = 0; i < M->len; i++)
    {
      mo = M->array[i];
      for (j = 0; j < mo->len; j++)
        {
          mon_part p = mo->array[j];
          if (!is_good(p)) continue;
          if (p->type != MO_WEIGHTS)
            nvars_so_far += p->nvars;
          else
            {
              /* Shift the weights over by nvars_so_far */
              mon_part q = mo_make(MO_WEIGHTS, nvars_so_far + p->nvars, nullptr);
              q->wts = getmemvectortype(int, q->nvars);
              for (j = 0; j < nvars_so_far; j++) q->wts[j] = 0;
              for (; j < q->nvars; j++) q->wts[j] = p->wts[j - nvars_so_far];
              p = q;
            }
          result->array[next++] = p;
        }
    }
  return result;
}

MonomialOrdering *rawProductMonomialOrdering(engine_RawMonomialOrderingArray M)
{
  MonomialOrdering *result;
  int i, j, sum, next, offset;
  sum = 0;
  for (i = 0; i < M->len; i++) sum += M->array[i]->len;

  result = make_mon_order(sum);
  next = 0;
  offset = 0;
  for (i = 0; i < M->len; i++)
    {
      int nvars = rawNumberOfVariables(M->array[i]);
      MonomialOrdering *mo = M2_mo_offset(M->array[i], offset);
      for (j = 0; j < mo->len; j++) result->array[next++] = mo->array[j];
      offset += nvars;
    }
  return result;
}

M2_string intarray_to_string(int len, int *p)
{
  int i;
  char s[200];
  M2_string result = M2_tostring("{");
  for (i = 0; i < len; i++)
    {
      if (i > 0) result = M2_join(result, M2_tostring(","));
      sprintf(s, "%d", p[i]);
      result = M2_join(result, M2_tostring(s));
    }
  result = M2_join(result, M2_tostring("}"));
  return result;
}

M2_string ones_to_string(int len)
{
  int i;
  char s[200];
  M2_string one;
  M2_string result = M2_tostring("{");
  sprintf(s, "1");
  one = M2_tostring(s);
  for (i = 0; i < len; i++)
    {
      if (i > 0) result = M2_join(result, M2_tostring(","));
      result = M2_join(result, one);
    }
  result = M2_join(result, M2_tostring("}"));
  return result;
}

unsigned int rawMonomialOrderingHash(const MonomialOrdering *mo)
{
  return mo->_hash;
}

M2_string IM2_MonomialOrdering_to_string(const MonomialOrdering *mo)
{
  int i;
  char s[200];
  int p_ones = 0;
  M2_string result = M2_tostring("MonomialOrder => {");
  for (i = 0; i < mo->len; i++)
    {
      mon_part p = mo->array[i];
      p_ones = 0;
      if (i == 0)
        result = M2_join(result, M2_tostring("\n    "));
      else
        result = M2_join(result, M2_tostring(",\n    "));
      switch (p->type)
        {
          case MO_LEX:
            sprintf(s, "Lex => %d", p->nvars);
            break;
          case MO_LEX2:
            sprintf(s, "LexSmall => %d", p->nvars);
            break;
          case MO_LEX4:
            sprintf(s, "LexTiny => %d", p->nvars);
            break;
          case MO_GREVLEX:
            sprintf(s, "GRevLex => ");
            p_ones = 1;
            break;
          case MO_GREVLEX2:
            sprintf(s, "GRevLexSmall => ");
            p_ones = 1;
            break;
          case MO_GREVLEX4:
            sprintf(s, "GRevLexTiny => ");
            p_ones = 1;
            break;
          case MO_GREVLEX_WTS:
            sprintf(s, "GRevLex => ");
            break;
          case MO_GREVLEX2_WTS:
            sprintf(s, "GRevLexSmall => ");
            break;
          case MO_GREVLEX4_WTS:
            sprintf(s, "GRevLexTiny => ");
            break;
          case MO_REVLEX:
            sprintf(s, "RevLex => %d", p->nvars);
            break;
          case MO_WEIGHTS:
            sprintf(s, "Weights => ");
            break;
          case MO_LAURENT:
            sprintf(s, "GroupLex => %d", p->nvars);
            break;
          case MO_LAURENT_REVLEX:
            sprintf(s, "GroupRevLex => %d", p->nvars);
            break;
          case MO_NC_LEX:
            sprintf(s, "NCLex => %d", p->nvars);
            break;
          case MO_POSITION_UP:
            sprintf(s, "Position => Up");
            break;
          case MO_POSITION_DOWN:
            sprintf(s, "Position => Down");
            break;
          default:
            sprintf(s, "UNKNOWN");
            break;
        }
      result = M2_join(result, M2_tostring(s));
      if (p->wts != nullptr)
        result = M2_join(result, intarray_to_string(p->nvars, p->wts));
      else if (p_ones)
        result = M2_join(result, ones_to_string(p->nvars));
    }
  result = M2_join(result, M2_tostring("\n    }"));
  return result;
}

#include "overflow.hpp"

std::vector<bool> laurentVariables(const MonomialOrder* mo)
{
  std::vector<bool> result;
  for (auto i = 0; i < mo->nvars; ++i)
    result.push_back(mo->is_laurent[i] == 1);
  return result;
}
/* TODO:
   -- negative exponent versions need to be included (at least for MO_LEX)
   -- non-commutative blocks should be added in
*/

static void mo_block_revlex(struct mo_block *b, int nvars)
{
  b->typ = MO_REVLEX;
  b->nvars = nvars;
  b->nslots = nvars;
  b->first_exp = 0;  /* will be set later */
  b->first_slot = 0; /* will be set later */
  b->nweights = 0;
  b->weights = nullptr;
}

static void mo_block_grevlex(struct mo_block *b, int nvars)
{
  b->typ = MO_GREVLEX;
  b->nvars = nvars;
  b->nslots = nvars;
  b->first_exp = 0;  /* will be set later */
  b->first_slot = 0; /* will be set later */
  b->nweights = 0;
  b->weights = nullptr;
}

static void mo_block_grevlex2(struct mo_block *b, int nvars)
{
  b->typ = MO_GREVLEX2;
  b->nvars = nvars;
  b->nslots = (nvars + 1) / 2; /* 2 per word */
  b->first_exp = 0;            /* will be set later */
  b->first_slot = 0;           /* will be set later */
  b->nweights = 0;
  b->weights = nullptr;
}

static void mo_block_grevlex4(struct mo_block *b, int nvars)
{
  b->typ = MO_GREVLEX4;
  b->nvars = nvars;
  b->nslots = (nvars + 3) / 4; /* 4 per word */
  b->first_exp = 0;            /* will be set later */
  b->first_slot = 0;           /* will be set later */
  b->nweights = 0;
  b->weights = nullptr;
}

static void mo_block_grevlex_wts(struct mo_block *b, int nvars)
{
  b->typ = MO_GREVLEX_WTS;
  b->nvars = nvars;
  b->nslots = nvars;
  b->first_exp = 0;  /* will be set later */
  b->first_slot = 0; /* will be set later */
  b->nweights = nvars;
  b->weights = nullptr; /* will be set later */
}

static void mo_block_grevlex2_wts(struct mo_block *b, int nvars)
{
  b->typ = MO_GREVLEX2_WTS;
  b->nvars = nvars;
  b->nslots = (nvars + 1) / 2; /* 2 per word */
  b->first_exp = 0;            /* will be set later */
  b->first_slot = 0;           /* will be set later */
  b->nweights = nvars;
  b->weights = nullptr; /* will be set later */
}

static void mo_block_grevlex4_wts(struct mo_block *b, int nvars)
{
  b->typ = MO_GREVLEX4_WTS;
  b->nvars = nvars;
  b->nslots = (nvars + 3) / 4; /* 4 per word */
  b->first_exp = 0;            /* will be set later */
  b->first_slot = 0;           /* will be set later */
  b->nweights = nvars;
  b->weights = nullptr; /* will be set later */
}

static void mo_block_lex(struct mo_block *b, int nvars)
{
  b->typ = MO_LEX;
  b->nvars = nvars;
  b->nslots = nvars;
  b->first_exp = 0;  /* will be set later */
  b->first_slot = 0; /* will be set later */
  b->nweights = 0;
  b->weights = nullptr;
}

static void mo_block_lex2(struct mo_block *b, int nvars)
{
  b->typ = MO_LEX2;
  b->nvars = nvars;
  b->nslots = (nvars + 1) / 2; /* 2 per word */
  b->first_exp = 0;            /* will be set later */
  b->first_slot = 0;           /* will be set later */
  b->nweights = 0;
  b->weights = nullptr;
}

static void mo_block_lex4(struct mo_block *b, int nvars)
{
  b->typ = MO_LEX4;
  b->nvars = nvars;
  b->nslots = (nvars + 3) / 4; /* 4 per word */
  b->first_exp = 0;            /* will be set later */
  b->first_slot = 0;           /* will be set later */
  b->nweights = 0;
  b->weights = nullptr;
}

static void mo_block_group_lex(struct mo_block *b, int nvars)
{
  b->typ = MO_LAURENT;
  b->nvars = nvars;
  b->nslots = nvars;
  b->first_exp = 0;  /* will be set later */
  b->first_slot = 0; /* will be set later */
  b->nweights = 0;
  b->weights = nullptr;
}

static void mo_block_group_revlex(struct mo_block *b, int nvars)
{
  b->typ = MO_LAURENT_REVLEX;
  b->nvars = nvars;
  b->nslots = nvars;
  b->first_exp = 0;  /* will be set later */
  b->first_slot = 0; /* will be set later */
  b->nweights = 0;
  b->weights = nullptr;
}

static void mo_block_wt_function(struct mo_block *b, int nvars, deg_t *wts)
{
  b->typ = MO_WEIGHTS;
  b->nvars = 0;
  b->nslots = 1;
  b->first_exp = 0;
  b->first_slot = 0; /* will be set later */
  b->nweights = nvars;
  b->weights = wts;
}

MonomialOrder *monomialOrderMake(const MonomialOrdering *mo)
{
  MonomialOrder *result;
  int i, j, nv, this_block;
  deg_t *wts = nullptr;
  /* Determine the number of variables, the number of blocks, and the location
     of the component */
  int nblocks = 0;
  int nvars = 0;
  int hascomponent = 0;
  for (i = 0; i < mo->len; i++)
    {
      struct mon_part_rec_ *t = mo->array[i];
      nblocks++;
      if (t->type == MO_POSITION_DOWN || t->type == MO_POSITION_UP)
        hascomponent++;
      else if (t->type == MO_NC_LEX)
        {
          // Currently, do nothing.
        }
      if (t->type != MO_WEIGHTS) nvars += t->nvars;
    }
  nblocks -= hascomponent;

  /* Now create the blocks, and fill them in. Also fill in the deg vector */
  result = getmemstructtype(MonomialOrder *);
  result->nvars = nvars;
  result->nslots = 0;
  result->nblocks = nblocks;
  result->blocks =
      (struct mo_block *)getmem(nblocks * sizeof(result->blocks[0]));
  result->degs = (deg_t *)getmem_atomic(nvars * sizeof(result->degs[0]));
  if (hascomponent == 0) result->nblocks_before_component = nblocks;

  this_block = 0;
  nvars = 0;
  for (i = 0; i < mo->len; i++)
    {
      struct mon_part_rec_ *t = mo->array[i];
      if (t->type != MO_WEIGHTS)
        {
          if (t->wts == nullptr)
            for (j = 0; j < t->nvars; j++) result->degs[nvars++] = 1;
          else
            for (j = 0; j < t->nvars; j++) result->degs[nvars++] = t->wts[j];
        }
      else
        {
          wts = (deg_t *)getmem_atomic(t->nvars * sizeof(wts[0]));
          for (j = 0; j < t->nvars; j++) wts[j] = t->wts[j];
        }
      switch (t->type)
        {
          case MO_REVLEX:
            mo_block_revlex(result->blocks + this_block++, t->nvars);
            break;
          case MO_GREVLEX:
            mo_block_grevlex(result->blocks + this_block++, t->nvars);
            break;
          case MO_GREVLEX2:
            mo_block_grevlex2(result->blocks + this_block++, t->nvars);
            break;
          case MO_GREVLEX4:
            mo_block_grevlex4(result->blocks + this_block++, t->nvars);
            break;
          case MO_GREVLEX_WTS:
            mo_block_grevlex_wts(result->blocks + this_block++, t->nvars);
            break;
          case MO_GREVLEX2_WTS:
            mo_block_grevlex2_wts(result->blocks + this_block++, t->nvars);
            break;
          case MO_GREVLEX4_WTS:
            mo_block_grevlex4_wts(result->blocks + this_block++, t->nvars);
            break;
          case MO_LEX:
            mo_block_lex(result->blocks + this_block++, t->nvars);
            break;
          case MO_LEX2:
            mo_block_lex2(result->blocks + this_block++, t->nvars);
            break;
          case MO_LEX4:
            mo_block_lex4(result->blocks + this_block++, t->nvars);
            break;
          case MO_WEIGHTS:
            // if extra weight values are given (more than "nvars", ignore the
            // rest.
            mo_block_wt_function(
                result->blocks + this_block++,
                (t->nvars <= result->nvars ? t->nvars : result->nvars),
                wts);
            break;
          case MO_LAURENT:
            mo_block_group_lex(result->blocks + this_block++, t->nvars);
            break;
          case MO_LAURENT_REVLEX:
            mo_block_group_revlex(result->blocks + this_block++, t->nvars);
            break;
          case MO_NC_LEX:
            /* MES */
            break;
          case MO_POSITION_UP:
            if (--hascomponent == 0)
              {
                // Set the information about the component
                result->component_up = 1;
                result->nblocks_before_component = this_block;
              }
            //  mo_block_position_up(result->blocks + this_block);
            break;
          case MO_POSITION_DOWN:
            if (--hascomponent == 0)
              {
                // Set the information about the component
                result->component_up = 0;
                result->nblocks_before_component = this_block;
              }
            //  mo_block_position_down(result->blocks + this_block);
            break;
        }
    }

  /* Go back and fill in the 'slots' information */
  /* Now fix the first_exp, first_slot values, and also result->{nslots,nvars};
   */
  nv = 0;
  result->nslots = 0;
  result->nslots_before_component = 0;
  for (i = 0; i < nblocks; i++)
    {
      enum MonomialOrdering_type typ = result->blocks[i].typ;

      result->blocks[i].first_exp = nv;
      result->blocks[i].first_slot = result->nslots;
      nv += result->blocks[i].nvars;
      result->nslots += result->blocks[i].nslots;

      if (typ == MO_WEIGHTS)
        {
          result->blocks[i].first_exp = 0;

          /* divide the wt vector by the degree vector */
          for (j = 0; j < result->blocks[i].nvars; j++)
            safe::div_by(result->blocks[i].weights[j], result->degs[j]);
          ;
        }
      else if (typ == MO_GREVLEX_WTS || typ == MO_GREVLEX2_WTS ||
               typ == MO_GREVLEX4_WTS)
        {
          result->blocks[i].weights =
              result->degs + result->blocks[i].first_exp;
        }

      if (i == result->nblocks_before_component - 1)
        {
          result->nslots_before_component = result->nslots;
        }
    }

  /* Set is_laurent */
  result->is_laurent = (int *)getmem_atomic(result->nvars * sizeof(int));
  for (i = 0; i < result->nvars; i++) result->is_laurent[i] = 0;

  for (i = 0; i < result->nblocks; i++)
    if (result->blocks[i].typ == MO_LAURENT ||
        result->blocks[i].typ == MO_LAURENT_REVLEX)
      {
        for (j = 0; j < result->blocks[i].nvars; j++)
          result->is_laurent[result->blocks[i].first_exp + j] = 1;
      }

  return result;
}

extern void monomialOrderFree(MonomialOrder *mo) {}
static void MO_pack4(int nvars, const int *expon, int *slots)
{
  int32_t i;
  if (nvars == 0) return;
  while (1)
    {
      i = safe::fits_7(*expon++) << 24;
      if (--nvars == 0) break;
      i |= safe::fits_7(*expon++) << 16;
      if (--nvars == 0) break;
      i |= safe::fits_7(*expon++) << 8;
      if (--nvars == 0) break;
      i |= safe::fits_7(*expon++);
      if (--nvars == 0) break;
      *slots++ = i;
    }
  *slots++ = i;
}

static void MO_pack2(int nvars, const int *expon, int *slots)
{
  int32_t i;
  if (nvars == 0) return;
  while (1)
    {
      i = safe::fits_15(*expon++) << 16;
      if (--nvars == 0) break;
      i |= safe::fits_15(*expon++);
      if (--nvars == 0) break;
      *slots++ = i;
    }
  *slots++ = i;
}

static void MO_unpack4(int nvars, const int *slots, int *expon)
{
  int32_t i;
  if (nvars == 0) return;
  while (1)
    {
      i = *slots++;
      *expon++ = (i >> 24);
      if (--nvars == 0) break;
      *expon++ = (i >> 16) & 0x7f;
      if (--nvars == 0) break;
      *expon++ = (i >> 8) & 0x7f;
      if (--nvars == 0) break;
      *expon++ = i & 0x7f;
      if (--nvars == 0) break;
    }
}

static void MO_unpack2(int nvars, const int *slots, int *expon)
{
  int32_t i;
  if (nvars == 0) return;
  while (1)
    {
      i = *slots++;
      *expon++ = i >> 16;
      if (--nvars == 0) break;
      *expon++ = i & 0x7fff;
      if (--nvars == 0) break;
    }
}

std::vector<int> MonomialOrder::overflow_flags() const
{
  std::vector<int> overflow;
  enum overflow_type flag;
  int i = 0;
  int k = 0;
  for (; i < nblocks; i++)
    {
      mo_block *b = &blocks[i];
      switch (b->typ)
        {
          case MO_REVLEX:
          case MO_WEIGHTS:
          case MO_LAURENT:
          case MO_LAURENT_REVLEX:
          case MO_NC_LEX:
            flag = OVER;
            goto fillin;
          case MO_POSITION_UP:
          case MO_POSITION_DOWN:
            ERROR(
                "internal error - MO_POSITION_DOWN or MO_POSITION_UP "
                "encountered");
            assert(0);
            break;
          case MO_LEX:
          case MO_GREVLEX:
          case MO_GREVLEX_WTS:
            flag = OVER1;
            goto fillin;
          case MO_LEX2:
          case MO_GREVLEX2:
          case MO_GREVLEX2_WTS:
            flag = OVER2;
            goto fillin;
          case MO_LEX4:
          case MO_GREVLEX4:
          case MO_GREVLEX4_WTS:
            flag = OVER4;
            goto fillin;
          fillin:
            assert(b->first_slot == k);
            for (int p = b->nslots; p > 0; p--)
              {
                assert(k < nslots);
                overflow.push_back(flag);
                k++;
              }
            break;
          default:
            ERROR("internal error - missing case");
            assert(0);
            break;
        }
    }
  assert(k == nslots);
  return overflow;
}



void monomialOrderEncodeFromActualExponents(const MonomialOrder *mo,
                                            const_exponents expon,
                                            monomial result_psums)
/* Given 'expon', compute the encoded partial sums value */
{
  if (mo == nullptr) return;
  int *tmpexp = static_cast<int *>(alloca((mo->nvars + 1) * sizeof(int)));
  int i, j, nvars, s;
  int *p1;
  deg_t *degs;
  struct mo_block *b = mo->blocks;
  int nblocks = mo->nblocks;
  const int *e = expon;
  int *p = result_psums;
  for (i = nblocks; i > 0; --i, b++) switch (b->typ)
      {
        case MO_LEX:
        case MO_LAURENT:
          nvars = b->nvars;
          for (j = 0; j < nvars; j++) *p++ = *e++;
          break;
        case MO_REVLEX:
        case MO_LAURENT_REVLEX:
          nvars = b->nvars;
          for (j = 0; j < nvars; j++) *p++ = safe::minus(*e++);
          break;
        case MO_GREVLEX:
          nvars = b->nvars;
          p += b->nslots;
          p1 = p;
          *--p1 = *e++;
          for (j = 1; j < nvars; j++)
            {
              --p1;
              *p1 = safe::add(*e++, p1[1]);
            }
          break;
        case MO_GREVLEX_WTS:
          nvars = b->nvars;
          degs = mo->degs + b->first_exp;
          p += b->nslots;
          p1 = p;
          *--p1 = safe::mult(*e++, *degs++);
          for (j = 1; j < nvars; j++)
            {
              --p1;
              int tmp = safe::mult(*e++, *degs++);
              *p1 = safe::add(tmp, p1[1]);
            }
          break;
        case MO_GREVLEX4:
          nvars = b->nvars;
          p1 = tmpexp + b->nvars;
          *--p1 = *e++;
          for (j = 1; j < nvars; j++)
            {
              --p1;
              *p1 = safe::add(*e++, p1[1]);
            }
          MO_pack4(nvars, p1, p);
          p += b->nslots;
          break;
        case MO_GREVLEX4_WTS:
          nvars = b->nvars;
          degs = mo->degs + b->first_exp;
          p1 = tmpexp + b->nvars;
          *--p1 = safe::mult(*e++, *degs++);
          for (j = 1; j < nvars; j++)
            {
              --p1;
              int tmp = safe::mult(*e++, *degs++);
              *p1 = safe::add(tmp, p1[1]);
            }
          MO_pack4(nvars, p1, p);
          p += b->nslots;
          break;
        case MO_GREVLEX2:
          nvars = b->nvars;
          p1 = tmpexp + b->nvars;
          *--p1 = *e++;
          for (j = 1; j < nvars; j++)
            {
              --p1;
              *p1 = safe::add(*e++, p1[1]);
            }
          MO_pack2(nvars, p1, p);
          p += b->nslots;
          break;
        case MO_GREVLEX2_WTS:
          nvars = b->nvars;
          degs = mo->degs + b->first_exp;
          p1 = tmpexp + b->nvars;
          *--p1 = safe::mult(*e++, *degs++);
          for (j = 1; j < nvars; j++)
            {
              --p1;
              int tmp = safe::mult(*e++, *degs++);
              *p1 = safe::add(tmp, p1[1]);
            }
          MO_pack2(nvars, p1, p);
          p += b->nslots;
          break;
        case MO_LEX4:
          nvars = b->nvars;
          MO_pack4(nvars, e, p);
          p += b->nslots;
          e += nvars;
          break;
        case MO_LEX2:
          nvars = b->nvars;
          MO_pack2(nvars, e, p);
          p += b->nslots;
          e += nvars;
          break;
        case MO_WEIGHTS:
          if (b->nweights == 0)
            {
              s = 0;
            }
          else
            {
              s = safe::mult(b->weights[0], expon[0]);
              for (j = 1; j < b->nweights; j++)
                s = safe::add(s, safe::mult(b->weights[j], expon[j]));
            }
          *p++ = s;
          break;
        case MO_POSITION_UP:
        case MO_POSITION_DOWN:
          /* nothing to do here */
          break;
        case MO_NC_LEX:
          /* nothing to do here */
          break;
      }
}

void monomialOrderDecodeToActualExponents(const MonomialOrder *mo,
                                          const_monomial psums,
                                          exponents_t expon)
{
  if (mo == nullptr) return;
  int *tmpexp = static_cast<int *>(alloca((mo->nvars + 1) * sizeof(int)));
  int i, j, nvars;
  deg_t *degs = mo->degs;
  deg_t *d;
  struct mo_block *b = mo->blocks;
  int nblocks = mo->nblocks;
  int *e = expon;
  const int *p = psums;
  for (i = nblocks; i > 0; --i, b++) switch (b->typ)
      {
        case MO_LEX:
        case MO_LAURENT:
          nvars = b->nvars;
          p = psums + b->first_slot;
          e = expon + b->first_exp;
          for (j = 0; j < nvars; j++) *e++ = *p++;
          break;
        case MO_REVLEX:
        case MO_LAURENT_REVLEX:
          nvars = b->nvars;
          p = psums + b->first_slot;
          e = expon + b->first_exp;
          for (j = 0; j < nvars; j++) *e++ = safe::minus(*p++);
          break;
        case MO_GREVLEX:
          nvars = b->nvars;
          p = psums + b->first_slot + nvars - 1;
          e = expon + b->first_exp;
          *e++ = *p--;
          for (j = nvars - 1; j >= 1; --j, --p) *e++ = safe::sub(*p, p[1]);
          break;
        case MO_GREVLEX_WTS:
          nvars = b->nvars;
          d = degs + b->first_exp;
          p = psums + b->first_slot + nvars - 1;
          e = expon + b->first_exp;
          *e++ = *p-- / *d++;
          for (j = nvars - 1; j >= 1; --j, --p)
            *e++ = safe::sub(*p, p[1]) / *d++;
          break;
        case MO_GREVLEX4:
          nvars = b->nvars;
          MO_unpack4(nvars, psums + b->first_slot, tmpexp);
          p = tmpexp + nvars - 1;
          e = expon + b->first_exp;
          *e++ = *p--;
          for (j = nvars - 1; j >= 1; --j, --p) *e++ = safe::sub(*p, p[1]);
          break;
        case MO_GREVLEX4_WTS:
          nvars = b->nvars;
          d = degs + b->first_exp;
          MO_unpack4(nvars, psums + b->first_slot, tmpexp);
          p = tmpexp + nvars - 1;
          e = expon + b->first_exp;
          *e++ = *p-- / *d++;
          for (j = nvars - 1; j >= 1; --j, --p)
            *e++ = safe::sub(*p, p[1]) / *d++;
          break;
        case MO_GREVLEX2:
          nvars = b->nvars;
          MO_unpack2(nvars, psums + b->first_slot, tmpexp);
          p = tmpexp + nvars - 1;
          e = expon + b->first_exp;
          *e++ = *p--;
          for (j = nvars - 1; j >= 1; --j, --p) *e++ = safe::sub(*p, p[1]);
          break;
        case MO_GREVLEX2_WTS:
          nvars = b->nvars;
          d = degs + b->first_exp;
          MO_unpack2(nvars, psums + b->first_slot, tmpexp);
          p = tmpexp + nvars - 1;
          e = expon + b->first_exp;
          *e++ = *p-- / *d++;
          for (j = nvars - 1; j >= 1; --j, --p)
            *e++ = safe::sub(*p, p[1]) / *d++;
          break;
        case MO_LEX4:
          nvars = b->nvars;
          e = expon + b->first_exp;
          MO_unpack4(nvars, psums + b->first_slot, e);
          break;
        case MO_LEX2:
          nvars = b->nvars;
          e = expon + b->first_exp;
          MO_unpack2(nvars, psums + b->first_slot, e);
          break;
        case MO_WEIGHTS:
          break;
        case MO_POSITION_UP:
        case MO_POSITION_DOWN:
          /* should not occur, but do nothing in any case */
          break;
        case MO_NC_LEX:
          /* nothing to do here */
          break;
      }
}

MonomialOrdering *MonomialOrderings::Lex(int nvars)
{
  return rawLexMonomialOrdering(nvars, 1);
}
MonomialOrdering *MonomialOrderings::Lex2(int nvars)
{
  return rawLexMonomialOrdering(nvars, 2);
}
MonomialOrdering *MonomialOrderings::Lex4(int nvars)
{
  return rawLexMonomialOrdering(nvars, 4);
}
MonomialOrdering *MonomialOrderings::GRevLex(int nvars)
{
  std::vector<int> w;
  for (int i = 0; i < nvars; i++) w.push_back(1);
  return GRevLex(w);
}
MonomialOrdering *MonomialOrderings::GRevLex2(int nvars)
{
  std::vector<int> w;
  for (int i = 0; i < nvars; i++) w.push_back(1);
  return GRevLex2(w);
}
MonomialOrdering *MonomialOrderings::GRevLex4(int nvars)
{
  std::vector<int> w;
  for (int i = 0; i < nvars; i++) w.push_back(1);
  return GRevLex4(w);
}
MonomialOrdering *MonomialOrderings::GRevLex(const std::vector<int> &wts)
{
  return GRevLex(wts, 1);
}
MonomialOrdering *MonomialOrderings::GRevLex2(const std::vector<int> &wts)
{
  return GRevLex(wts, 2);
}
MonomialOrdering *MonomialOrderings::GRevLex4(const std::vector<int> &wts)
{
  return GRevLex(wts, 4);
}
MonomialOrdering *MonomialOrderings::RevLex(int nvars)
{
  return rawRevLexMonomialOrdering(nvars);
}
MonomialOrdering *MonomialOrderings::Weights(const std::vector<int> &wts)
{
  mon_part p = mo_make(MO_WEIGHTS, wts.size(), wts.data());
  MonomialOrdering *result = make_mon_order(1);
  result->array[0] = p;
  return result;
}
MonomialOrdering *MonomialOrderings::GroupLex(int nvars)
{
  return rawGroupLexMonomialOrdering(nvars);
}
MonomialOrdering *MonomialOrderings::GroupRevLex(int nvars)
{
  return rawGroupRevLexMonomialOrdering(nvars);
}
MonomialOrdering *MonomialOrderings::PositionUp()
{
  return rawPositionMonomialOrdering(true);
}
MonomialOrdering *MonomialOrderings::PositionDown()
{
  return rawPositionMonomialOrdering(false);
}

MonomialOrdering *MonomialOrderings::GRevLex(const std::vector<int> &degs,
                                             int packing)
{
  MonomialOrdering *result;
  mon_part p;
  enum MonomialOrdering_type typ;
  const int *wts;
  bool all_one = true;
  for (int i = 0; i < degs.size(); i++)
    if (degs[i] <= 0)
      {
        ERROR("grevlex: expected all degrees to be positive");
        return nullptr;
      }
    else if (degs[i] > 1)
      all_one = false;

  if (all_one)
    {
      if (packing == 2)
        typ = MO_GREVLEX2;
      else if (packing == 4)
        typ = MO_GREVLEX4;
      else
        typ = MO_GREVLEX;
      wts = nullptr;
    }
  else
    {
      if (packing == 2)
        typ = MO_GREVLEX2_WTS;
      else if (packing == 4)
        typ = MO_GREVLEX4_WTS;
      else
        typ = MO_GREVLEX_WTS;
      wts = degs.data();
    }

  p = mo_make(typ, degs.size(), wts);
  result = make_mon_order(1);
  result->array[0] = p;
  return result;
}

MonomialOrdering *MonomialOrderings::join(
    const std::vector<MonomialOrdering *> &M)
{
  MonomialOrdering *result;
  const MonomialOrdering *mo;
  int i, j, sum, next;
  int nvars_so_far = 0;

  sum = 0;
  for (i = 0; i < M.size(); i++)
    {
      mo = M[i];
      for (j = 0; j < mo->len; j++)
        if (is_good(mo->array[j])) sum++;
    }

  result = make_mon_order(sum);
  next = 0;
  for (i = 0; i < M.size(); i++)
    {
      mo = M[i];
      for (j = 0; j < mo->len; j++)
        {
          mon_part p = mo->array[j];
          if (!is_good(p)) continue;
          if (p->type != MO_WEIGHTS)
            nvars_so_far += p->nvars;
          else
            {
              /* Shift the weights over by nvars_so_far */
              mon_part q = mo_make(MO_WEIGHTS, nvars_so_far + p->nvars, nullptr);
              q->wts = getmemvectortype(int, q->nvars);
              for (j = 0; j < nvars_so_far; j++) q->wts[j] = 0;
              for (; j < q->nvars; j++) q->wts[j] = p->wts[j - nvars_so_far];
              p = q;
            }
          result->array[next++] = p;
        }
    }
  return result;
}

MonomialOrdering *MonomialOrderings::product(
    const std::vector<MonomialOrdering *> &M)
{
  MonomialOrdering *result;
  int i, j, sum, next, offset;
  sum = 0;
  for (i = 0; i < M.size(); i++) sum += M[i]->len;

  result = make_mon_order(sum);
  next = 0;
  offset = 0;
  for (i = 0; i < M.size(); i++)
    {
      int nvars = rawNumberOfVariables(M[i]);
      MonomialOrdering *mo = M2_mo_offset(M[i], offset);
      for (j = 0; j < mo->len; j++) result->array[next++] = mo->array[j];
      offset += nvars;
    }
  return result;
}


// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e "
// indent-tabs-mode: nil
// End:

