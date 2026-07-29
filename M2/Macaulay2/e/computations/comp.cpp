// Copyright 2004 Michael E. Stillman

#include "computations/comp.hpp"

#include <stdio.h>

#include "buffer.hpp"
#include "error.h"
#include "exceptions.hpp"
#include "finalize.hpp"

namespace {
std::atomic<std::size_t> nextThreadTag{1};
}

std::size_t Computation::currentThreadTag()
{
  static thread_local std::size_t tag =
      nextThreadTag.fetch_add(1, std::memory_order_relaxed);
  return tag;
}

Computation::ThreadGuard::ThreadGuard(Computation *C)
    : mComputation(C), mOtherThread(0), mOk(false)
{
  std::size_t self = currentThreadTag();
  std::size_t owner = 0;  // expected value: unowned
  if (C->mOwner.compare_exchange_strong(
          owner, self, std::memory_order_acq_rel, std::memory_order_acquire))
    {
      // we just claimed it
      C->mOwnerDepth = 1;
      mOk = true;
    }
  else if (owner == self)
    {
      // already ours: a nested engine call on the same thread
      C->mOwnerDepth++;
      mOk = true;
    }
  else
    {
      // another thread is inside this computation
      mOtherThread = owner;
    }
}

Computation::ThreadGuard::~ThreadGuard()
{
  if (!mOk) return;
  if (--mComputation->mOwnerDepth == 0)
    mComputation->mOwner.store(0, std::memory_order_release);
}

void Computation::ThreadGuard::reportConflict() const
{
  ERROR(
      "engine computation already in use by thread %lu (this is thread %lu): "
      "Groebner basis and resolution computations are not thread safe",
      static_cast<unsigned long>(mOtherThread),
      static_cast<unsigned long>(currentThreadTag()));
}

Computation /* or null */ *Computation::set_stop_conditions(
    M2_bool always_stop,
    M2_arrayint degree_limit,
    int basis_element_limit,
    int syzygy_limit,
    int pair_limit,
    int codim_limit,
    int subring_limit,
    M2_bool just_min_gens,
    M2_arrayint length_limit)
{
  stop_.always_stop = always_stop;
  stop_.stop_after_degree = (degree_limit != nullptr && degree_limit->len > 0);
  stop_.degree_limit = degree_limit;
  stop_.basis_element_limit = basis_element_limit;
  stop_.syzygy_limit = syzygy_limit;
  stop_.pair_limit = pair_limit;
  stop_.use_codim_limit = (codim_limit >= 0);
  stop_.codim_limit = codim_limit;
  stop_.subring_limit = subring_limit;
  stop_.just_min_gens = just_min_gens;
  stop_.length_limit = length_limit;

  if (stop_conditions_ok())
    return this;
  else
    return nullptr;
}

Computation::Computation()
{
  computation_status = COMP_NOT_STARTED;

  stop_.always_stop = false;
  stop_.stop_after_degree = false;
  stop_.degree_limit = nullptr;
  stop_.basis_element_limit = 0;
  stop_.syzygy_limit = 0;
  stop_.pair_limit = 0;
  stop_.use_codim_limit = false;
  stop_.codim_limit = 0;
  stop_.subring_limit = 0;
  stop_.just_min_gens = false;
  stop_.length_limit = nullptr;
}

Computation::~Computation() {}

void Computation::text_out(buffer &o) const { o << "-- computation --"; }

void Computation::show() const
{
  printf("No show method available for this computation type\n");
}

enum ComputationStatusCode Computation::set_status(enum ComputationStatusCode c)
{
  switch (computation_status)
    {
      case COMP_OVERFLOWED:
        // if (computation_status == COMP_NEED_RESIZE) break;
        throw(exc::internal_error(
            "attempted to reset status of a computation that overflowed"));
      default:
        return computation_status = c;
    }
}

// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e "
// indent-tabs-mode: nil
// End:
