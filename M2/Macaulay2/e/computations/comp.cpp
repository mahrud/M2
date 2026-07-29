// Copyright 2004 Michael E. Stillman

#include "computations/comp.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <chrono>

#include "buffer.hpp"
#include "error.h"
#include "exceptions.hpp"
#include "finalize.hpp"
#include "interrupted.hpp"

namespace {
std::atomic<std::size_t> nextThreadTag{1};
}

std::size_t Computation::currentThreadTag()
{
  static thread_local std::size_t tag =
      nextThreadTag.fetch_add(1, std::memory_order_relaxed);
  return tag;
}

// Record that the current thread now holds C->mMutex.  Only the outermost
// acquisition sets mOwner; nested ones just deepen the count.  Both fields are
// only ever touched while holding the mutex.
static void noteAcquired(std::atomic<std::size_t> &owner, int &depth)
{
  if (depth++ == 0)
    owner.store(Computation::currentThreadTag(), std::memory_order_relaxed);
}

static void noteReleased(std::atomic<std::size_t> &owner, int &depth)
{
  if (--depth == 0) owner.store(0, std::memory_order_relaxed);
}

bool Computation::failFastOnThreadConflict()
{
  static const bool failFast = [] {
    const char *s = getenv("M2_ENGINE_THREAD_CONFLICT");
    return s != nullptr && strcmp(s, "error") == 0;
  }();
  return failFast;
}

Computation::ThreadLock::ThreadLock(Computation *C)
    : mComputation(C), mOtherThread(0), mOk(false), mInterrupted(false)
{
  // The mutex is recursive, so try_lock always succeeds for the thread that
  // already holds this computation; only a genuinely concurrent thread waits.
  if (C->mMutex.try_lock())
    {
      noteAcquired(C->mOwner, C->mOwnerDepth);
      mOk = true;
      return;
    }

  mOtherThread = C->mOwner.load(std::memory_order_relaxed);
  if (failFastOnThreadConflict()) return;

  // Wait for the other thread, but in slices, so that a thread parked behind a
  // long Groebner basis computation still answers ^C instead of hanging.
  while (!C->mMutex.try_lock_for(std::chrono::milliseconds(100)))
    {
      if (system_interrupted())
        {
          mInterrupted = true;
          return;
        }
      mOtherThread = C->mOwner.load(std::memory_order_relaxed);
    }
  noteAcquired(C->mOwner, C->mOwnerDepth);
  mOk = true;
}

Computation::ThreadLock::~ThreadLock()
{
  if (!mOk) return;
  noteReleased(mComputation->mOwner, mComputation->mOwnerDepth);
  mComputation->mMutex.unlock();
}

void Computation::ThreadLock::reportConflict() const
{
  if (mInterrupted)
    ERROR("interrupted while waiting for an engine computation held by thread %lu",
          static_cast<unsigned long>(mOtherThread));
  else
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
