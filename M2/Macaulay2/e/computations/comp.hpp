// Copyright 2004 Michael E. Stillman

#ifndef M2_COMPUTATIONS_COMP_HPP_
#define M2_COMPUTATIONS_COMP_HPP_

#include <atomic>
#include <cstddef>
#include <mutex>

#include "interface/computation.h"
#include "hash.hpp"

class GBComputation;
class ResolutionComputation;

class buffer;

/**
    @ingroup computations
*/

class Computation : public MutableEngineObject
{
 private:
  enum ComputationStatusCode computation_status;

  // Engine computations (gbA in particular) carry mutable internal state and
  // are not thread safe.  Two interpreter threads driving the same computation
  // object -- which happens when parallel tasks force the same cached Groebner
  // basis or resolution -- corrupt that state, and the symptom is a SIGSEGV or
  // "pure virtual method called" deep inside the engine, far from the cause.
  //
  // Every interface routine that enters a computation therefore takes mMutex
  // first, via ThreadLock.  Concurrent callers serialize rather than race: the
  // second thread waits, and by the time it gets in the work is usually already
  // done, so sharing a computation across tasks is merely slow, not fatal.
  //
  // The mutex is recursive so nested engine calls on one thread do not
  // self-deadlock, and timed so a thread waiting on a long computation can
  // still be interrupted.  mOwner is the tag of the thread holding it (0 when
  // free), used only to name that thread when reporting a conflict.
  std::recursive_timed_mutex mMutex;
  std::atomic<std::size_t> mOwner{0};
  int mOwnerDepth{0};

 protected:
  StopConditions stop_;

  Computation();

  enum ComputationStatusCode set_status(enum ComputationStatusCode);

  virtual bool stop_conditions_ok() = 0;
  // If the stop conditions in stop_ are inappropriate,
  // return false, and use ERROR(...) to provide an error message.

  virtual ~Computation();

 public:
  // A unique small integer identifying the calling thread, for error messages.
  static std::size_t currentThreadTag();

  // True when the engine has been asked to report thread conflicts instead of
  // serializing through them -- set M2_ENGINE_THREAD_CONFLICT=error in the
  // environment.  Useful for finding the M2 code that shares a computation
  // between tasks, which is otherwise invisible once the mutex hides it.
  static bool failFastOnThreadConflict();

  // RAII claim of exclusive use of a computation by the current thread.
  //
  // Normally this blocks until whichever thread is inside the computation is
  // done, so ok() is true and callers can ignore it.  It returns false in two
  // cases, and then the caller must bail out without entering the computation:
  //
  //   - fail-fast mode (see above), where a concurrent thread is a reportable
  //     error rather than something to wait for;
  //   - the wait was interrupted, so we stop waiting and let the interpreter
  //     unwind.
  //
  // Either way reportConflict() issues the appropriate top-level message.
  // Re-entrant: the thread already holding the computation is let straight
  // back in, so nested engine calls on one thread cannot self-deadlock.
  class ThreadLock
  {
   private:
    Computation *mComputation;
    std::size_t mOtherThread;
    bool mOk;
    bool mInterrupted;

   public:
    explicit ThreadLock(Computation *C);
    ~ThreadLock();

    ThreadLock(const ThreadLock &) = delete;
    ThreadLock &operator=(const ThreadLock &) = delete;

    bool ok() const { return mOk; }
    // tag of the thread that is inside the computation; only meaningful
    // when ok() is false
    std::size_t otherThread() const { return mOtherThread; }
    // issues the appropriate top-level error; only call when !ok()
    void reportConflict() const;
  };

  Computation /* or null */ *set_stop_conditions(M2_bool always_stop,
                                                 M2_arrayint degree_limit,
                                                 int basis_element_limit,
                                                 int syzygy_limit,
                                                 int pair_limit,
                                                 int codim_limit,
                                                 int subring_limit,
                                                 M2_bool just_min_gens,
                                                 M2_arrayint length_limit);
  // returns NULL if there is a general problem with one of the stop
  // conditions.

  enum ComputationStatusCode status() const { return computation_status; }
  virtual int complete_thru_degree() const = 0;
  // This is computation specific information.  However, for homogeneous
  // GB's, the GB coincides with the actual GB in degrees <= the returned value.
  // For resolutions of homogeneous modules, the resolution
  // coincides with the actual one in (slanted) degrees <= the returned value.

  virtual void start_computation() = 0;
  // Do the computation as specified by the stop conditions.
  // This routine should set the status of the computation.

  virtual GBComputation *cast_to_GBComputation() { return nullptr; }
  virtual ResolutionComputation *cast_to_ResolutionComputation() { return nullptr; }
  virtual void text_out(buffer &o) const;

  virtual void show() const;  // debug display of some computations
};

#endif

// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e "
// indent-tabs-mode: nil
// End:
