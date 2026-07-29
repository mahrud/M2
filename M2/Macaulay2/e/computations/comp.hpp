// Copyright 2004 Michael E. Stillman

#ifndef M2_COMPUTATIONS_COMP_HPP_
#define M2_COMPUTATIONS_COMP_HPP_

#include <atomic>
#include <cstddef>

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
  // Rather than race, the interface routines claim the computation with a
  // ThreadGuard and report a top-level error when a second thread shows up.
  //
  // mOwner is the tag of the thread currently inside this computation, or 0
  // when nobody holds it; mOwnerDepth is touched only by the owning thread.
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

  // RAII claim of exclusive use of a computation by the current thread.
  // Re-entrant for the thread that already holds the claim, so nested engine
  // calls on one thread are fine.  If another thread holds it, ok() is false
  // and the caller must bail out instead of proceeding into the computation.
  class ThreadGuard
  {
   private:
    Computation *mComputation;
    std::size_t mOtherThread;
    bool mOk;

   public:
    explicit ThreadGuard(Computation *C);
    ~ThreadGuard();

    ThreadGuard(const ThreadGuard &) = delete;
    ThreadGuard &operator=(const ThreadGuard &) = delete;

    bool ok() const { return mOk; }
    // tag of the thread that is already inside the computation; only
    // meaningful when ok() is false
    std::size_t otherThread() const { return mOtherThread; }
    // issues the standard top-level error message; only call when !ok()
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
