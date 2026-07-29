// (c) 1997 Michael E. Stillman

// This stuff would not be necessary if exceptions worked well,
// unfortunately, there are internal errors so that they fail
// at least in the presence of template code (even in version 2.7.2).
// So, we use the following simple, if unpleasant, scheme.
//
// Upon finding an error, call ERROR any number of times.  Each time,
// The message is added to the end (but, it won't write over the end of
// the buffer).
//
// Then, it is unfortunately necessary to check whether an error has
// occurred, using 'error'.
//
// Finally, the error flag is cleared upon giving control back to the front end.

#include "error.h"
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

// The flag and the message are per-thread: ERROR is called by the engine on
// whichever thread made the call, and error()/error_message() are read back by
// the D layer on that same thread (see EngineError in d/engine.dd).  With a
// single shared buffer, two threads erroring at once clobber each other and one
// of them ends up reporting "unknown engine error" instead of the real message.
#if defined(_MSC_VER)
#define M2_THREADLOCAL __declspec(thread)
#else
#define M2_THREADLOCAL __thread
#endif

#define MAXERROR 200
static M2_THREADLOCAL int iserror = 0;
static M2_THREADLOCAL char errmsg[MAXERROR] = {'\0'};

void ERROR(const char *s, ...)
{
  va_list ap;
  if (iserror) fprintf(stderr, "--error message bumped: %s\n", errmsg);
  iserror = 1;
  va_start(ap, s);
  vsprintf(errmsg, s, ap);
  va_end(ap);
}

void INTERNAL_ERROR(const char *s, ...)
{
  char buf[MAXERROR];
  buf[0] = 0;
  va_list ap;
  va_start(ap, s);
  vsprintf(buf, s, ap);
  va_end(ap);
  fprintf(stderr, "--internal error: %s\n", buf);
  abort(); /* we should exit after an internal error, to avoid crashing */
}

int error() { return iserror; }
const char *error_message()
{
  if (iserror == 0) return "";
  iserror = 0;
  return errmsg;
}

/*
// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e "
// indent-tabs-mode: nil
// End:
*/
