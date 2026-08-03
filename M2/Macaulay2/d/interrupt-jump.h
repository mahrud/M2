#ifndef M2_INTERRUPT_JUMP_H
#define M2_INTERRUPT_JUMP_H

/* Set this jump and the flag below if the handler should always jump;
   e.g., for interrupting a slow third-party or system library routine. */
#include <setjmp.h>
#ifdef _POSIX_C_SOURCE
# define JMPBUF sigjmp_buf
# define SETJMP(env) sigsetjmp(env, 1)
# define LONGJUMP(env) siglongjmp(env, 1)
#else
# define JMPBUF jmp_buf
# define SETJMP(env) setjmp(env)
# define LONGJUMP(env) longjmp(env, 1)
#endif

struct JumpCell
{
  JMPBUF addr;
  bool is_set;
};

#endif
