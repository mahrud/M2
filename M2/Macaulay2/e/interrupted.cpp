#include "tasks/supervisor.hpp"
#include "tasks/supervisorinterface.h"
#define interrupted() \
  test_Field(THREADLOCAL(interrupts_interruptedFlag, struct atomic_field))

bool system_interrupted() { return interrupted(); }
