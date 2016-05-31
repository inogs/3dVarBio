#include "c_exit.h"

#include <stdlib.h>

void C_EXIT(int * exit_code) {
  //printf("c:exit_code=%d\n", *exit_code);
  exit(*exit_code);
}
