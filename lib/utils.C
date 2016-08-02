#include "utils.h"

bool power_of_two(unsigned long x)
{
  return (x != 0) && ((x & (x - 1)) == 0);
}



