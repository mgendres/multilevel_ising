#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include <math.h>

// Global parameters
#include "utils.h"


int main(void)
{

  cout << setprecision(15);

  int count = 0;
  for (int l=0; l<10000000; ++l) {
    if (power_of_two(l)) {
      std::cout << count << " : " <<  l << std::endl;
      count++;
    }
  }

  return(0);

}

