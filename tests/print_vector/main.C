#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include <math.h>
#include "utils.h"



int main(void)
{

  cout << setprecision(15);

  std::vector<int> v;

  v.push_back(42);
  v.push_back(13);
  v.push_back(7);
  std::cout << "Printing a vector<int> :" << std::endl;
  printcoll(v);

  std::vector<double> w;

  w.push_back(4.2);
  w.push_back(1.3);
  w.push_back(.7);
  std::cout << "Printing a vector<double> :" << std::endl;
  printcoll(w);



  return(0);

}

