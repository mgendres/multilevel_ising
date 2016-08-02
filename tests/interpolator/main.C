#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
#include <stdlib.h>
#include <stdio.h>
#include <ctime>

// Global parameters
#include "interpolator.h"


class rg_func {
    interpolator* r;
  public:
    rg_func();
    ~rg_func();
    double operator()(double);
};

rg_func::rg_func()
{
  r = new interpolator("migdal_kadanoff.dat");
}

rg_func::~rg_func()
{
  delete r;
}

double rg_func::operator()(double K)
{
  if (K < 0.01) { 
    return 2.718311365715351*K*K - 3.6951050599195727*K*K*K*K + 12002.80348471685*K*K*K*K*K*K;
  }
  if (K > 18.0) {
    return 2.0*K;
  }
  return r->operator()(K);
}


int main(void)
{

  cout << setprecision(15);
  rg_func R;

  if (0) {
    for (double K=0.0; K<10.0; K+=0.01) {
      std::cout << K << " " <<  R(K) << std::endl;
    }
  }

  if (1) {
    double Kx;
    for (double K=0.0; K<2.0; K+=0.01) {
      Kx = K;
      for (int j=0; j<7; ++j) {
        std::cout << Kx << " ";
        Kx = R(Kx);
      }
      std::cout << std::endl;
    }
  }

  return(0);

}

