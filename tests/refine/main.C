#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include <math.h>

// Global parameters
#include "lattice.h"
#include "mt19937.h"

#define NX 7
#define NY 13
#define K_CRIT 0.44068679350977

double R(double K) {
  return 3.15/8.0*log(cosh(4.0*K));
}

int main(void)
{

  cout << setprecision(15);
  rng.init_genrand( time(0) ); // Reseed generator before starting anything

  std::cout << "Coarse lattice :" << std::endl;
  lattice isingC(NX,NY); 
  isingC.init(DISORDERED); 
  isingC.print_lattice();

  std::cout << "Fine lattice :" << std::endl;
  lattice isingF(2*NX,2*NY); 
  isingF.init(ORDERED); 
  isingF.print_lattice();

  std::cout << "Refine (K0 = 0, K1 = 0) :" << std::endl;
  refine(isingC, isingF,0,0);
  isingF.print_lattice();

  std::cout << "Refine (K0 = Kc, K1 = Kc) :" << std::endl;
  refine(isingC, isingF, K_CRIT, K_CRIT);
  isingF.print_lattice();

  std::cout << "Refine (K0 = 99.9, K1 = 99.9) :" << std::endl;
  refine(isingC, isingF, 99.9, 99.9);
  isingF.print_lattice();

  std::cout << "Failures check :" << std::endl;
  lattice isingX(NX+1,2*NY);
  lattice isingY(2*NX,2*NY+1);
  std::cout << decimate(isingX, isingC,0,0) << std::endl;
  std::cout << decimate(isingY, isingC,0,0) << std::endl;

  return(0);

}

