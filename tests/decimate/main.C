#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
#include <stdlib.h>
#include <stdio.h>
#include <ctime>

// Global parameters
#include "lattice.h"
#include "mt19937.h"


#define NX 4
#define NY 4

int main(void)
{

  cout << setprecision(15);
  rng.init_genrand( time(0) ); // Reseed generator before starting anything

  std::cout << "Fine lattice :" << std::endl;
  lattice isingF(2*NX,2*NY); 
  isingF.init(DISORDERED); 
  isingF.print_lattice();

  std::cout << "Coarse lattice :" << std::endl;
  lattice isingC(NX,NY); 
  isingC.print_lattice();

  std::cout << "Decimate (0,0) :" << std::endl;
  decimate(isingF, isingC,0,0);
  isingC.print_lattice();

  std::cout << "Decimate (1,0) :" << std::endl;
  decimate(isingF, isingC,1,0);
  isingC.print_lattice();

  std::cout << "Decimate (0,1) :" << std::endl;
  decimate(isingF, isingC,0,1);
  isingC.print_lattice();

  std::cout << "Decimate (1,1) :" << std::endl;
  decimate(isingF, isingC,1,1);
  isingC.print_lattice();

  std::cout << "Failures check :" << std::endl;
  std::cout << decimate(isingF, isingC,2,1) << std::endl;
  std::cout << decimate(isingF, isingC,1,2) << std::endl;
  std::cout << decimate(isingF, isingC,-1,1) << std::endl;
  std::cout << decimate(isingF, isingC,1,-1) << std::endl;
  lattice isingX(NX+1,2*NY);
  lattice isingY(2*NX,2*NY+1);
  std::cout << decimate(isingX, isingC,0,0) << std::endl;
  std::cout << decimate(isingY, isingC,0,0) << std::endl;

  return(0);

}

