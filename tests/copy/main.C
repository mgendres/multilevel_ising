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


#define NX 7
#define NY 13

int main(void)
{

  cout << setprecision(15);
  rng.init_genrand( time(0) ); // Reseed generator before starting anything

  lattice ising0(NX,NY); 
  ising0.init(DISORDERED); 
  std::cout << "Lattice 0 :" << std::endl;
  ising0.print_lattice();

  lattice ising1(NX,NY); 
  ising1.init(ORDERED); 
  std::cout << "Lattice 1 :" << std::endl;
  ising1.print_lattice();

  lattice ising2(NX,NY+1); 
  ising2.init(ORDERED); 
  std::cout << "Lattice 2 :" << std::endl;
  ising2.print_lattice();

  lattice ising3(NX+1,NY); 
  ising3.init(ORDERED); 
  std::cout << "Lattice 3 :" << std::endl;
  ising3.print_lattice();

  if ( copy(ising0,ising1) ) {
    std::cout << "Copy lattice 0 to lattice 1 : SUCCESSFUL" << std::endl;
  } else {
    std::cout << "Copy lattice 0 to lattice 1 : FAILED" << std::endl;
  }

  if ( copy(ising0,ising2) ) {
    std::cout << "Copy lattice 0 to lattice 1 : SUCCESSFUL" << std::endl;
  } else {
    std::cout << "Copy lattice 0 to lattice 1 : FAILED" << std::endl;
  }

  if ( copy(ising0,ising3) ) {
    std::cout << "Copy lattice 0 to lattice 1 : SUCCESSFUL" << std::endl;
  } else {
    std::cout << "Copy lattice 0 to lattice 1 : FAILED" << std::endl;
  }

  std::cout << "Lattice 0 :" << std::endl;
  ising0.print_lattice();

  std::cout << "Lattice 1 :" << std::endl;
  ising1.print_lattice();

  std::cout << "Lattice 2 :" << std::endl;
  ising2.print_lattice();

  std::cout << "Lattice 3 :" << std::endl;
  ising3.print_lattice();

  return(0);

}

