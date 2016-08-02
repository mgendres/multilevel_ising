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


#define NX 23
#define NY 13

int main(void)
{

  cout << setprecision(15);
  rng.init_genrand( time(0) ); // Reseed generator before starting anything

  lattice ising(NX,NY); 
  ising.init(DISORDERED); 
  ising.print_lattice();

  std::cout << "magnetization : " << ising.magnetization() << std::endl;
  std::cout << "hamiltonian : " << ising.hamiltonian() << std::endl;
  std::cout << "magnetization (0)  : " << ising.magnetization(0) << std::endl;
  std::cout << "hamiltonian (0) : " << ising.hamiltonian(0) << std::endl;
  std::cout << "magnetization (1) : " << ising.magnetization(1) << std::endl;
  std::cout << "hamiltonian (1) : " << ising.hamiltonian(1) << std::endl;

  return(0);

}

