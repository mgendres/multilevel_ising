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
#define NY 41

int main(void)
{

  cout << setprecision(15);
  rng.init_genrand( time(0) ); // Reseed generator before starting anything

  lattice ising(NX,NY); 
  ising.init(ORDERED); 
  cout << "Magnetization : " << ising.magnetization() << endl;
  cout << "Hamiltonian: " << ising.hamiltonian() << endl;
  ising.print_lattice();

  ising.init(DISORDERED); 
  cout << "Magnetization : " << ising.magnetization() << endl;
  cout << "Hamiltonian: " << ising.hamiltonian() << endl;
  ising.print_lattice();

  return(0);

}

