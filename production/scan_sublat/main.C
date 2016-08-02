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

// Must be an even number, otherwise pbcs connect all sites!!
#define NX 58
#define NY 58
#define CB 0

#define THERM 1000
#define CONFS 100000
#define CHKPT 100
#define K_START 0.0
#define K_END 1.0
#define K_INC 0.01
#define K_CRIT 0.44068679350977

#define HEAT_BATH true
#define WOLFF_CLUSTER true

int main(void)
{

  cout << setprecision(15);
  rng.init_genrand( time(0) ); // Reseed generator before starting anything

  lattice ising(NX,NY); 
  ising.init(ORDERED); 

  ofstream m0_file;
  ofstream h0_file;
  ofstream c0_file;
  m0_file.open("data/magnetization_0.dat");
  h0_file.open("data/hamiltonian_0.dat");
  c0_file.open("data/correlator_0.dat");
  m0_file << setprecision(15);
  h0_file << setprecision(15);
  c0_file << setprecision(15);

  ofstream m1_file;
  ofstream h1_file;
  ofstream c1_file;
  m1_file.open("data/magnetization_1.dat");
  h1_file.open("data/hamiltonian_1.dat");
  c1_file.open("data/correlator_1.dat");
  m1_file << setprecision(15);
  h1_file << setprecision(15);
  c1_file << setprecision(15);

  for (double K=K_START; K<=K_END; K+=K_INC) {
    std::cout << "K : " << K << std::endl;

    for (int i=0; i<THERM; ++i) { // Thermalization
      if (HEAT_BATH) ising.heat_bath(K,CB); 
      if (WOLFF_CLUSTER) ising.wolff_cluster(K,CB); 
    }

    for (int i=0; i<CONFS; ++i) {
      if (HEAT_BATH) ising.heat_bath(K,CB); 
      if (WOLFF_CLUSTER) ising.wolff_cluster(K,CB); 
      if (i%CHKPT==0) {

        m0_file << K << " " << ising.magnetization(0) << std::endl;
        h0_file << K << " " << ising.hamiltonian(0) << std::endl;;
        c0_file << K;
        for (int k=0; k<NX; ++k) { c0_file << " " << 2.0*ising.correlator(k,k) - 1.0; } c0_file << std::endl;

        m1_file << K << " " << ising.magnetization(1) << std::endl;
        h1_file << K << " " << ising.hamiltonian(1) << std::endl;;
        c1_file << K;
        for (int k=0; k<NX; ++k) { c1_file << " " << 2.0*ising.correlator(k,k) - 1.0; } c1_file << std::endl;

      }
    }
  }

  m0_file.close();
  h0_file.close();
  c0_file.close();

  m1_file.close();
  h1_file.close();
  c1_file.close();

  return(0);

}

