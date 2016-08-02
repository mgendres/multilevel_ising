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

#define NX 41
#define NY 41

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

  ofstream m_file;
  ofstream h_file;
  ofstream c_file;
  m_file.open("data/magnetization.dat");
  h_file.open("data/hamiltonian.dat");
  c_file.open("data/correlator.dat");
  m_file << setprecision(15);
  h_file << setprecision(15);
  c_file << setprecision(15);

  for (double K=K_START; K<=K_END; K+=K_INC) {
    std::cout << "K : " << K << std::endl;

    for (int i=0; i<THERM; ++i) { // Thermalization
      if (HEAT_BATH) ising.heat_bath(K); 
      if (WOLFF_CLUSTER) ising.wolff_cluster(K); 
    }

    for (int i=0; i<CONFS; ++i) {
      if (HEAT_BATH) ising.heat_bath(K); 
      if (WOLFF_CLUSTER) ising.wolff_cluster(K); 
      if (i%CHKPT==0) {
        m_file << K << " " << ising.magnetization() << std::endl;
        h_file << K << " " << ising.hamiltonian() << std::endl;;
        c_file << K;
        for (int k=0; k<NX; ++k) { c_file << " " << ising.correlator(0,k); } c_file << std::endl;
      }
    }
  }

  m_file.close();
  h_file.close();
  c_file.close();

  return(0);

}

