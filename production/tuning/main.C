#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include <sstream>


// Global parameters
#include "lattice.h"
#include "mt19937.h"

#define NX 32
#define NY 32
#define CB 1

#define THERM 1000
#define THERM_SUBLAT 100
#define CONFS 100000
#define CHKPT 100
#define K_START 0.0
#define K_END 1.0
#define K_INC 0.1
#define K_CRIT 0.44068679350977

#define HEAT_BATH false
#define WOLFF_CLUSTER true

int main(void)
{

  cout << setprecision(15);
  rng.init_genrand( time(0) ); // Reseed generator before starting anything

  lattice ising(NX,NY); 
  lattice isingC(NX,NY); 
  lattice isingUC(NX/2,NY/2); 
  ising.init(DISORDERED); 

  ofstream* m_files;
  ofstream* h_files;
  m_files = new ofstream[4];
  h_files = new ofstream[4];
  std::stringstream ss;
  for (int i=0; i<4; ++i) {
    ss.str(std::string());
    ss << "data/magnetization_" << i <<".dat";
    m_files[i].open(ss.str().c_str());

    ss.str(std::string());
    ss << "data/hamiltonian_" << i <<".dat";
    h_files[i].open(ss.str().c_str());
  }

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
        copy(ising,isingC);
        for (int j=0; j<THERM_SUBLAT; ++j) {
          if (HEAT_BATH) isingC.heat_bath(K,CB); 
          if (WOLFF_CLUSTER) isingC.wolff_cluster(K,CB); 
        }
        decimate(isingC,isingUC,0,CB);
        // Observables on fine lattice
        m_files[0] << K << " " << ising.magnetization() << std::endl;
        h_files[0] << K << " " << ising.hamiltonian() << std::endl;
        // Observables on coarse lattice, pre-therm
        m_files[1] << K << " " << isingC.magnetization(!CB) << std::endl;
        h_files[1] << K << " " << isingC.hamiltonian(!CB) << std::endl;
        // Observables on coarse lattice, post-therm
        m_files[2] << K << " " << isingC.magnetization(CB) << std::endl;
        h_files[2] << K << " " << isingC.hamiltonian(CB) << std::endl;
        // Observables on ultra-coarse lattice, pre-therm
        m_files[3] << K << " " << isingUC.magnetization() << std::endl;
        h_files[3] << K << " " << isingUC.hamiltonian() << std::endl;
      }
    }
  }

  for (int i=0; i<4; ++i) {
    m_files[i].close();
    h_files[i].close();
  }
  delete[] m_files;
  delete[] h_files;

  return(0);

}

