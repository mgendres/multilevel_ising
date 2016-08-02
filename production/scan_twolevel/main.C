#include <iostream>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream> 

using namespace std;
#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include "math.h"

// Global parameters
#include "lattice.h"
#include "multilevel.h"
#include "mt19937.h"

#define NX 32
#define NY 32

#define THERM 1000
#define CONFS 100000
#define CHKPT 100
#define K_START 0.0
#define K_END 1.0
#define K_INC 0.01
#define K_CRIT 0.44068679350977
#define HITS 5

#define HEAT_BATH true
#define WOLFF_CLUSTER true

double R(double K) {
  if (K<K_CRIT) {
    return 3.2090432488719616/8.0*log(cosh(4.0*K));
  } else {
    return K;
  }
}

int main(void)
{

  cout << setprecision(15);
  rng.init_genrand( time(0) ); // Reseed generator before starting anything

  lattice ising(NX,NY);
  ising.init(DISORDERED); 

  multilevel ml_ising(NX,NY); 
  for (int i=0; i<ml_ising.levels(); ++i) { ml_ising[i].init(ORDERED); }

  ofstream m_file;
  ofstream h_file;
  ofstream c_file;
  m_file.open("data/magnetization.dat");
  h_file.open("data/hamiltonian.dat");
  c_file.open("data/correlator.dat");
  m_file << setprecision(15);
  h_file << setprecision(15);
  c_file << setprecision(15);

  ofstream* m_files;
  ofstream* h_files;
  ofstream* c_files;
  m_files = new ofstream[HITS];
  h_files = new ofstream[HITS];
  c_files = new ofstream[HITS];
  std::stringstream ss;
  for (int i=0; i<HITS; ++i) {
    ss.str(std::string());
    ss << "data/magnetization_" << i <<".dat";
    m_files[i].open(ss.str().c_str());

    ss.str(std::string());
    ss << "data/hamiltonian_" << i <<".dat";
    h_files[i].open(ss.str().c_str());

    ss.str(std::string());
    ss << "data/correlator_" << i <<".dat";
    c_files[i].open(ss.str().c_str());
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
        decimate(ising, ml_ising[1],0,0);
        vector<double> K1;
        vector<double> K0;
        K1.push_back(K);
        K0.push_back(R(K));
        //K1.push_back(0.0);
        //K0.push_back(0.0);
        ml_ising[0].init(ORDERED);
        ml_ising.rg_refine(K0,K1);
        m_file << K << " " << ising.magnetization() << std::endl;
        h_file << K << " " << ising.hamiltonian() << std::endl;;
        c_file << K;
        for (int k=0; k<NX; ++k) { c_file << " " << ising.correlator(0,k); } c_file << std::endl;
        for (int j=0; j<HITS; ++j) {
          m_files[j] << K << " " << ml_ising[0].magnetization() << std::endl;
          h_files[j] << K << " " << ml_ising[0].hamiltonian() << std::endl;
          c_files[j] << K;
          for (int k=0; k<NX; ++k) { c_files[j] << " " << ml_ising[0].correlator(0,k); } c_files[j] << std::endl;
          ml_ising[0].heat_bath(K);
        }
      }
    }
  }

  m_file.close();
  h_file.close();
  c_file.close();
  for (int i=0; i<HITS; ++i) {
    m_files[i].close();
    h_files[i].close();
    c_files[i].close();
  }
  delete[] m_files;
  delete[] h_files;
  delete[] c_files;

  return(0);

}

