#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace std;
#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include <math.h>

// Global parameters
#include "multilevel.h"
#include "mt19937.h"
#include "utils.h"

#define NX 32
#define NY 32

#define THERM 1000
#define CONFS 100000
#define CHKPT 100
#define K_START 0.0
#define K_END 1.0
#define K_INC 0.01
#define K_CRIT 0.44068679350977
#define HITS 10

#define START_LEVEL 1
// Example: 32^2 lattice has 6 levels:
// level volume
// 0     32
// 1     16
// 2      8
// 3      4      <-- START_LEVEL = 3
// 4      2
// 5      1

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

  multilevel ising(NX,NY);
  std::cout << "levels : " << ising.levels() << std::endl;
  for (int i=0; i<ising.levels(); ++i) {
    ising[i].init(ORDERED);
    std::cout << "Level : " << i << std::endl;
    ising[i].print_lattice();
  }

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

  int start_level;
  if (START_LEVEL < 0) {
    start_level = ising.levels() -1;
  } else {
    start_level = START_LEVEL;
  }

  for (double K=K_START; K<=K_END; K+=K_INC) {
    std::cout << "K : " << K << std::endl;

    vector<double> K0;
    vector<double> K1;

    double Kx = K;
    for (int i=0; i<start_level; ++i) {
      K1.push_back( Kx );
      Kx = R(Kx);
      K0.push_back( Kx );
      Kx = R(Kx);
    }

    std::cout << "K1 : ";
    printcoll(K1);
    std::cout << "K0 : "; 
    printcoll(K0);
    std::cout << "Kx : " << Kx << std::endl; 


    for (int i=0; i<THERM; ++i) { // Thermalization
      if (HEAT_BATH) ising[start_level].heat_bath(Kx);
      if (WOLFF_CLUSTER) ising[start_level].wolff_cluster(Kx);
    }

    for (int i=0; i<CONFS; ++i) {

      if (HEAT_BATH) ising[start_level].heat_bath(Kx);
      if (WOLFF_CLUSTER) ising[start_level].wolff_cluster(Kx);
      if (i%CHKPT==0) {
        ising.rg_refine(K0,K1);
        for (int j=0; j<HITS; ++j) {
          m_files[j] << K << " " << ising[0].magnetization() << std::endl;
          h_files[j] << K << " " << ising[0].hamiltonian() << std::endl;
          c_files[j] << K;
          for (int k=0; k<NX; ++k) { c_files[j] << " " << ising[0].correlator(0,k); } c_files[j] << std::endl;
          ising[0].heat_bath(K);
          ising[0].wolff_cluster(K);
        }
      }
    }
  }

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

