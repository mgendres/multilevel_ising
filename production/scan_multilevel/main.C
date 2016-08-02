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

#define NX 128
#define NY 128

#define CONFS 100
#define K_START 0.0
#define K_END 1.0
#define K_INC 0.01
#define K_CRIT 0.44068679350977
#define HITS 5

double R(double K) {
  return 3.2090432488719616/8.0*log(cosh(4.0*K));
}

int main(void)
{

  cout << setprecision(15);
  rng.init_genrand( time(0) ); // Reseed generator before starting anything

  multilevel ising(NX,NY);

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

    vector<double> K0;
    vector<double> K1;

    double Kx = K;
    for (int i=0; i<ising.levels()-1; ++i) {
      K1.push_back( Kx );
      Kx = R(Kx);
      K0.push_back( Kx );
      Kx = R(Kx);
    } 

//  for (int i=0; i<ising.levels(); ++i) {
//    std::cout << "K0[" << i <<"]" << " :" << K0[i] << std::endl;
//    std::cout << "K1[" << i <<"]" << " :" << K1[i] << std::endl;
//  }

    for (int i=0; i<CONFS; ++i) {
      ising[ising.levels()-1].init(DISORDERED);
      ising.rg_refine(K0,K1);

      for (int j=0; j<HITS; ++j) {
        m_files[j] << K << " " << ising[0].magnetization() << std::endl;
        h_files[j] << K << " " << ising[0].hamiltonian() << std::endl;
        c_files[j] << K;
        for (int k=0; k<NX; ++k) { c_files[j] << " " << ising[0].correlator(0,k); } c_files[j] << std::endl;
        ising[0].heat_bath(K);
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

