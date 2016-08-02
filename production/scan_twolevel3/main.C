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
#include "interpolator.h"

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

#define HEAT_BATH true
#define WOLFF_CLUSTER false


class rg_func {
    interpolator* r;
  public:
    rg_func();
    ~rg_func();
    double operator()(double);
};

rg_func::rg_func()
{
  r = new interpolator("/Users/endres/Documents/MIT/research/bdeop/endres/ising/code_v0.7/Mathematica/migdal_kadanoff.dat");
}

rg_func::~rg_func()
{
  delete r;
}

double rg_func::operator()(double K)
{
  if (K < 0.01) { 
    return 2.718311365715351*K*K - 3.6951050599195727*K*K*K*K + 12002.80348471685*K*K*K*K*K*K;
  }
//  if (K > 18.0) {
//    return 2.0*K;
//  }
  if (K > K_CRIT) {
    return K;
  }


  return r->operator()(K);
}

int main(void)
{

  cout << setprecision(15);
  rng.init_genrand( time(0) ); // Reseed generator before starting anything

//  lattice ising(NX,NY);
//  ising.init(DISORDERED); 

  multilevel ising(NX,NY); 
  for (int i=0; i<ising.levels(); ++i) { ising[i].init(ORDERED); }

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

  rg_func R;
  double RK;
  double RRK;
  for (double K=K_START; K<=K_END; K+=K_INC) {

    RK = R(K);
    RRK = R(R(K));

    std::cout << "K : " << K << std::endl;
    std::cout << "RK : " << RK << std::endl;
    std::cout << "RRK : " << RRK << std::endl;

    for (int i=0; i<THERM; ++i) { // Thermalization
      if (HEAT_BATH) ising[1].heat_bath(RRK); 
      if (WOLFF_CLUSTER) ising[1].wolff_cluster(RRK); 
    }

    for (int i=0; i<CONFS; ++i) {
      if (HEAT_BATH) ising[1].heat_bath(RRK); 
      if (WOLFF_CLUSTER) ising[1].wolff_cluster(RRK); 
      if (i%CHKPT==0) {
        vector<double> K1;
        vector<double> K0;
        K1.push_back(K);
        K0.push_back(RK);
        //K1.push_back(0.0);
        //K0.push_back(0.0);
        ising[0].init(ORDERED);
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

