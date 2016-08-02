#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include <math.h>

// Global parameters
#include "multilevel.h"
#include "mt19937.h"

#define NX 64
#define NY 64

#define K_CRIT 0.44068679350977
#define P 0.99

double R(double K) {
  return 3.2090432488719616/8.0*log(cosh(4.0*K));
}

int main(void)
{

  cout << setprecision(15);
  rng.init_genrand( time(0) ); // Reseed generator before starting anything

  multilevel ising(NX,NY);
  std::cout << "levels : " << ising.levels();
  std::cout << std::endl;

  if (1) {

    std::cout << std::endl;
    std::cout << "**********" << std::endl;
    std::cout << "* Test 1 *" << std::endl;
    std::cout << "**********" << std::endl;

    vector<double> K0;
    vector<double> K1;

    double K=P*K_CRIT;
    for (int i=0; i<ising.levels()-1; ++i) {
      K1.push_back( K );
      K = R(K);
      K0.push_back( K );
      K = R(K);
    }

    for (int i=0; i<K0.size(); ++i) {
      std::cout << "K1[" << i <<"]" << " :" << K1[i] << std::endl;
      std::cout << "K0[" << i <<"]" << " :" << K0[i] << std::endl;
    }

    ising.rg_refine(K0,K1);

    for (int i=0; i<ising.levels(); ++i) {
      ising[ising.levels()-i-1].print_lattice();
    }

  }

  if (1) {

    std::cout << std::endl;
    std::cout << "**********" << std::endl;
    std::cout << "* Test 2 *" << std::endl;
    std::cout << "**********" << std::endl;

    for (int i=0; i< ising.levels(); ++i) ising[i].init(ORDERED);

    vector<double> K0;
    vector<double> K1;

    double K = P*K_CRIT;
    K1.push_back( K );
    K = R(K);
    K0.push_back( K );
    K = R(K);

    for (int i=0; i<1000; ++i) {
      ising[1].wolff_cluster(K);
    }

    ising.rg_refine(K0,K1);

    for (int i=0; i<ising.levels(); ++i) {
      ising[ising.levels()-i-1].print_lattice();
    }

  }

  if (1) {

    std::cout << std::endl;
    std::cout << "**********" << std::endl;
    std::cout << "* Test 3 *" << std::endl;
    std::cout << "**********" << std::endl;


    for (int i=0; i< ising.levels(); ++i) ising[i].init(ORDERED);

    vector<double> K0;
    vector<double> K1;

    double K = P*K_CRIT;
    K1.push_back( 0.0 );
    K0.push_back( 0.0 );
    K1.push_back( 0.0 );
    K0.push_back( 0.0 );

    ising.rg_refine(K0,K1,4);

    for (int i=0; i<ising.levels(); ++i) {
      ising[ising.levels()-i-1].print_lattice();
    }

  }


  if (1) {
    std::cout << std::endl;
    std::cout << "************" << std::endl;
    std::cout << "* Failures *" << std::endl;
    std::cout << "************" << std::endl;

    vector<double> K0;
    vector<double> K1;

    K0.push_back( 0.0 );
    K0.push_back( 0.0 );
    K1.push_back( 0.0 );

    std::cout << ising.rg_refine(K0, K1) << std::endl;
    std::cout << ising.rg_refine(K1, K0) << std::endl;

    std::cout << ising.rg_refine(K1, K0, K0.size()-1) << std::endl;
    std::cout << ising.rg_refine(K0, K1, K0.size()-1) << std::endl;

    std::cout << ising.rg_refine(K0, K0, ising.levels() ) << std::endl;
    std::cout << ising.rg_refine(K0, K0, K0.size()-1) << std::endl;

  }



  return(0);

}

