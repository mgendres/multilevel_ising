#include <sstream>
#include <string>
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

#define NX 64
#define NY 64

#define THERM 1000
#define CONFS 100000
#define CHKPT 100
//#define K 0.44068679350977
//#define K 0.4
#define K 0.5

int main(void)
{

  cout << setprecision(15);
  rng.init_genrand( time(0) ); // Reseed generator before starting anything

  lattice ising(NX,NY); 
  ising.init(ORDERED); 

  ofstream file;
  std::stringstream ss;
  ss << setprecision(15);
  ss << "data/cfgs_" << K << ".dat";
  std::string s = ss.str();
  file.open(s);
  file << setprecision(15);
  std::cout << "K : " << K << std::endl;

  for (int i=0; i<THERM; ++i) { // Thermalization
    ising.wolff_cluster(K); 
  }

  for (int i=0; i<CONFS; ++i) {
    ising.wolff_cluster(K); 
    if (i%CHKPT==0) {
      std::cout << "Config: " << i << endl;
      // Write to disK
      //ising.print_lattice();
      for (int x=0; x<NX; ++x) {
        for (int y=0; y<NX; ++y) {
          file << ising.elem(x,y) << " ";
        }
      }
      file << endl;
    }
  }

  file.close();

  return(0);

}

