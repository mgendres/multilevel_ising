#include "multilevel.h"
#include <iostream>
#include "mt19937.h"
#include <math.h>
#include "utils.h"


multilevel::multilevel(int arg_nx, int arg_ny)
{

  int nx = arg_nx;
  int ny = arg_ny;

  lattice* lat_p;
  while (power_of_two(nx)&&power_of_two(ny)) {
      lat_p = new lattice(nx,ny);
      lat.push_back( lat_p );
    nx /= 2;
    ny /= 2;
  }


}

multilevel::~multilevel() {

  for(int i=0; i < lat.size(); ++i){
    delete lat[i];
  }

}

lattice & multilevel::operator[](int i)
{
  return *lat[i];
}

int multilevel::levels()
{
  return lat.size();
}


bool multilevel::rg_refine(const std::vector<double>K0, const std::vector<double>K1)
{

  lat[lat.size()-1]->init(DISORDERED);
  int levels = K0.size();
  if ( !rg_refine(K0, K1, levels) ) return false;
  return true;
}

bool multilevel::rg_refine(const std::vector<double>K0, const std::vector<double>K1, int start_level)
// Levels labeled: 0, 1, ... levels-1
// start_level should be one of these values
{

  int levels = lat.size();

  if ( K0.size() != K1.size() ) {
    std::cout << "ERROR : K0 and K1 must be of equal length" << std::endl;
    return false;
  }

  if ( start_level >= levels ) {
    std::cout << "ERROR : start_level must be less than levels" << std::endl;
    return false;
  }

  int refinements = K0.size();
  if ( refinements > start_level ) {
    std::cout << "ERROR : number of refinements must be of length start_level or less" << std::endl;
    return false;
  }

  int j;
  for (int i=start_level; i>start_level-refinements; --i) {
    j = refinements+i-start_level-1;
//    std::cout << "level : " << i << std::endl;
//    std::cout << "refinement to level " << i-1 << " using :" <<  std::endl;
//    std::cout << "K0 : " << K0[j] << std::endl;
//    std::cout << "K1 : " << K1[j] << std::endl;
    if ( !refine(*lat[i], *lat[i-1], K0[j], K1[j]) ) return false;
  }

  return true;

}

