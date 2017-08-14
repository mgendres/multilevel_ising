#include "lattice.h"
#include <iostream>
#include "mt19937.h"
#include <math.h>


lattice::lattice(int arg_nx, int arg_ny)
: nx(arg_nx),
  ny(arg_ny)
{

  lat = new int* [nx];
  for (int i=0; i<nx; ++i) {
    lat[i] = new int [ny];
  }

  cluster = new bool* [nx];
  for (int i=0; i<nx; ++i) {
    cluster[i] = new bool [ny];
  }

  init(ORDERED); // default

}

lattice::~lattice() {
  for (int i=0; i<nx; ++i) {
    delete [] lat[i];
    delete [] cluster[i];
  }
  delete [] lat;
  delete [] cluster;
}


void lattice::init(enum lat_state state) {

  if (state==ORDERED)
  for (int i=0; i<nx; ++i)
  for (int j=0; j<ny; ++j) {
    lat[i][j] = 1;
  }

  if (state==DISORDERED)
  for (int i=0; i<nx; ++i)
  for (int j=0; j<ny; ++j) {
    lat[i][j] = 2*rng.genrand_int53()-1;
  }

}


// K == beta*J
void lattice::heat_bath(double K) {
  double p[9];
  for (int s=0; s<9; s+=2) {
     p[s] = 1.0;
     p[s] += tanh( K*(s-4) ); 
     p[s] *= 0.5;
  }

  heat_bath_eo(p, 0);
  heat_bath_eo(p, 1);

}

void lattice::heat_bath(double K, double KK) {
  double* p[9];
  for (int s=0; s<9; s+=2) {
    p[s] = new double[9];
    for (int ss=0; ss<9; ss+=2) {
        p[s][ss] = 1.0;
        p[s][ss] += tanh( K*(s-4) + KK*(ss-4) ); 
        p[s][ss] *= 0.5;
    }
  }

  heat_bath_eo((double**)p, 0);
  heat_bath_eo((double**)p, 1);

  for (int s=0; s<9; s+=2) { delete [] p[s]; }

}


// If cb = 0, update even sites
// If cb = 1, update odd sites
void lattice::heat_bath_eo(double *p, int cb) {

  int s;
  int iPrev, jPrev, iNext, jNext;
  for (int i=0; i<nx; ++i)
  for (int j=0; j<ny; ++j)
  if( (i+j)%2==cb ) {
    iPrev = i-1;
    jPrev = j-1;
    iNext = i+1;
    jNext = j+1;
 
    s = 0;
    s += lat[iNext==nx ? 0 : iNext][j];
    s += lat[i==0 ? nx-1 : iPrev][j];
    s += lat[i][jNext==ny ? 0 : jNext];
    s += lat[i][j==0 ? ny-1 : jPrev];
    if ( p[s+4] > rng.genrand_res53() ) {
      lat[i][j] = 1;
    } else {
      lat[i][j] = -1;
    }
  }

}

// If cb = 0, update even sites
// If cb = 1, update odd sites
void lattice::heat_bath_eo(double **p, int cb) {

  int s, ss;
  int iPrev, jPrev, iNext, jNext;
  for (int i=0; i<nx; ++i)
  for (int j=0; j<ny; ++j)
  if( (i+j)%2==cb ) {
    iPrev = i-1;
    jPrev = j-1;
    iNext = i+1;
    jNext = j+1;

    s = 0;
    s += lat[iNext==nx ? 0 : iNext][j];
    s += lat[i==0 ? nx-1 : iPrev][j];
    s += lat[i][jNext==ny ? 0 : jNext];
    s += lat[i][j==0 ? ny-1 : jPrev];

    ss = 0;
    ss += lat[iNext==nx ? 0 : iNext][jNext==ny ? 0 : jNext];
    ss += lat[i==0 ? nx-1 : iPrev][jNext==ny ? 0 : jNext];
    ss += lat[iNext==nx ? 0 : iNext][j==0 ? ny-1 : jPrev];
    ss += lat[i==0 ? nx-1 : iPrev][j==0 ? ny-1 : jPrev];
    
    if ( p[s+4][ss+4] > rng.genrand_res53() ) {
      lat[i][j] = 1;
    } else {
      lat[i][j] = -1;
    }
  }

}



// K == beta*J
void lattice::heat_bath(double K, int cb) {

  double p[9];
  for (int s=0; s<9; s+=2) {
     p[s] = 1.0;
     p[s] += tanh( K*(s-4) ); 
     p[s] *= 0.5;
  }

  if (!cb) {
    heat_bath_eeoo(p, 0);
    heat_bath_eeoo(p, 1);
  } else {
    heat_bath_eooe(p, 0);
    heat_bath_eooe(p, 1);
  }
}

void lattice::heat_bath(double K, double KK, int cb) {

  double* p[9];
  for (int s=0; s<9; s+=2) {
    p[s] = new double[9];
    for (int ss=0; ss<9; ss+=2) {
        p[s][ss] = 1.0;
        p[s][ss] += tanh( K*(s-4) + KK*(ss-4) ); 
        p[s][ss] *= 0.5;
    }
  }

  if (!cb) {
    heat_bath_eeoo((double**)p, 0);
    heat_bath_eeoo((double**)p, 1);
  } else {
    heat_bath_eooe((double**)p, 0);
    heat_bath_eooe((double**)p, 1);
  }

  for (int s=0; s<9; s+=2) { delete [] p[s]; }

}



// If cb = 0, update (even,even) sites
// If cb = 1, update (odd,odd) sites
// Updates are performed bases on nearest diagonal neighbors
void lattice::heat_bath_eeoo(double *p, int cb) {

  int s;
  int iPrev, jPrev, iNext, jNext;
  for (int i=0; i<nx; ++i)
  for (int j=0; j<ny; ++j)
  if( (i%2==cb)&&(j%2==cb) ) {
    iPrev = (i==0 ? nx-1 : i-1);
    jPrev = (j==0 ? ny-1 : j-1);
    iNext = (i==nx-1 ? 0 : i+1);
    jNext = (j==ny-1 ? 0 : j+1);
 
    s = 0;
    s += lat[iPrev][jPrev];
    s += lat[iNext][jPrev];
    s += lat[iPrev][jNext];
    s += lat[iNext][jNext];
    if ( p[s+4] > rng.genrand_res53() ) {
      lat[i][j] = 1;
    } else {
      lat[i][j] = -1;
    }
  }
}

// If cb = 0, update (even,odd) sites
// If cb = 1, update (odd,even) sites
// Updates are performed bases on nearest diagonal neighbors
void lattice::heat_bath_eooe(double *p, int cb) {

  int s;
  int iPrev, jPrev, iNext, jNext;
  for (int i=0; i<nx; ++i)
  for (int j=0; j<ny; ++j)
  if( (i%2==cb)&&(j%2==!cb) ) {
    iPrev = (i==0 ? nx-1 : i-1);
    jPrev = (j==0 ? ny-1 : j-1);
    iNext = (i==nx-1 ? 0 : i+1);
    jNext = (j==ny-1 ? 0 : j+1);
 
    s = 0;
    s += lat[iPrev][jPrev];
    s += lat[iNext][jPrev];
    s += lat[iPrev][jNext];
    s += lat[iNext][jNext];
    if ( p[s+4] > rng.genrand_res53() ) {
      lat[i][j] = 1;
    } else {
      lat[i][j] = -1;
    }
  }
}

void lattice::heat_bath_eeoo(double **p, int cb) {

  int s, ss;
  int iPrev, jPrev, iNext, jNext;
  int iPrevPrev, jPrevPrev, iNextNext, jNextNext;
  for (int i=0; i<nx; ++i)
  for (int j=0; j<ny; ++j)
  if( (i%2==cb)&&(j%2==cb) ) {
    iPrev = (i==0 ? nx-1 : i-1);
    jPrev = (j==0 ? ny-1 : j-1);
    iNext = (i==nx-1 ? 0 : i+1);
    jNext = (j==ny-1 ? 0 : j+1);

    iPrevPrev = (i<2 ? i-2+nx : i-2);
    jPrevPrev = (j<2 ? j-2+ny : j-2);
    iNextNext = (i+2>nx-1 ? i+2-nx : i+2);
    jNextNext = (j+2>ny-1 ? j+2-ny : j+2);

    s = 0;
    s += lat[iPrev][jPrev];
    s += lat[iNext][jPrev];
    s += lat[iPrev][jNext];
    s += lat[iNext][jNext];

    ss = 0;
    ss += lat[iPrevPrev][j];
    ss += lat[iNextNext][j];
    ss += lat[i][jPrevPrev];
    ss += lat[i][jNextNext];

    if ( p[s+4][ss+4] > rng.genrand_res53() ) {
      lat[i][j] = 1;
    } else {
      lat[i][j] = -1;
    }
  }
}

void lattice::heat_bath_eooe(double **p, int cb) {

  int s, ss;
  int iPrev, jPrev, iNext, jNext;
  int iPrevPrev, jPrevPrev, iNextNext, jNextNext;
  for (int i=0; i<nx; ++i)
  for (int j=0; j<ny; ++j)
  if( (i%2==cb)&&(j%2==!cb) ) {
    iPrev = (i==0 ? nx-1 : i-1);
    jPrev = (j==0 ? ny-1 : j-1);
    iNext = (i==nx-1 ? 0 : i+1);
    jNext = (j==ny-1 ? 0 : j+1);

    iPrevPrev = (i<2 ? i-2+nx : i-2);
    jPrevPrev = (j<2 ? j-2+ny : j-2);
    iNextNext = (i+2>nx-1 ? i+2-nx : i+2);
    jNextNext = (j+2>ny-1 ? j+2-ny : j+2);

    s = 0;
    s += lat[iPrev][jPrev];
    s += lat[iNext][jPrev];
    s += lat[iPrev][jNext];
    s += lat[iNext][jNext];

    ss = 0;
    ss += lat[iPrevPrev][j];
    ss += lat[iNextNext][j];
    ss += lat[i][jPrevPrev];
    ss += lat[i][jNextNext];

    if ( p[s+4][ss+4] > rng.genrand_res53() ) {
      lat[i][j] = 1;
    } else {
      lat[i][j] = -1;
    }
  }
}

// K == beta*J
void lattice::wolff_cluster(double K) {

  add_prob = 1.0 - exp(-2.0*K);
  for (int i=0; i<nx; ++i)
  for (int j=0; j<ny; ++j) {
    cluster[i][j] = false;
  }

  int i = int(rng.genrand_res53() * nx);
  int j = int(rng.genrand_res53() * ny);
  grow_cluster(i, j, lat[i][j]);

}

void lattice::grow_cluster(int i, int j, int s)
{
  cluster[i][j] = true;
  lat[i][j] = -lat[i][j];

  int iPrev = i == 0 ? nx-1 : i-1;
  int iNext = i == nx-1 ? 0 : i+1;
  int jPrev = j == 0 ? ny-1 : j-1;
  int jNext = j == ny-1 ? 0 : j+1;

  if (!cluster[iPrev][j]) try_add(iPrev, j, s);
  if (!cluster[iNext][j]) try_add(iNext, j, s);
  if (!cluster[i][jPrev]) try_add(i, jPrev, s);
  if (!cluster[i][jNext]) try_add(i, jNext, s);

}

void lattice::try_add(int i, int j, int s)
{
  if (lat[i][j] == s)
  if (rng.genrand_res53() < add_prob)
    grow_cluster(i, j, s);
}

void lattice::wolff_cluster(double K, int cb) {

  add_prob = 1.0 - exp(-2.0*K);
  for (int i=0; i<nx; ++i)
  for (int j=0; j<ny; ++j) {
    cluster[i][j] = false;
  }

  int i, j;
  do {
    i = int(rng.genrand_res53() * nx);
    j = int(rng.genrand_res53() * ny);
  } while ( (i+j)%2!=cb );
  grow_cluster_eo(i, j, lat[i][j]);

}

void lattice::grow_cluster_eo(int i, int j, int s)
{
  cluster[i][j] = true;
  lat[i][j] = -lat[i][j];

  int iPrev = i == 0 ? nx-1 : i-1;
  int iNext = i == nx-1 ? 0 : i+1;
  int jPrev = j == 0 ? ny-1 : j-1;
  int jNext = j == ny-1 ? 0 : j+1;

  if (!cluster[iPrev][jPrev]) try_add_eo(iPrev, jPrev, s);
  if (!cluster[iNext][jPrev]) try_add_eo(iNext, jPrev, s);
  if (!cluster[iPrev][jNext]) try_add_eo(iPrev, jNext, s);
  if (!cluster[iNext][jNext]) try_add_eo(iNext, jNext, s);

}

void lattice::try_add_eo(int i, int j, int s)
{
  if (lat[i][j] == s)
  if (rng.genrand_res53() < add_prob)
    grow_cluster_eo(i, j, s);
}


double lattice::magnetization() {
  int ans=0;
  for (int i=0; i<nx; ++i)
  for (int j=0; j<ny; ++j) {
    ans += lat[i][j];
  }
  return 1.0*ans / (double) (nx*ny);
}

double lattice::magnetization(int cb) {
  int ans=0;
  for (int i=0; i<nx; ++i)
  for (int j=0; j<ny; ++j)
  if( (i+j)%2==cb ) {
    ans += lat[i][j];
  }
  return 2.0*ans / (double) (nx*ny);
}

double lattice::hamiltonian() {
  int ans=0;
  int iNext, jNext;
  for (int i=0; i<nx; ++i)
  for (int j=0; j<ny; ++j) {
    iNext = i+1;
    jNext = j+1;
    ans += lat[i][j]*lat[iNext==nx ? 0 : iNext][j];
    ans += lat[i][j]*lat[i][jNext==ny ? 0 : jNext];
  }
  return ans / (double) (2*nx*ny);
}

double lattice::hamiltonian(int cb) {
  int ans=0;
  int iPrev, jPrev, iNext, jNext;
  for (int i=0; i<nx; ++i)
  for (int j=0; j<ny; ++j)
  if( (i+j)%2==cb ) {

    iPrev = (i==0 ? nx-1 : i-1);
    jPrev = (j==0 ? ny-1 : j-1);
    iNext = (i==nx-1 ? 0 : i+1);
    jNext = (j==ny-1 ? 0 : j+1);

    ans += lat[i][j]*lat[iNext][jNext];
    ans += lat[i][j]*lat[iNext][jPrev];
  }
  return ans / (double) (nx*ny);
}



double lattice::correlator(int dx, int dy)
{

  double corr(0.0);
  
  int k, l;
  for (int i=0; i<nx; ++i)
  for (int j=0; j<ny; ++j) {
    k = i + dx;
    k = (k >= nx ? k-nx : k);
    l = j + dy;
    l = (l >= ny ? l-ny : l);
    corr += lat[i][j]*lat[k][l];
  }

  return corr / (nx*ny); 
}

void lattice::print_lattice()
{
  print_lattice('+',' ','-',':');
}

std::ostream& operator<< (std::ostream & out, lattice & lat)
{
  lat.print_lattice('+',' ','-',':');
  return out;
}

void lattice::print_lattice(char up, char dn, char tb, char lr)
{

  for (int j=0; j<2*ny+3; ++j) { std::cout << tb; } // Top 
  std::cout << std::endl;

  for (int i=0; i<nx; ++i) {
     std::cout << lr << " "; // Left
    for (int j=0; j<ny; ++j) {
      std::cout << ( lat[i][j] == 1 ? up : dn ) << " ";
    }
    std::cout << lr; // Right
    std::cout << std::endl;
  }

  for (int j=0; j<2*ny+3; ++j) { std::cout << tb; } // Bottom 
  std::cout << std::endl;

}

int & lattice::elem(int i, int j) {
  return lat[i][j];
}
 

void lattice::print_cluster()
{
  for (int i=0; i<nx; ++i) {
    for (int j=0; j<ny; ++j) {
      std::cout << ( cluster[i][j] == true ? "+" : " " ) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}
 

bool decimate(const lattice& source, lattice& target, int a, int b)
{

  if (source.nx != 2*target.nx) {
    std::cout << "ERROR : source.nx != 2*target.nx" << std::endl;
    return false;
  }

  if (source.ny != 2*target.ny) {
    std::cout << "ERROR : source.ny != 2*target.ny" << std::endl;
    return false;
  }

  if (a<0||a>1) {
    std::cout << "ERROR : choose a in [0,1]" << std::endl;
    return false;
  }

  if (b<0||b>1) {
    std::cout << "ERROR : choose b in [0,1]" << std::endl;
    return false;
  }

  for (int i=a; i<source.nx; i+=2)
  for (int j=b; j<source.ny; j+=2) {
    target.lat[i/2][j/2] = source.lat[i][j];
  }

  return true;

}

bool refine(const lattice& source, lattice& target, double K0, double K1)
{

  if (2*source.nx != target.nx) {
    std::cout << "ERROR : 2*source.nx != target.nx" << std::endl;
    return false;
  }

  if (2*source.ny != target.ny) {
    std::cout << "ERROR : 2*source.ny != target.ny" << std::endl;
    return false;
  }

  // Compute probabilities for each stage of refinement
  double p0[9];
  for (int s=0; s<9; s+=2) {
     p0[s] = 1.0;
     p0[s] += tanh( K0*(s-4) ); 
     p0[s] *= 0.5;
  }

  double p1[9];
  for (int s=0; s<9; s+=2) {
     p1[s] = 1.0;
     p1[s] += tanh( K1*(s-4) ); 
     p1[s] *= 0.5;
  }

  // Map the ultra-coarse sites first
  for (int i=0; i<source.nx; ++i)
  for (int j=0; j<source.ny; ++j) {
    target.lat[2*i][2*j] = source.lat[i][j];
  }

  // Then perform ultra-coarse to coarse
  target.heat_bath_eeoo(p0, 1);

  // Finally perform coarse to fine
  target.heat_bath_eo(p1, 1);

  return true;

}

bool copy(const lattice& source, lattice& target)
{

  if (source.nx != target.nx) {
    std::cout << "ERROR : source.nx != target.nx" << std::endl;
    return false;
  }

  if (source.ny != target.ny) {
    std::cout << "ERROR : source.ny != target.ny" << std::endl;
    return false;
  }

  for (int i=0; i<source.nx; ++i)
  for (int j=0; j<source.ny; ++j) {
    target.lat[i][j] = source.lat[i][j];
  }

  return true;
}


