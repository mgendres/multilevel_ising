#include <vector>
#include "lattice.h"

#ifndef INCLUDED_MULTILEVEL
#define INCLUDED_MULTILEVEL


class multilevel
{
  private:

    std::vector<lattice*> lat;

    multilevel& operator=(const multilevel&);
    multilevel(const multilevel&);

  public:
    explicit multilevel(int, int);
    ~multilevel();

    bool rg_refine(const std::vector<double>, const std::vector<double>); 
    bool rg_refine(const std::vector<double>, const std::vector<double>, int); 
    lattice & operator[](int);
    int levels();

};

#endif
