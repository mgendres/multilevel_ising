#include <string>
#include <vector>

#ifndef INCLUDED_INTERPOLATOR
#define INCLUDED_INTERPOLATOR

class interpolator
{
  private:

    std::vector<double> f;
    std::vector<double> x;

    interpolator& operator=(const interpolator&);
    interpolator(const interpolator&);

  public:
    explicit interpolator(std::string);
    ~interpolator();

    double operator()(double);

};

#endif
