#include <ostream>

#ifndef INCLUDED_LATTICE
#define INCLUDED_LATTICE

enum lat_state { ORDERED, DISORDERED };

class lattice
{
  private:

    int nx;
    int ny;
    int vol;
    int** lat;
    bool** cluster;

    lattice& operator=(const lattice&);
    lattice(const lattice&);

    // Heat bath functions
    void heat_bath_eo(double *, int);
    void heat_bath_eo(double **, int);

    void heat_bath_eeoo(double *, int); 
    void heat_bath_eooe(double *, int);
    void heat_bath_eeoo(double **, int); 
    void heat_bath_eooe(double **, int);

    // Wolff cluster functions
    double add_prob;
    void grow_cluster(int, int, int); 
    void try_add(int, int, int);

    void grow_cluster_eo(int, int, int); 
    void try_add_eo(int, int, int);


  public:
    explicit lattice(int, int);
    ~lattice();

    // Functions that act on entire lattice
    void init(enum lat_state); 
    void heat_bath(double);
    void heat_bath(double, double);
    void wolff_cluster(double);

    // Observables, etc.
    double magnetization();
    double hamiltonian();
    double correlator(int,int);
    void print_lattice();
    
    void print_lattice(char, char, char, char);
    void print_cluster();

    // Functions that act on a e/o sublattices
    void heat_bath(double,int);
    void heat_bath(double, double, int);
    void wolff_cluster(double, int);
    double magnetization(int);
    double hamiltonian(int);

    int &elem(int i, int j);


    friend std::ostream& operator<< (std::ostream &, lattice &);

    friend bool copy(const lattice&, lattice&);
    friend bool decimate(const lattice&, lattice&, int, int);
    friend bool refine(const lattice&, lattice&, double, double);

};

#endif
