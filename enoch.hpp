#ifndef ENOCH
#define ENOCH

#include "util/structure.hpp"

class Enoch {
public:
  // array of pointers to nodes
  Node** nodes;

  // array of pointers to elements
  Element **elements;

  // 2D array of nodal loads
  double **loads;

  /* creates an instance of Enoch that points to n, e and loads
   *   The instance points to the existing memory, it does not
   * make a copy.
   */
  Enoch(Node **n, Element **e, double **loads,
	int nNodes, int nElements, int nTimeSteps,
	double dt, double alphaR, double betaR,
	double alphaN, double betaN);


  int run();

  inline int getNElements() {return nElements;}

private:
  double dt, alphaRayleigh, betaRayleigh;
  double alphaNewark, betaNewark;

  double **K, **Ksac, **Klin;
  double **C, **Co;
  double  *m;

  int nElements;
  int nNodes;
  int nDoF;
  int nTimeSteps;  

};

#endif
