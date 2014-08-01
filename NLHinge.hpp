#ifndef NLHINGE
#define NLHINGE

#include "util/structure.hpp"

class NLHinge: public Element {
public:
  NLHinge (Node **nodes, int nNodes, double K1, double K2, double My);
  int getElementK(double** K);
  
protected:
  double K1, K2, dPpos, dPneg, dPo, Po, P;
};

#endif
