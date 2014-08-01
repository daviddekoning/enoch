#include "util/structure.hpp"

class LinearBeam: public Element {
public:
  LinearBeam(Node **nodes, int nNodes, double A, double E, double I);
  
protected:
  double A, E, I, L;
};
