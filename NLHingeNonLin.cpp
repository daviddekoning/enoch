#include <math.h>

#include "NLHingeNonLin.hpp"

NLHingeNonLin::NLHingeNonLin(Node **nodes, int nNodes,
		       double K1, double K2, double My) : Element(nodes, nNodes) {
  if(nNodes != 2) {
    // BORK!!
    cerr << "Non-linear hinge objects can only be associated w/ two nodes" 
<< endl;
    exit(1);
  }

  for(int i = 0; i < nNodes; i++) {
    if(nodes[i]->nDoF != 3) {
      // BORK!!
      cerr << "Non-linear hinges must connect nodes w/ 3 DoFs!" << endl;
      exit(2);
    }
  }

  this->K1 = K1;
  this->K2 = K2;
  this->dP = (My/K1);


  // sanity check
  if(nDoF != 6) {
    // BORK !!!
    cerr << "NLHingeNonLin really ought to have 6 DoFs" << endl;
    exit(1);
  }

  double dx = nodes[1]->Xo[0] - nodes[0]->Xo[0];
  cerr << "dx = " << dx << endl;
  double dy = nodes[1]->Xo[1] - nodes[0]->Xo[1];
  cerr << "dy = " << dy << endl;
  double L = sqrt(dx*dx + dy*dy);
  cerr << "length of element is " << L << endl;

  K[0][0] = K[1][1] = K[3][3] = K[4][4] = 1000000000000.0;
  K[0][3] = K[3][0] = K[1][4] = K[4][1] = 1000000000000.0;
  K[2][2] = K[5][5] =  K1;
  K[2][5] = K[5][2] = -K1;

  cerr << K[0][0] << endl;


  for(int i = 0; i < nDoF; i++)
    for(int j = 0; j < nDoF; j++)
      Ko[i][j] = K[i][j];



  for(int i = 0; i < nDoF; i++) {
    for(int j = 0; j < nDoF; j++)
      cerr << Ko[i][j] << ", ";
    cerr << endl;
  }


}

int NLHingeNonLin::getElementK(double** K) {

  P = nodes[1]->X[2] - nodes[0]->X[2];

  if (P > (Po + dP)){
    Po = P - dP;
    this->K[2][2] = this->K[5][5] =  K2;
    this->K[2][5] = this->K[5][2] = -K2;
  }else if(P < (Po - dP)){
    Po = P + dP;
    this->K[2][2] = this->K[5][5] =  K2;
    this->K[2][5] = this->K[5][2] = -K2;
  }else{
    this->K[2][2] = this->K[5][5] =  K1;
    this->K[2][5] = this->K[5][2] = -K1;
  }

  cout << ".";

  for(int i = 0; i < nDoF; i++)
    for(int j = 0; j < nDoF; j++)
      K[i][j] = this->K[i][j];

  return 0;
}

bool NLHingeNonLin::isLinear() {
  return false;
}



