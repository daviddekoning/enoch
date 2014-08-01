#include <math.h>

#include "util/structure.hpp"
#include "beam.hpp"

LinearBeam::LinearBeam(Node **nodes, int nNodes,
		       double A, double E, double I) : Element(nodes, nNodes) {
  if(nNodes != 2) {
    // BORK!!
    cerr << "Linear beam objects can only be associated w/ two nodes" << endl;
    exit(1);
  }
  
  for(int i = 0; i < nNodes; i++) {
    if(nodes[i]->nDoF != 3) {
      // BORK!!
      cerr << "Linear beams must connect nodes w/ 3 DoFs!" << endl;
      cerr << "  problem w/ node " << i << endl;
      exit(2);
    }
  }

  this->A = A;
  this->E = E;
  this->I = I;

  // sanity check
  if(nDoF != 6) {
    // BORK !!!
    cerr << "LinearBeam really ought to have 6 DoFs" << endl;
    exit(1);
  }
  
  double dx = nodes[1]->Xo[0] - nodes[0]->Xo[0];
  double dy = nodes[1]->Xo[1] - nodes[0]->Xo[1];
  double L = sqrt(dx*dx + dy*dy);
  
  K[0][0] = K[3][3] = E*A/L;
  K[0][3] = K[3][0] = -1*K[0][0];
  K[1][1] = K[4][4] = 12*E*I/(L*L*L);
  K[1][4] = K[4][1] = -1*K[1][1];
  K[2][5] = K[5][2] = 2*E*I/L;
  K[2][2] = K[5][5] = 2*K[2][5];
  K[2][1] = K[1][2] = K[5][1] = K[1][5] = 6*E*I/(L*L);
  K[2][4] = K[4][2] = K[4][5] = K[5][4] = -1*K[1][2];


  cerr << K[0][0] << endl;

  double cos = dx / L;
  double sin = dy / L;

  double** gamma = new double*[nDoF];
  for(int i = 0; i < nDoF; i++)
    gamma[i] = new double[nDoF];

  for(int i = 0; i < nDoF; i++)
    for(int j = 0; j < nDoF; j++)
      gamma[i][j] = 0.0;

  gamma[0][0] = gamma[1][1] = gamma[3][3] = gamma[4][4] = cos;
  gamma[0][1] = gamma[3][4] = sin;
  gamma[1][0] = gamma[4][3] = -1*sin;
  gamma[2][2] = gamma[5][5] = 1.0;

  for(int i = 0; i < nDoF; i++)
    for(int j = 0; j < nDoF; j++)
      Ko[i][j] = 0.0;


  // transform K into global coordinates
  // this is a double sum (j & k) on each matrix term (i & l)
  for(int i = 0; i < nDoF; i++)
    for(int j = 0; j < nDoF; j++)
      for(int k = 0; k < nDoF; k++)
	for(int l = 0; l < nDoF; l++)
	  Ko[i][l] += gamma[j][i]*K[j][k]*gamma[k][l];

  for(int i = 0; i < nDoF; i++) {
    for(int j = 0; j < nDoF; j++)
      cerr << Ko[i][j] << ", ";
    cerr << endl;
  }


  // Ko is now in global coordinates
  // copy it over to K
  for(int i = 0; i < nDoF; i++) {
    for(int j = 0; j < nDoF; j++) {
      K[i][j] = Ko[i][j];
    }
  }


  
  
  
}

