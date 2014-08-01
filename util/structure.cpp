#include <string.h>
#include <iostream>

#include "structure.hpp"


Node::Node(double *Xo, double *Vo, double *Ao, int nDoF,
	   bool *constrainType, double *constrainValue,
	   double *mass, int number) {
  
  this->number = number;
  this->nDoF = nDoF;

  this->Xo = new double[nDoF];
  this->Vo = new double[nDoF];
  this->Ao = new double[nDoF];

  this->X = new double[nDoF];
  this->V = new double[nDoF];
  this->A = new double[nDoF];
  
  this->constrainType  = new bool[nDoF];
  this->constrainValue = new double[nDoF];
  this->mass           = new double[nDoF];

  this->map = new int[nDoF];

  for(int i = 0; i < nDoF; i++) {
    this->Xo[i] = Xo[i];
    this->Vo[i] = Vo[i];
    this->Ao[i] = Ao[i];

    X[i] = Xo[i];
    V[i] = Vo[i];
    A[i] = Ao[i];

    this->constrainType[i]  = constrainType[i];
    this->constrainValue[i] = constrainValue[i];
    this->mass[i] = mass[i];
  }
}


Node::Node(Node* node) {

  this->number = node->number;
  this->nDoF   = node->nDoF;

  Xo = new double[nDoF];
  Vo = new double[nDoF];
  Ao = new double[nDoF];

  X = new double[nDoF];
  V = new double[nDoF];
  A = new double[nDoF];

  constrainType  = new bool[nDoF];
  constrainValue = new double[nDoF];
  mass           = new double[nDoF];
  
  this->map = new int[nDoF];

  for(int i = 0; i < nDoF; i++) {
    Xo[i] = node->Xo[i];
    Vo[i] = node->Vo[i];
    Ao[i] = node->Ao[i];

    X[i] = node->X[i];
    V[i] = node->V[i];
    A[i] = node->A[i];

    constrainType[i]  = node->constrainType[i];
    constrainValue[i] = node->constrainValue[i];
    this->mass[i] = node->mass[i];
  }
    
}
 
Element::Element(Node **nodes, int nNodes) {

  // does not allocate memory for an array of pointers to nodes
  this->nodes = new Node*[nNodes];
  this->nNodes = nNodes;

  nDoF = 0;

  // count the DoFs
  for(int i = 0; i < nNodes; i++) {
    this->nodes[i] = nodes[i];
    nDoF += nodes[i]->nDoF;
  }

  this->map = new int[nDoF];

  Ko = (double**) malloc(nDoF*sizeof(double*));
  Ko[0] = (double*) malloc(nDoF*nDoF*sizeof(double));

  K = (double**) malloc(nDoF*sizeof(double*));
  K[0] = (double*) malloc(nDoF*nDoF*sizeof(double));
  
  for(int i = 0; i < nDoF; i++) {
    K[i] = K[0] + i*nDoF;
    Ko[i] = Ko[0] + i*nDoF;
    for(int j = 0; j < nDoF; j++) {
      if(j==i)
	K[i][j] = Ko[i][j] = 1.0;
      else
	K[i][j] = Ko[i][j] = 0.0;
    }
  }

}

int Element::getElementK(double** K) {

  for(int i = 0; i < nDoF; i++)
    for(int j = 0; j < nDoF; j++)
      K[i][j] = this->K[i][j];

  return 0;
}

int Element::getElementKo(double** K) {
  memcpy(Ko, this->Ko, nDoF*nDoF*sizeof(double));
  return 0;
}

int Element::getElementForces(double* F) {
  for(int i = 0; i < nDoF; i++) {
    F[i] = 1.0;
  }
  return 0;
}

bool Element::isLinear() {
  return true;
}
