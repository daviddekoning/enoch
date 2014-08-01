/***************************************************************************
 *   Copyright (C) 2005 by David de Koning                                 *
 *   david@dekoning.ca                                                     *
 *                                                                         *
 *   Permission is hereby granted, free of charge, to any person obtaining *
 *   a copy of this software and associated documentation files (the       *
 *   "Software"), to deal in the Software without restriction, including   *
 *   without limitation the rights to use, copy, modify, merge, publish,   *
 *   distribute, sublicense, and/or sell copies of the Software, and to    *
 *   permit persons to whom the Software is furnished to do so, subject to *
 *   the following conditions:                                             *
 *                                                                         *
 *   The above copyright notice and this permission notice shall be        *
 *   included in all copies or substantial portions of the Software.       *
 *                                                                         *
 *   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,       *
 *   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF    *
 *   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.*
 *   IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR     *
 *   OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, *
 *   ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR *
 *   OTHER DEALINGS IN THE SOFTWARE.                                       *
 ***************************************************************************/
#ifndef STRUCTURE
#define STRUCTURE

#include <string.h>
#include <iostream.h>

class Node {
public:
  Node(double *Xo, double *Vo, double *Ao, int nDoF,
       bool *constrainType, double *constrainValue,
       double *mass, int number);
  Node(Node * node);
  double *Xo, *Vo, *Ao, *X, *V, *A;
  int nDoF;
  double *mass;

  /*
   * Array holding the mapping of this node to the
   * global DoFs.
   */
  int *map;

  /*
   * The node's number.  If the node is stored in an array, this should
   * be the same as the node's index in the array.
   */
  int number;
  bool   *constrainType;
  double *constrainValue;
  static const bool FORCE = false, DISP = true;
};

class Element {
public:
  Element(Node **nodes, int nNodes);
  int nNodes;
  Node **nodes;

  /*
   * Array holding the mapping of this element to the
   * global DoFs.
   */
  int *map;

  /*
   * The elements identification number.  Useful for postprocessing.
   */
  int number;

  virtual int getElementK(double** K);
  virtual int getElementKo(double** Ko);
  virtual int getElementForces(double* F);


  virtual bool isLinear();

  int getDoFs() {return nDoF;};

protected:
  double **Ko;
  double **K;
  int nDoF;
};


#endif 
 
 
 
 
 
 
 
