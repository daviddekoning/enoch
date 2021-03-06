#include "enoch.hpp"
#include "util/solve.cpp"

Enoch::Enoch(Node **n, Element **e, double **l,
	     int nNodes, int nElements, int nTimeSteps,
	     double dt, double alphaR, double betaR,
	     double alphaN, double betaN) {

  this->nodes    = n;
  this->elements = e;
  this->loads    = l;

  this->nNodes = nNodes;
  this->nElements = nElements;
  this->nTimeSteps = nTimeSteps;

  nDoF = 0;
  for(int i = 0; i < nNodes; i++) {
    nDoF += nodes[i]->nDoF;
  }
  
  K = (double**) malloc(nDoF*sizeof(double*));
  K[0] = (double*) malloc(nDoF*nDoF*sizeof(double));

  Ksac = (double**) malloc(nDoF*sizeof(double*));
  Ksac[0] = (double*) malloc(nDoF*nDoF*sizeof(double));

  Klin = (double**) malloc(nDoF*sizeof(double*));
  Klin[0] = (double*) malloc(nDoF*nDoF*sizeof(double));

  C = (double**) malloc(nDoF*sizeof(double*));
  C[0] = (double*) malloc(nDoF*nDoF*sizeof(double));

  Co = (double**) malloc(nDoF*sizeof(double*));
  Co[0] = (double*) malloc(nDoF*nDoF*sizeof(double));

  for(int i = 0; i < nDoF; i++) {
    K[i]    = K[0] + i*nDoF;
    Ksac[i] = Ksac[0] + i*nDoF;
    Klin[i] = Klin[0] + i*nDoF;
    C[i]    = C[0] + i*nDoF;
    Co[i]   = Co[0] + i*nDoF;
    for(int j = 0; j < nDoF; j++) {
      K[i][j]    = 0.0;
      Ksac[i][j] = 0.0;
      Klin[i][j] = 0.0;
      C[i][j]    = 0.0;
      Co[i][j]   = 0.0;
    }
  }

  m = (double*) malloc(nDoF*sizeof(double));

  
  this->dt = dt;

  alphaRayleigh = alphaR;
  betaRayleigh  = betaR;

  alphaNewark = alphaN;
  betaNewark  = betaN;

}

int Enoch::run() {
  /*
   * Determine how non-linear our model is
   */
  int nNonLinElements = 0;
  
  for(int i = 0; i < nElements; i++)
    if( !elements[i]->isLinear() ) 
      nNonLinElements++;

  // if nNonLinElements > 0, we will perform a nonlinear analysis
  
  /*
   * build maps to translate between local and global coordinates
   */

  /*
   * nodes[i]->map[j] holds the number of the global DoF
   *   corresponding to the jth DoF of node i.
   */
  int globalDoF = 0, localDoF = 0;

  for(int i = 0; i < nNodes; i++)
    for(int j = 0; j < nodes[i]->nDoF; j++)
      nodes[i]->map[j] = globalDoF++;


  /*
   * elements[i]->map[j] holds the number of the global DoF
   *   corresponding to the jth DoF in element i.
   */
  int maxDoFe = 0;
  for(int e = 0, localDoF = 0; e < nElements; e++, localDoF = 0) {

    if(elements[e]->getDoFs() > maxDoFe)
      maxDoFe = elements[e]->getDoFs();
    
    for(int j = 0; j < elements[e]->nNodes; j++)
      for(int k = 0; k < elements[e]->nodes[j]->nDoF; k++)
	elements[e]->map[localDoF++] = elements[e]->nodes[j]->map[k];
  }

  // allocate memory for element [K] matrix
  double **Kelement = (double**) malloc(maxDoFe*sizeof(double*));
  Kelement[0] = (double*) malloc(maxDoFe*maxDoFe*sizeof(double));
  for(int i = 1; i < maxDoFe; i++) {
    Kelement[i] = Kelement[0] + i*maxDoFe;
  }

  /*
   * Compute linear (constant) matrix
   */

  // Klin
  // for each element...
  for(int e = 0; e < nElements; e++) {

    // ignore the nonlinear elements (for now)...
    if(!elements[e]->isLinear()) 
      continue;

    // get the element's K matrix...
    elements[e]->getElementK(Kelement);

    // and add is to the global K matrix
    for(int i = 0; i < elements[e]->getDoFs(); i++)
      for(int j = 0; j < elements[e]->getDoFs(); j++)
	Klin[elements[e]->map[i]][elements[e]->map[j]]
	  += Kelement[i][j];
    
  }

  /*
   * Compute the mass matrix, [m].  Since we will have a diagonal matrix,
   * store it in a vector.
   *
   *  WE ARE ASSUMING that the order in which the nodes appear in
   * the array nodes is the order they are arranged in the vectors
   * of displacement, velocity and acceleration.
   */
  int c = 0;
  for(int i = 0; i < nNodes; i++)
    for(int j = 0; j < nodes[i]->nDoF; j++)
      m[c++] = nodes[i]->mass[j];

  /*
   * Compute the linear damping matrix, [Co].  Use Rayleigh's Method:
   *   [Co] = aphlaRayleigh*[m] + betaRayleigh*[K]
   */
  for(int i = 0; i < nDoF; i++)
    for(int j = 0; j <nDoF; j++)
      if(i == j)
	Co[i][j] = alphaRayleigh*m[i] + betaRayleigh*Klin[i][j];
      else
	Co[i][j] = betaRayleigh*Klin[i][j];




  // variables used during time-history integration
  double **x = (double**) malloc(nDoF*sizeof(double*));
  double **v = (double**) malloc(nDoF*sizeof(double*));
  double **a = (double**) malloc(nDoF*sizeof(double*));
  double **L = (double**) malloc(nDoF*sizeof(double*));
  x[0] = (double*) malloc(2*nDoF*sizeof(double));
  v[0] = (double*) malloc(2*nDoF*sizeof(double));
  a[0] = (double*) malloc(2*nDoF*sizeof(double));
  L[0] = (double*) malloc(2*nDoF*sizeof(double));
  for(int i = 1; i < nDoF; i++) {
    x[i] = x[0] + i*2;
    v[i] = v[0] + i*2;
    a[i] = a[0] + i*2;
    L[i] = L[0] + i*2;
  }

  double *R  = new double[nDoF];
  double *dx = new double[nDoF];
  double *dv = new double[nDoF];
  double *da = new double[nDoF];

  for(int i = 0; i < nDoF; i++)
    for(int j = 0; j < 2; j++) {
      x[i][j] = 0.0;
      v[i][j] = 0.0;
      a[i][j] = 0.0;
      L[i][j] = 0.0;
    }

  int t0 = 0, t1 = 1;

  // output files
  FILE* disp = fopen("disp.csv","w");
  FILE* velo = fopen("velo.csv","w");
  FILE* accel = fopen("accel.csv","w");

  // initialize x[][t], v[][t] and a[][t]
  globalDoF = 0;
  for(int i = 0; i < nNodes; i++)
    for(int j = 0; j < nodes[i]->nDoF; j++) {
      x[globalDoF][t0] = nodes[i]->Xo[j];
      v[globalDoF][t0] = nodes[i]->Vo[j];
      a[globalDoF][t0] = nodes[i]->Ao[j];
      globalDoF++;
    }


  // the vector tells us whether the force or the displacement is contrained
  bool constraint[nDoF];
  
  int forces = 0, disps = 0;
  double F[nDoF], D[nDoF];

  cout << "Processing applied loads and displacements" << endl;
  
  int dof = 0;
  for(int i = 0; i < nNodes; i++) {
    for(int j = 0; j < nodes[i]->nDoF; j++) {
      constraint[dof] = nodes[i]->constrainType[j];

      if(constraint[dof] == Node::FORCE) {
	F[dof] = nodes[i]->constrainValue[j];
	D[dof] = 0;
	forces++;
      } else {
	D[dof] = nodes[i]->constrainValue[j];
	F[dof] = 0;
	disps++;
      }
      dof++;
    }
  }
  
  
  // submatrix of elements w/ constrained force
  double **Kf = (double**) malloc(forces*sizeof(double*));
  Kf[0] = (double*) malloc(forces*forces*sizeof(double));
  for(int i = 0; i < forces; i++) {
    Kf[i] = Kf[0] + i*forces;
  }
  
  // translation vector used to assemble the matrix
  int forceTrans[forces];
  double *subForces = new double[forces];
  double *subDisp = new double[forces];
  
  int forceI = 0, forceJ = 0;
  
  for(int i = 0; i < nDoF; i++) {
    if( constraint[i] == Node::FORCE ) {
      forceTrans[forceI++] = i;
    }
  }
  
  
  if( (forces - forceI) != 0 ) {
    cerr << "something screwy is going down wrt the constrained forces." << endl;
    return(6);
  }
  
    // swap t and dt so that they will get swapped back into the right
  // configuration to start w/
  t0 = t1;
  t1 = 1 - t0;

  /*
   * ******* ******* *     * ***** 
   *    *       *    **   ** *      
   *    *       *    * * * * *      
   *    *       *    *  *  * ****  - history loop starts here 
   *    *       *    *     * *     
   *    *       *    *     * *     
   *    *    ******* *     * ***** 
   */


  //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // start of time-history integration loop
  for(int timeStep = 0; timeStep < nTimeSteps; timeStep++) {


    // 0. Swap time variables
    t0 = t1;
    t1 = 1 - t0;

    // 1. Build current [K]
    //  a) start w/ linear [K]
    for(int i = 0; i < nDoF; i++)
      for(int j = 0; j < nDoF; j++) {
        K[i][j] = Klin[i][j];

	// and the linear [C]
	C[i][j] = Co[i][j];
      }
    
    //  b) add the nonlinear element matrices
    for(int e = 0; e < nElements; e++) {
      // only consider nonlinear elements...
      if(elements[e]->isLinear()) {
	continue;
      }

      // get the matrix...
      //  (nonlinear considerations are deal with here)
      elements[e]->getElementK(Kelement);

      // add it to the global matrix
      for(int i = 0; i < elements[e]->getDoFs(); i++)
	for(int j = 0; j < elements[e]->getDoFs(); j++) {
	  K[elements[e]->map[i]][elements[e]->map[j]] +=
	    Kelement[i][j];

	  // and while we're here we can add it to [c], too
	  C[elements[e]->map[i]][elements[e]->map[j]] +=
	    betaRayleigh*Kelement[i][j];
	}
    }

    // we could maybe consider geometric nonlinearities
    // (but not yet...)
    
        
    // lump LHS of eq'n (1) into Kf[i][j],
    //   partioning as we go
    for(int i = 0; i < forces; i++)
      for(int j = 0; j < forces; j++) {
	Kf[i][j] = K[forceTrans[i]][forceTrans[j]]
		   + C[forceTrans[i]][forceTrans[j]]*alphaNewark/(dt*betaNewark);
	if(i == j)
	  Kf[i][j] += m[forceTrans[i]]/(betaNewark*dt*dt);
      }
	
    // lump RHS of eq'n (1) into R[i],
    //   partitioning as we go
    for(int i = 0; i < forces; i++) {
      R[i] = loads[forceTrans[i]][timeStep] +
	m[forceTrans[i]]*
	( v[forceTrans[i]][t0]/(betaNewark*dt) +
	  a[forceTrans[i]][t0]/(2*betaNewark)    );

      for(int j = 0; j < forces; j++)
	R[i] += C[forceTrans[i]][forceTrans[j]]*
	  (  (alphaNewark/betaNewark)*v[forceTrans[j]][t0] +
	     (alphaNewark/(2*betaNewark) - 1)*dt*a[forceTrans[j]][t0] );
    }

    // solve [Kf]{dx} = {R}

    solve_lup(forces, Kf, subDisp, R);
    
    // reset the dx vector
    // DEAL WITH IMPOSED DISPLACEMENTS !!!!!
    //  --> the imposed displacements need to be accounted for
    //        before the unconstrained displacements are found
    for(int i = 0; i < nDoF; i++)
      dx[i] = 0.0;

    // fill in the unconstrained values from this time step
    for(int i = 0; i < forces; i++)
      dx[forceTrans[i]] = subDisp[i];
   

    for(int i = 0; i < nDoF; i++) {
      // add {dx} to the nodal displacements (x[][t1])
      x[i][t1] = x[i][t0] + dx[i];

      // calc v & a

      dv[i] = (alphaNewark/betaNewark)*dx[i]/dt - 
	(alphaNewark/betaNewark)*v[i][t0] -
	( alphaNewark/(2*betaNewark) - 1)*dt*a[i][t0];

      v[i][t1] = dv[i] + v[i][t0];

      da[i] = loads[i][timeStep];
      for(int j = 0; j < nDoF; j++) {
	da[i] -= C[i][j]*dv[j] + K[i][j]*dx[j];
      }
      da[i] = da[i] / m[i];
      a[i][t1] = da[i] + a[i][t0];
      
      
    }

    // update displacements and velocities in nodes
    for(int i = 0; i < nNodes; i++) {
      for(int j = 0; j < nodes[i]->nDoF; j++) {
	nodes[i]->X[j] = x[nodes[i]->map[j]][t1];
	nodes[i]->V[j] = v[nodes[i]->map[j]][t1];
	nodes[i]->A[j] = a[nodes[i]->map[j]][t1];
      }
    }

    // output data, if necessary
    fprintf(disp, "%lf,",timeStep*dt);
    fprintf(velo, "%lf,",timeStep*dt);
    fprintf(accel,"%lf,",timeStep*dt);
    
    for(int i = 0; i < nDoF; i++) {
      fprintf(disp, "%lf,",x[i][t1]);
      fprintf(velo, "%lf,",v[i][t1]);
      fprintf(accel,"%lf,",a[i][t0]);

    }
    fprintf(disp, "\n");
    fprintf(velo, "\n");
    fprintf(accel,"\n");

  } // end of time-history integration loop

  fclose(disp);
  fclose(velo);
  fclose(accel);

  //  TODO!!
  // should probably free up the memory here, too...

  return 0; 
}
