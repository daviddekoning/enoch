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
 *   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. *
 *   IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR     *
 *   OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, *
 *   ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR *
 *   OTHER DEALINGS IN THE SOFTWARE.                                       *
 ***************************************************************************/

 #include <iostream.h>
 #include <math.h>
 
 /* Solve the system of linear equations F = KD, using the LUP
  * factorization method, (with pivoting).
  */
 int solve_lup(int size, double **K, double *D, double *F) {
   
  /* Data Structures
   * These are  vectors of length size.
  */
  double y[size];
  int pi[size];


  // variables for switching array terms
  double ftemp, p;
  int itemp, i, j, k, kprime;

  // this is a compact representation of the permutation matrix
  // we start off with no permutations
  for(i = 0; i < size; i++) {
    pi[i] = i;
  }

  /*  LUP Decomposition:
      Algorithm from "Introduction to Algorithms" (1st Ed.)
      Cormen, Leiserson & Rivest, pp. 749-761
  */
  
  //time1 = time(0);
  // print the partitioned matrix to the file "Kf.csv"

  FILE* out = fopen("Kf.csv","w");
  for(int a = 0; a < size; a++) {
    for(int b = 0; b < size; b++) {
      fprintf(out,"%lf,",K[a][b]);
    }
    fprintf(out,"\n");
  }
  fclose(out);
  
  for (k = 0; k < size - 1; k++) {

    // start pivoting
    p = 0;
    for (i = k; i < size; i++) {
      if( fabs(K[i][k]) > p ) {
	//printf("%d,%d\n",i,k);
	p = fabs(K[i][k]);
	kprime = i;
      }
    }

    // check for matrix singularity
    if( fabs(p) < 0.000000000001 ) {
      cerr << "Matrix is singular.  oops." << endl;
      cerr << "Singularity encountered @ [" << i << ", " << k << "]" << endl;
      out = fopen("singular.csv","w");
      for(int a = 0; a < size; a++) {
	for(int b = 0; b < size; b++) {
	  fprintf(out,"%lf,",K[a][b]);
	}
	fprintf(out,"\n");
      }
      fclose(out);
      return(1);
    }
    

    // exchange pi[k] & pi[kprime]
    itemp = pi[k];
    pi[k] = pi[kprime];
    pi[kprime] = itemp;
    
    // exchange rows k & kprime of matrix
    for(i = 0; i < size; i++) {
      ftemp = K[k][i];
      K[k][i] = K[kprime][i];
      K[kprime][i] = ftemp;
    }
    
    // finish pivoting

    // modify matrix terms:
    for(i = k + 1; i < size; i++) {
      K[i][k] = K[i][k]/K[k][k];
      for(j = k + 1; j < size; j++) {
	K[i][j] -= K[i][k]*K[k][j];
      }
    }
  }
  
  /* Forward substitution to find intermidiate sol'n y:
     Ly = Pb
  */
  for(i = 0; i < size; i++) {
    y[i] = F[ pi[i] ];
    for(j = 0; j < i; j++) {
      y[i] -= K[i][j]*y[j];
    }
  }
  
  /* Back substitution to find final sol'n D:
     Ux = y
  */
  for(i = size - 1; i >= 0; i--) {
    ftemp = y[i];
    for(j = i+1; j < size; j++) {
      ftemp -= K[i][j]*D[pi[j]];
    }
    D[pi[i]] = ftemp / K[i][i];
    //cout << " | " << D[i] << " | ";
  }
  
  return 0;  

 }
 
 /* Solve the system of linear equations F = KD, using the LU
  * factorization method, (without pivoting).
  *
  * If you know that your matrix is +ve definite (I think that's the
  * criteria) then this function will be marginally faster than 
  * the solve_lup function.
  */
 int solve_lu
 (int size, double **K, double *D, double *F) {
 
 
  /* Data Structures
     These are all SIZE x SIZE arrays.
  */
  double *y;

  // variables for switching array terms
  double ftemp;
  
  // Allocate Memory:
  y = new double[size];

  if( y == NULL) {
    printf("Not enougth memory available to solve system.\n");
    return 1; 
  }

  /*  LU Decomposition:
      Algorithm from "Introduction to Algorithms" (1st Ed.)
      Cormen, Leiserson & Rivest, pp. 749-761
      
  */
  
  for (int k = 0; k < size - 1; k++) {

    // modify matrix terms:
    for(int i = k + 1; i < size; i++) {
      K[i][k] = K[i][k]/K[k][k];
      for(int j = k + 1; j < size; j++) {
	K[i][j] -= K[i][k]*K[k][j];
      }
    }
  }
  
  
  /* Forward substitution to find intermidiate sol'n y:
     Ly = Pb
  */
  for(int i = 0; i < size; i++) {
    y[i] = F[ i ];
    for(int j = 0; j < i; j++) {
      y[i] -= K[i][j]*y[j];
    }
  }
  
  /* Back substitution to find final sol'n x:
     Ux = y
  */
  for(int i = size - 1; i >= 0; i--) {
    ftemp = y[i];
    for(int j = i+1; j < size; j++) {
      ftemp -= K[i][j]*D[j];
    }
    D[i] = ftemp / K[i][i];
  }

  delete y;
  
   return 0;
 }
 
 /* Solve the system of linear equations F = KD, using the Cholesky
  * factorization method, without pivoting.
  *
  * The matrix K must be positive definite for this function to work.
  */
 int solve_cholesky(int size, double **K, double *D, double *F) {
 
 
  /* Data Structures
     These are all SIZE x SIZE arrays.
  */
  double *y;

  // variables for switching array terms
  double ftemp;
  
  // Allocate Memory:
  y = new double[size];

  
  if( y == NULL) {
    printf("Not enougth memory available to solve system.\n");
    return 1; 
  }

  /*  Cholesky decomposition
   *
   */
  
  for (int i = 0; i < size; i++) {
    // off-diagonal terms...
    for(int j = i + 1; j < size; j++) {
      //K[j][i] = K[j][i];
      for(int k = 0; k < i; k++) {
	K[j][i] -= K[j][k]*K[i][k];
      }
      K[j][i] = K[j][i]/K[i][i];
    }
    
    // diagonal terms...
    //K[i][i] = K[i][i];
    for(int k = 0; k < i; k++) {
      K[i][i] -= K[i][k]*K[i][k];
    }
    K[i][i] = sqrt(K[i][i]);
  }
  
  // fill upper triangle w/ transpose(L)
  for(int i = 0; i < size; i++) {
    for(int j = i+1; j < size; j++) {
      K[i][j] = K[j][i];
    }
  }
  
  /* Forward substitution to find intermidiate sol'n y:
     Ly = Pb
  */
  for(int i = 0; i < size; i++) {
    y[i] = F[ i ];
    for(int j = 0; j < i; j++) {
      y[i] -= K[i][j]*y[j];
    }
    y[i] = y[i] / K[i][i];
  }
  
  /* Back substitution to find final sol'n x:
     Ux = y
  */
  for(int i = size - 1; i >= 0; i--) {
    ftemp = y[i];
    for(int j = i+1; j < size; j++) {
      ftemp -= K[i][j]*D[j];
    }
    D[i] = ftemp / K[i][i];
  }

  delete y;
  
   return 0;
 }
 
 /* Solve the system of linear equations F = KD, using Gaussian
  * elimination
  */
 int solve_gauss(int size, double **K, double *D, double *F) {
   
  double **temp;
  double ratio;

  // allocate memory for 'temp'
  temp = (double**) malloc(size*sizeof(double*));
  temp[0] = (double*) malloc(size*size*sizeof(double));
  for(int i = 1; i < size; i++) {
    temp[i] = temp[0] + (size+1)*i;
  }


  // populate 'temp'  
  for(int i = 0; i < size; i++) {
    temp[i][size] = F[i];
    for(int j = 0; j < size; j++)
      temp[i][j] = K[i][j];
  }
   
   
   // reduce matrix to upper triangular form, keeping the vector B consistent
   for(int j = 0; j < size; j++) {

     for(int i = j; i < size; i++) {

       if(i == j) {
	 ratio = 1/temp[i][j];
	 for(int j2 = 0; j2 < size+1; j2++)
	   temp[i][j2] = temp[i][j2]*ratio;
       } else {
	 ratio = temp[i][j]/temp[j][j];
	 for(int j2 = 0; j2 < size+1; j2++)
	   temp[i][j2] = temp[i][j2] - temp[j][j2]*ratio;	   
       }
     }
   }

   // back substitute to find x[i]
   for (int i = size - 1; i >= 0; i--) {
     D[i] = temp[i][size];
     for(int j = size-1; j > i; j--)
       D[i] -= temp[i][j]*D[j];
   }

   free(temp[0]);
   free(temp);
   
    
   return 0;
 }
  
 /* Solve the system of linear equations F = KD, using Gauss-Jordan
  * elimination.
  */
 int solve_gaussjordan(int size, double **K, double *D, double *F) {
 
 
   double **temp;
   double ratio;

   temp = (double**) malloc(size*sizeof(double*));
   temp[0] = (double*) malloc(size*(size+1)*sizeof(double));
   for(int i = 1; i < size; i++)
     temp[i] = temp[0] +(size+1)*i;

   // Set up vectors
   for(int i = 0; i < size; i++) {
     temp[i][size] = F[i];
   }

   // set up matixes
   for(int i = 0; i < size; i++)
     for(int j = 0; j < size; j++)
       temp[i][j] = K[i][j];

   // solve for x w/ Gauss-Jordan elimination
   for(int j = 0; j < size; j++) {

     for(int i = 0; i < size; i++) {

       if(i == j) {
	 ratio = 1/temp[i][j];
	 for(int j2 = 0; j2 < size+1; j2++)
	   temp[i][j2] = temp[i][j2]*ratio;
       } else {
	 ratio = temp[i][j]/temp[j][j];
	 for(int j2 = 0; j2 < size+1; j2++)
	   temp[i][j2] = temp[i][j2] - temp[j][j2]*ratio;	   
       }
     }
   }

   for(int i = 0; i < size; i++) {
     D[i] = temp[i][size];
   }

   free(temp[0]);
   free(temp);
   
 
   return 0;
 }
 
  /*
  * Solve the system of linear eq'ns F = KD, where one of each F[i] and D[i]
  * is known, as specified by the vector constraints[i].  If constraints[i]
  * is false, then F[i] is known and D[i] is unknown.  If constraints[i]
  * is true, then F[i] is unknown and D[i] is known.   K is the parameter used by the Paynes &
  * Irons algorithm.  This method used the solve_gauss function to find a sol'n.
  *
  * Note that the K matrix is modified during this process in such a way that the results
  * contained in F once the function returns will be wrong.  Keep another copy of K and F,
  * and use them, along with the computed values of D to find the nodal forces.
  */         
 int solve_paynesnirons(int size, double **K, double *D, double *F, bool *constraints, double Kpi) {

  int returnValue;
  
  double **KK = (double**) malloc(size*sizeof(double*));
  KK[0] = (double*)  malloc(size*size*sizeof(double));
  for(int i = 0; i < size; i++) {
    KK[i] = KK[0] + size*i;
  }
  
  
  for(int i = 0; i < size; i++)
    for(int j = 0; j < size; j++)
      KK[i][j] = K[i][j];


 
    for(int i = 0; i < size; i++) {
      if(constraints[i] == true) {
        // F[i] is unknown, D[i] known
        KK[i][i] *= Kpi;
        F[i] = D[i]*Kpi;
      }
    }
    
    returnValue =  solve_gauss(size, KK, D, F);
    
    free(KK[0]);
    free(KK);
    
    return returnValue;
    
 }
 
 
 /*
  * This method will solve the system F = KD using the iterative Gauss-Seidel
  * method.  It assumes that matrices F & K are known and D is unknown.
  *
  *  It continues until the relative change in each element of D is less than maxError,
  * and uses the relaxation factor w.
  */
int solve_gaussseidel(int size, double **K, double *D, double *F, double w, double maxError) {

  bool done = false;
  double *error = (double*) malloc(size*sizeof(double));
  double ainverse = 1, previous = 0;
  int iteration = 0;
  
  // arbitrary starting point for D[i]
  // start w/ error greater than the cuttoff (maxError)
  for(int i = 0; i < size; i++) {
    D[i] = 1;
    error[i] = maxError + 1;
  }
  
  while(!done) {
    iteration++;
    //cerr << "iteration: " << iteration++ << endl;
    for(int i = 0; i < size; i++) {
      previous = D[i];
      ainverse = 1/K[i][i];
      
      // find the new value of D[i]
      D[i] = F[i]*ainverse;
      for(int j = 0; j < size; j++)
        if(i != j) {
          D[i] -= K[i][j]*ainverse*D[j];
        }
      
      // over relax the value of D[i]
      D[i] = previous*(1-w) + D[i]*w;
        
      // calulate the relative change in D[i]
      error[i] = (D[i] - previous)/D[i];
      
      // is this term okay?
      if( error[i] > maxError)
        continue;
      
      // if so, are the rest ok?
      done = true;
      for(int j = 0; j < size; j++)
        if(error[j] > maxError) {
          done = false;
          break;
        }
      
      if(done)
        break;
      
    }
  
  }
  
  cerr << iteration << " iterations" << endl;
  
  free(error);
  
  return 0;
  
}


