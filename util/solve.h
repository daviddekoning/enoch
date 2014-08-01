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

 /* Solve the system of linear equations F = KD, using the Cholesky
  * factorization method, with pivoting.
  */
 int solve_choleskyP(int size, double **K, double *D, double *F);
 
 /* Solve the system of linear equations F = KD, using the Cholesky
  * factorization method, without pivoting.
  */
 int solve_cholesky(int size, double **K, double *D, double *F);
 
 /* Solve the system of linear equations F = KD, using Gaussian
  * elimination
  */
 int solve_gauss(int size, double **K, double *D, double *F);
  
 /* Solve the system of linear equations F = KD, using Gauss-Jordan
  * elimination.
  */
 int solve_gaussjordan(int size, double **K, double *D, double *F);

 
 /*
  * Solve the system of linear eq'ns F = KD, where one of each F[i] and D[i]
  * is known, as specified by the vector constraints[i].  If constraints[i]
  * is false, then F[i] is known and D[i] is unknown.  If constraints[i]
  * is true, then D[i] is unknown.   K is the parameter used by the Paynes &
  * Irons algorithm.  This method used the solve_choleskyP function to find a sol'n.
  */         
 int solve_paynesnirons(int size, double **K, double *D, double *F, bool constraints, double Kpi);
 
 
 
 