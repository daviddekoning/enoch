#include <argtable2.h>
#include <math.h>

#include "util/structure.hpp"
#include "NLHingeNonLin.hpp"
#include "enoch.hpp"

int main(int argc, char *argv[]) {

  struct arg_dbl *_dt = arg_dbl0("t",
			       "time-step",
			       NULL,
			       "time step for time history integration");

  struct arg_dbl *_d = arg_dbl0("d",
			       "duration",
			       NULL,
			       "the duration of the analysis");

  struct arg_dbl *_a = arg_dbl0("a",
			       "alpha",
			       NULL,
			       "Reynolds damping alpha");

  struct arg_dbl *_b = arg_dbl0("b",
			       "beta",
			       NULL,
			       "Reynolds damping beta");

  struct arg_dbl *_k1 = arg_dbl0("K",
			       "initial-stiffness",
			       NULL,
			       "the initial stiffness");

  struct arg_dbl *_k2 = arg_dbl0("k",
			       "post-yield-stiffness",
			       NULL,
			       "the member's cross-sectional stiffness");

  struct arg_dbl *_My = arg_dbl0("y",
			       "yield-moment",
			       NULL,
			       "the hinge's yield moment");

  struct arg_dbl *_m = arg_dbl0("m",
			       "mass",
			       NULL,
			       "the mass at the top of the structure");

  struct arg_dbl *_M = arg_dbl0("M",
			       "moment",
			       NULL,
			       "the moment applied at the top of the hinge");

  struct arg_file *_eqfile = arg_file1("f",
				       "eq-file",
				       NULL,
				       "the file that contains the earthquake acceleration record");


  _dt->dval[0] = 0.1;
  _d->dval[0] = 10;
  _a->dval[0] = 0.0115;
  _b->dval[0] = 0.0115;
  _k1->dval[0] = 50000000;
  _k2->dval[0] = 3000000;
  _My->dval[0] = 25000;
  _m->dval[0] = 1;
  _M->dval[0] = 40000;


  struct arg_end *end = arg_end(10);
  void* argtable[] = {_dt,_d,_a,_b,_k1,_k2,_My,_m,_M,_eqfile,end};

  const char* progname = "sdof-step";
  int exitcode = 0, returnvalue = 0;
  int nerrors = arg_nullcheck(argtable);


  if(nerrors != 0) {
    cerr << "Not enough memory to proceed" << endl;
    arg_print_errors(stderr,end,progname);
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
    return nerrors;
  }

  nerrors = arg_parse(argc, argv, argtable);

  if(nerrors != 0) {
    arg_print_errors(stderr,end,progname);
    cerr << "Usage: " << endl;
    arg_print_syntaxv(stderr,argtable,"\n");
    return(nerrors);
  }

  FILE* param = fopen("parameters.txt","w");
  fprintf(param,"Input Parameters:\n");
  fprintf(param,"Time step:  %lf\n",_dt->dval[0]);
  fprintf(param,"Duration:   %lf\n",_d->dval[0]);
  fprintf(param,"alpha:      %lf\n",_a->dval[0]);
  fprintf(param,"beta:       %lf\n",_b->dval[0]);
  fprintf(param,"My:         %lf\n",_My->dval[0]);
  fprintf(param,"k1:         %lf\n",_k1->dval[0]);
  fprintf(param,"k2:         %lf\n",_k2->dval[0]);
  fprintf(param,"mass:       %lf\n",_m->dval[0]);
  fprintf(param,"Mmax:       %lf\n",_M->dval[0]);
  fprintf(param,"E/Q record: %s\n",_eqfile->filename[0]);
  fclose(param);

  double* x = new double[3];
  double* v = new double[3];
  double* a = new double[3];
  bool* constrainType = new bool[3];
  double* constrainValue = new double[3];
  double* mass = new double[3];

  Node** nodes = new Node*[2];
  
  x[0] = x[1] = x[2] = 0.0;
  v[0] = v[1] = v[2] = a[0] = a[1] = a[2] = 0.0;

  constrainType[0] = Node::DISP;
  constrainType[1] = Node::DISP;
  constrainType[2] = Node::DISP;

  constrainValue[0] = 0.0;
  constrainValue[1] = 0.0;
  constrainValue[2] = 0.0;

  mass[0] = _m->dval[0];
  mass[1] = _m->dval[0];
  mass[2] = _m->dval[0];

  nodes[0] = new Node(x, v, a, 3, constrainType, constrainValue,
		      mass, 0);

  // set y coord of second node

  constrainType[0] = Node::DISP;
  constrainType[1] = Node::DISP;
  constrainType[2] = Node::FORCE;

  cerr << "about to make the second node" << endl;

  nodes[1] = new Node(x, v, a, 3, constrainType, constrainValue,
		      mass, 1);

  cerr << "made 2nd node" << endl;

  Element** elements = new Element*[1];
  
  elements[0] = new NLHingeNonLin(nodes, 
				  2,
				  _k1->dval[0],
				  _k2->dval[0],
				  _My->dval[0]);
  
  cerr << "generated elements" << endl;

  double dt = _dt->dval[0];

  int nSteps;

  FILE *eq = fopen(_eqfile->filename[0],"r");
  fscanf(eq,"%d",&nSteps);

  double time = dt*nSteps;


  cerr << "Number of Time Steps: " << nSteps << endl;
  cerr << "Duration: " << time << endl;

  double **loads;

  loads = new double*[6];
  for(int i = 0; i < 6; i++) {
    loads[i] = new double[nSteps];
  }

  double accel;

  delete [] mass;

  mass = new double[6];

  mass[0] = mass[1] = mass[2] = mass[3] = mass[4] = 0;
  mass[5] = _m->dval[0];

  for(int j = 0; j < 6; j++)
    loads[j][0] = 0;

  for(int i = 1; i < nSteps; i++) {
    fscanf(eq,"%lf",&accel);
    for(int j = 0; j < 6; j++)
      loads[j][i] = -1*mass[j]*accel - loads[j][i-1];
  }

  
  FILE* load = fopen("load.csv", "w");
  for(int i = 0; i < nSteps; i++) {
    for(int j = 0; j < 6; j++)
      fprintf(load,"%lf,",loads[j][i]);
    fprintf(load,"\n");
  }

  fclose(load);
  fclose(eq);

  Enoch* enoch = new Enoch(nodes,elements,loads,2,1,nSteps, dt,
			   _a->dval[0], _b->dval[0], 0.5, 0.25);

  enoch->run();
}
