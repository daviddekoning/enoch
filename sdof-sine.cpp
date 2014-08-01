#include <argtable2.h>
#include <math.h>

#include "util/structure.hpp"
#include "beam.hpp"
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

  struct arg_dbl *_E = arg_dbl0("E",
			       "youngs-modulus",
			       NULL,
			       "the material stiffness");

  struct arg_dbl *_A = arg_dbl0("A",
			       "area",
			       NULL,
			       "the member's cross-sectional stiffness");

  struct arg_dbl *_I = arg_dbl0("I",
			       "moment-of-interia",
			       NULL,
			       "the member's moment of interia");

  struct arg_dbl *_m = arg_dbl0("m",
			       "mass",
			       NULL,
			       "the mass at the top of the member");

  struct arg_dbl *_F = arg_dbl0("F",
			       "force",
			       NULL,
			       "the force applied at the top of the member");

  struct arg_dbl *_w = arg_dbl0("w",
			       "forcing-frequency",
			       NULL,
			       "the forcing frequency, in rads/sec");


  _dt->dval[0] = 0.1;
  _d->dval[0] = 50;
  _a->dval[0] = 0.0115;
  _b->dval[0] = 0.0115;
  _A->dval[0] = 0.0002;
  _E->dval[0] = 200000000;
  _I->dval[0] = 0.0000000133333;
  _m->dval[0] = 10;
  _F->dval[0] = 40;
  _w->dval[0] = 3.45;


  struct arg_end *end = arg_end(9);
  void* argtable[] = {_dt,_d,_a,_b,_A,_E,_I,_m,_F,_w,end};

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
  fprintf(param,"Time step: %lf\n",_dt->dval[0]);
  fprintf(param,"Duration:  %lf\n",_d->dval[0]);
  fprintf(param,"alpha:     %lf\n",_a->dval[0]);
  fprintf(param,"beta:      %lf\n",_b->dval[0]);
  fprintf(param,"area:      %lf\n",_A->dval[0]);
  fprintf(param,"I:         %lf\n",_I->dval[0]);
  fprintf(param,"E:         %lf\n",_E->dval[0]);
  fprintf(param,"mass:      %lf\n",_m->dval[0]);
  fprintf(param,"F:         %lf\n",_F->dval[0]);
  fprintf(param,"wf:        %lf\n",_w->dval[0]);
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

  // second node has a y coordinate of 2 m 
  x[1] = 2;

  constrainType[0] = Node::FORCE;
  constrainType[1] = Node::FORCE;
  constrainType[2] = Node::FORCE;

  cerr << "about to make the second node" << endl;

  nodes[1] = new Node(x, v, a, 3, constrainType, constrainValue,
		      mass, 1);

  cerr << "made 2nd node" << endl;

  Element** elements = new Element*[1];
  
  elements[0] = new LinearBeam(nodes, 
			       2,
			       _A->dval[0],
			       _E->dval[0],
			       _I->dval[0]);

  cerr << "generated elements" << endl;

  double dt = _dt->dval[0];

  double time = _d->dval[0];
  int nSteps = (int) (time / dt);
  cerr << "Number of Time Steps: " << nSteps << endl;
  double wf = 2;
  cerr << "bo" << endl;
  double **loads;
  cerr << "bi" << endl;
  loads = new double*[6];
  for(int i = 0; i < 6; i++) {
    loads[i] = new double[nSteps];
  }
  cerr << "ba" << endl;

  
  for(int i = 0; i < nSteps; i++) {
    for(int j = 0; j < 6; j++)
      loads[j][i] = 0.0;
    loads[3][i] = _F->dval[0]*(sin(_w->dval[0]*dt*i) - sin(_w->dval[0]*dt*(i-1)));
  }
  
  /*
  loads[3][1] = 0.2;
  loads[3][2] = 0.2;
  */

  
  FILE* load = fopen("load.csv", "w");
  for(int i = 0; i < nSteps; i++) {
    for(int j = 0; j < 6; j++)
      fprintf(load,"%lf,",loads[j][i]);;
    fprintf(load,"\n");;
  }

  Enoch* enoch = new Enoch(nodes,elements,loads,2,1,nSteps, dt,
			   _a->dval[0], _b->dval[0], 1.0, 1.0);

  enoch->run();
}
