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

  struct arg_dbl *_h = arg_dbl0("h",
			       "height",
			       NULL,
			       "the height of the member");

  struct arg_int * _n = arg_int0("n",
				 "number-of-storeys",
				 NULL,
				 "the number of storeys to be modelled");

  _dt->dval[0] = 0.1;
  _d->dval[0] = 10;
  _a->dval[0] = 0.0115;
  _b->dval[0] = 0.0115;
  _A->dval[0] = 18737;
  _E->dval[0] = 200;
  _I->dval[0] = 4213000000.0;
  _m->dval[0] = 1;
  _F->dval[0] = 30;
  _h->dval[0] = 4000;
  _n->ival[0] = 2;

  struct arg_end *end = arg_end(11);
  void* argtable[] = {_dt,_d,_a,_b,_A,_E,_I,_m,_F,_h,_n,end};

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
  fprintf(param,"h:         %lf\n",_h->dval[0]);
  fprintf(param,"F:         %lf\n",_F->dval[0]);
  fprintf(param,"n:         %d\n",_n->ival[0]);
  fclose(param);

  double* x = new double[3];
  double* v = new double[3];
  double* a = new double[3];
  bool* constrainType = new bool[3];
  double* constrainValue = new double[3];
  double* mass = new double[3];

  x[0] = x[1] = x[2] = 0.0;
  v[0] = v[1] = v[2] = a[0] = a[1] = a[2] = 0.0;

  constrainType[0] = Node::FORCE;
  constrainType[1] = Node::DISP;
  constrainType[2] = Node::FORCE;

  constrainValue[0] = 0.0;
  constrainValue[1] = 0.0;
  constrainValue[2] = 0.0;

  mass[0] = _m->dval[0];
  mass[1] = _m->dval[0];
  mass[2] = _m->dval[0]/4;

  Node** nodes = new Node*[_n->ival[0]];
  for(int i = 0 ; i < _n->ival[0]; i++) {
    x[1] = i*_h->dval[0];
    nodes[i] = new Node(x, v, a, 3, constrainType, constrainValue,
			mass, i);
  }  

  nodes[0]->constrainType[0] = Node::DISP;
  nodes[0]->constrainType[2] = Node::DISP;


  Element** elements = new Element*[_n->ival[0] - 1];
  Node** elementNodes = new Node*[2];
  for(int i = 0 ; i < _n->ival[0]-1; i++) {
    cout << "generating element " << i << endl;
    elementNodes[0] = nodes[i];
    elementNodes[1] = nodes[i+1];
    elements[i] = new LinearBeam(elementNodes,
				 2,
				 _A->dval[0],
				 _E->dval[0],
				 _I->dval[0]);
    cerr << "elements[" << i << "] goes from (";
    cerr << elements[i]->nodes[0]->X[0] << ",";
    cerr << elements[i]->nodes[0]->X[1] << ") to (";
    cerr << elements[i]->nodes[1]->X[0] << ",";
    cerr << elements[i]->nodes[1]->X[1] << ")" << endl;;
  }  
  


  double dt = _dt->dval[0];

  double time = _d->dval[0];
  int nSteps = (int) (time / dt);
  double wf = 2;
  double **loads;
  loads = new double*[_n->ival[0]*3];
  for(int i = 0; i < _n->ival[0]*3; i++) {
    loads[i] = new double[nSteps];
  }

  double E = _E->dval[0];
  double I = _I->dval[0];
  double L = _h->dval[0];

  double k = 3*E*I/(L*L*L);

  double w = sqrt(k/_m->dval[0]);
  double T = 2*M_PI/w;

  cerr << "Period is " << T << endl;

  for(int i = 0; i < nSteps; i++)
    for(int j = 0; j < _n->ival[0]*3; j++)
      loads[j][i] = 0.0;

  cerr << "apply load over " << 0.1/dt << " time steps." << endl;

  int applySteps =(int) trunc(0.1/dt);

  for(int i = 0; i < applySteps; i++)
    for(int j = 0; j < _n->ival[0]*3 ; j += 3)
      loads[j][i] = _F->dval[0]/(double) applySteps;

  
  FILE* load = fopen("load.csv", "w");
  for(int i = 0; i < nSteps; i++) {
    for(int j = 0; j < _n->ival[0]*3; j++)
      fprintf(load,"%lf,",loads[j][i]);
    fprintf(load,"\n");;
  }

  Enoch* enoch = new Enoch(nodes,elements,loads,_n->ival[0],_n->ival[0]-1,nSteps, dt,
			   _a->dval[0], _b->dval[0], 0.5, 0.25);

  enoch->run();
}
