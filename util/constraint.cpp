#include "constraint.hpp"

#include <math.h>

StaticConstraint::StaticConstraint(bool t, double v) : Constraint(t) {
  value = v;
}


ConstantDynamicConstraint::ConstantDynamicConstraint(bool type, double To,
						     double v) : 
  DynamicConstraint(type) {
  
  value = v;
  startTime = To;
  
}


double ConstantDynamicConstraint::getValue(double t) {
  if(t > startTime)  return value;
  else return 0;
}


double ConstantDynamicConstraint::getdValue(double t, double dt) {
  // if the start time occurs between t and t + dt, there is a change in
  // the constraint value
  if( fabs(startTime - t) <= dt ) return value;
  else return 0;
}

SinusoidalConstraint::SinusoidalConstraint(bool type, double A,
						  double w, double phi) : 
  DynamicConstraint(type) {
  
  this->A = A;
  this->w = w;
  this->phi = phi;
}

double SinusoidalConstraint::getValue(double t) {
  return A*sin(w*t + phi);
}

double SinusoidalConstraint::getdValue(double t, double dt) {
  return A*( sin(w*(t+dt) + phi) - sin(w*t + phi) );
}
