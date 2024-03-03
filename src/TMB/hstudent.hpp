#include "LocalCop/student.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template <class Type>
Type hstudent(objective_function<Type> *obj) {
  // R inputs
  DATA_VECTOR(u1);
  DATA_VECTOR(u2);
  DATA_VECTOR(weights)
  PARAMETER_VECTOR(theta);
  PARAMETER_VECTOR(nu);
  // output
  vector<Type> lpart = LocalCop::hstudent(u1, u2, theta, nu, 1);
  lpart.array() *= weights.array();
  return -lpart.sum();
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
