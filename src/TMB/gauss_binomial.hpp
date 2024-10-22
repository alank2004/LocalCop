#include "LocalCop/gaussian.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type gauss_binomial(objective_function<Type> *obj)
{
  DATA_VECTOR(x);
  DATA_VECTOR(n);
  DATA_MATRIX(A);
  PARAMETER_VECTOR(b);
  vector<Type> mu = A * b;
  PARAMETER(logsd);

  Type sd = exp(logsd);
  Type ans = 0;
  Type tiny = 0.0; // Set to 1e-12 for robustness
  for(int i=0; i < x.size(); i++) {
    ans -= log( LocalCop::GaussBinomial(x(i), n(i), mu(i), sd) + tiny );
  }
  return ans;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
