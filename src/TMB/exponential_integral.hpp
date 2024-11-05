#include "LocalCop/gaussian.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type exponential_integral(objective_function<Type> *obj)
{
  DATA_VECTOR(lambda); // TODO: Provided code suggests it be PARAMETER_VECTOR(lambda) and etc.
  DATA_VECTOR(lower);
  DATA_VECTOR(upper);

  Type ans = 0;
  Type tiny = 0.0; // Set to 1e-12 for robustness
  for(int i=0; i < lambda.size(); i++) {
    ans -= log( LocalCop::ExponentialIntegral(lambda(i), lower(i), upper(i)) + tiny );
  }
  return ans;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
