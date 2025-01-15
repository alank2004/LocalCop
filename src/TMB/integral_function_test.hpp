#include "LocalCop/gaussian.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type integral_function_test(objective_function<Type> *obj)
{
  PARAMETER(x);
  PARAMETER(theta);

  Type tiny = 0.0; // Set to 1e-12 for robustness
  //Type ans = LocalCop::IntegralFunctionTest(x, theta) + tiny; // Swap with next line if using tiny_ad code from gaussian.hpp
  Type ans = LocalCop::exponential_evaluator(x, theta) + tiny;
  return ans;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
