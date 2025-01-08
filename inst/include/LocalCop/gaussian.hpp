/// @file gaussian.hpp

#ifndef LOCALCOP_GAUSSIAN_HPP
#define LOCALCOP_GAUSSIAN_HPP

// this is where RefVector_t etc. is defined
#include "config.hpp"

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI  0.918938533204672741780329736406
/* log(sqrt(2*pi))
== log(2*pi)/2 */
#endif

namespace LocalCop {

  // TODO: This code is a TMB example, copied here to test how it works. Need to replace it with pgaussian.
  /*
    This class evaluates the integral of the exponential density with rate parameter lambda from lower to upper
  */
  template<class Float>
  struct Exponential {
    typedef Float Scalar; // Required by integrate
    Float x, y, theta;         // Parameters 
    // Evaluate z * y * theta
    Float operator() (Float z) {
      Float ans = z * y * theta;
      return ans;
    }
    // Integrate latent variable z out
    Float integrate() {
      using gauss_kronrod::integrate;
      Float ans =
        integrate(*this, 0, x);
      return ans;
    }
  };
  // ****** How to use it in TMB:
  // 1. Create an evaluator 'eval' for previous class
  template<class Float>
  Float eval(Float x, Float y, Float theta) {
    Exponential<Float> f = {x, y, theta};
    return f.integrate();
  }
  // 2. Run 'eval' through tiny_ad and obtain an atomic function
  //    'func'.  The '111' tells tiny_ad that we need a derivative
  //    with respect to every variable so that we can test the gradient.
  TMB_BIND_ATOMIC(func, 111, eval(x[0], x[1], x[2]))
  // 3. Create a more user-friendly version ('func' takes vector
  //    arguments and there's a final invisible argument that
  //    corresponds to the derivative order)
  template<class Type>
  Type IntegralFunctionTest(Type x, Type y, Type theta) {
    vector<Type> args(4); // Last index reserved for derivative order
    args << x, y, theta, 0;
    return LocalCop::func(CppAD::vector<Type>(args))[0];
  }



  /// Calculate Gaussian copula partial derivative with respect to u1.
  ///
  /// @param[in] u1 First uniform variable.
  /// @param[in] u2 Second uniform variable. 
  /// @param[in] theta Parameter of the Gaussian copula with the range $(-1, 1)$.
  /// @param give_log Whether or not to return on the log scale. 
  ///
  /// @return Value of the h-function.  
  template <class Type>
  Type hgaussian(Type u1, Type u2, Type theta, int give_log=0) {
    Type z1 = qnorm(u1);
    Type z2 = qnorm(u2);
    Type determinant = Type(1.0) - theta * theta;
    Type ans = pnorm((z2 - theta * z1) / sqrt(determinant));
    if(give_log) return log(ans); else return ans;
  }
  VECTORIZE4_ttti(hgaussian)
      
  /// Calculate Gaussian copula PDF.
  ///
  /// @param[in] u1 First uniform variable.
  /// @param[in] u2 Second uniform variable. 
  /// @param[in] theta Parameter of the Gaussian copula with the range $(-1, 1)$.
  /// @param give_log Whether or not to return on the log scale. 
  ///
  /// @return Value of the copula PDF. 
  template <class Type>
  Type dgaussian(Type u1, Type u2, Type theta, int give_log=0) {
    // normal quantiles
    Type z1 = qnorm(u1);
    Type z2 = qnorm(u2);
    Type det = 1.0 - theta*theta;
    Type ans = theta*theta * (z1*z1 + z2*z2) - 2.0*theta * z1*z2;
    ans = -.5 * (ans / det + log(det));
    if(give_log) return ans; else return exp(ans);
  }
  VECTORIZE4_ttti(dgaussian)

} // end namespace LocalCop

#endif // LOCALCOP_GAUSSIAN_HPP
