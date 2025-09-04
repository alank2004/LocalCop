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

  /*
    This class evaluates the integral of the bivariate normal distribution via the direct method.
    Since we already have pnorm which is fast and accurate, we can avoid doing 2D integration by
    conditioning on one of the two variables, in our case on X1 and its associated variable b1.
  */
  template<class Float>
  struct BVNIntegrand {
    typedef Float Scalar; // Required by integrate
    Float b1, b2, rho;         // Parameters 
    // Evaluate conditional CDF
    Float operator() (Float x) {
      Float loc = rho * x;
      Float scale = sqrt(1 - rho * rho);

      // Replace Float(0.0) by adding mu to parameters above to control mean
      // Replace the first Float(1.0) by adding sigma1 to parameters above to control the first s.d.
      // Replace the second Float(1.0) by adding sigma2 to parameters above to control the second s.d.
      Float ans = pnorm((b2 - loc) / scale, Float(0.0), Float(1.0)) *
                  dnorm(x, Float(0.0), Float(1.0), false);


      // Float ans = pnorm((b2 - loc) / scale, Float(0.0), Float(1.0), false, true) *
      //             dnorm(x, Float(0.0), Float(1.0), false);
      return ans;
    }
    // Integrate conditional CDF into the CDF
    Float integrate() {
      using gauss_kronrod::integrate;
      Float ans =
        integrate(*this, Float(-INFINITY), b1);
      return ans;
    }
  };


  // An externally available integration evaluator
  template<class Float>
  Float pbvn(Float b1, Float b2, Float rho) {
    BVNIntegrand<Float> f = {b1, b2, rho};
    return f.integrate();
  }

  VECTORIZE3_ttt(pbvn)

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
