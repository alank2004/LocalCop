// Generated by TMBtools: do not edit by hand

#define TMB_LIB_INIT R_init_LocalCop_TMBExports
#include <TMB.hpp>
#include "dclayton.hpp"
#include "dfrank.hpp"
#include "dgaussian.hpp"
#include "dgumbel.hpp"
#include "dstudent.hpp"
#include "hclayton.hpp"
#include "hfrank.hpp"
#include "hgaussian.hpp"
#include "hgumbel.hpp"
#include "hstudent.hpp"
#include "integral_function_test.hpp"
#include "LocalLikelihood.hpp"
#include "pclayton.hpp"
#include "pfrank.hpp"
#include "pgumbel.hpp"
#include "pt.hpp"
#include "qt.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model);
  if(model == "dclayton") {
    return dclayton(this);
  } else if(model == "dfrank") {
    return dfrank(this);
  } else if(model == "dgaussian") {
    return dgaussian(this);
  } else if(model == "dgumbel") {
    return dgumbel(this);
  } else if(model == "dstudent") {
    return dstudent(this);
  } else if(model == "hclayton") {
    return hclayton(this);
  } else if(model == "hfrank") {
    return hfrank(this);
  } else if(model == "hgaussian") {
    return hgaussian(this);
  } else if(model == "hgumbel") {
    return hgumbel(this);
  } else if(model == "hstudent") {
    return hstudent(this);
  } else if(model == "integral_function_test") {
    return integral_function_test(this);
  } else if(model == "LocalLikelihood") {
    return LocalLikelihood(this);
  } else if(model == "pclayton") {
    return pclayton(this);
  } else if(model == "pfrank") {
    return pfrank(this);
  } else if(model == "pgumbel") {
    return pgumbel(this);
  } else if(model == "pt") {
    return pt(this);
  } else if(model == "qt") {
    return qt(this);
  } else {
    Rf_error("Unknown model.");
  }
  return 0;
}
