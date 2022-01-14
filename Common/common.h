#ifndef SWITCH
#define SWITCH 1 // v mode 1 pocita s boost multiprecision, v mode 0 s built-in double presnostou
#endif

#include <cmath>
#include <fstream>
#include <omp.h>

#if SWITCH == 1
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/eigen.hpp>
#include <boost/multiprecision/mpc.hpp>
#include <boost/multiprecision/mpfr.hpp>
#endif

#include <complex>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <iostream>
#include <vector>

#if SWITCH == 1
#define COMPLEX_TYPE boost::multiprecision::number<boost::multiprecision::backends::mpc_complex_backend<50>>
#define FLOAT_TYPE boost::multiprecision::number<boost::multiprecision::backends::mpfr_float_backend<50>>
#define PI boost::math::constants::pi<FLOAT_TYPE>()
using namespace boost::multiprecision;
#endif

#if SWITCH == 0
#define COMPLEX_TYPE complex<double>
#define FLOAT_TYPE double
#define PI M_PI
#endif

const COMPLEX_TYPE i(0, 1);

Eigen::Matrix<COMPLEX_TYPE, 2, 2> transfermatrix(FLOAT_TYPE, FLOAT_TYPE,
                                                 FLOAT_TYPE, FLOAT_TYPE, FLOAT_TYPE);
Eigen::Matrix<COMPLEX_TYPE, 2, 2> intermatrix(FLOAT_TYPE, FLOAT_TYPE,
                                              FLOAT_TYPE, FLOAT_TYPE, FLOAT_TYPE);

Eigen::Matrix<COMPLEX_TYPE, 2, 2> transfermatrix(FLOAT_TYPE permittivity1,
                                                 FLOAT_TYPE permittivity2,
                                                 FLOAT_TYPE theta,
                                                 FLOAT_TYPE omega,
                                                 FLOAT_TYPE eps_one)
{

      FLOAT_TYPE inter_1 = permittivity1 - eps_one * sin(theta) * sin(theta);
      COMPLEX_TYPE k_1 = omega * sqrt(COMPLEX_TYPE(inter_1, 0));
      FLOAT_TYPE inter_2 = permittivity2 - eps_one * sin(theta) * sin(theta);
      COMPLEX_TYPE k_2 = omega * sqrt(COMPLEX_TYPE(inter_2, 0));
      COMPLEX_TYPE chi_te = k_1 / k_2;
      Eigen::Matrix<COMPLEX_TYPE, 2, 2> Out;
      Out << 1. + chi_te, 1. - chi_te, 1. - chi_te, 1. + chi_te;
      Out /= 2.;
      return Out;
}

Eigen::Matrix<COMPLEX_TYPE, 2, 2> intermatrix(FLOAT_TYPE permittivity,
                                              FLOAT_TYPE width,
                                              FLOAT_TYPE theta,
                                              FLOAT_TYPE omega,
                                              FLOAT_TYPE eps_one)
{
      FLOAT_TYPE inter = permittivity - eps_one * sin(theta) * sin(theta);
      COMPLEX_TYPE k = omega * sqrt(COMPLEX_TYPE(inter, 0));
      COMPLEX_TYPE phi = i * k * width;
      Eigen::Matrix<COMPLEX_TYPE, 2, 2> Out;
      Out << exp(phi), 0, 0, exp(-phi);
      return Out;
}
