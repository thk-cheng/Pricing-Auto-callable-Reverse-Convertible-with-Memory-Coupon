#ifndef _COS_h_
#define _COS_h_

#include <stdio.h>
#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include <list>
#define pi 3.14159265358979

using namespace std;

int factorial(int n); //calculate the factorial
double norm_s_cdf(double x); //calculate the standard normal cdf 
double m_n(double r, double sigma_s, double freq, double mu_j, double var_j, double T, int N); //component of closed-form call price
double d1(double S_0,double K,double m_n,double T,double sigma_s,double var_j,int N); //component of closed-form call price
double d2(double d1,double sigma_s,double T,double var_j,int N); //component of closed-form call price
double Jump_option_price(double T,double sigma_s,double freq, double mu_j, double var_j, double r,double q,double S_0,double K,int N); //Calculate closed-form call price
#endif
