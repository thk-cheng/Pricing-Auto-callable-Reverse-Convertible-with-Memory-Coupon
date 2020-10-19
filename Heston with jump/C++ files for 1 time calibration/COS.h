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

complex<double> D (double omega,double lamda,double eta,double rho);    //define the function D
complex<double> G (double omega,double lamda,double eta,double rho);    //define the function G
//define the characteristic function of the log-asset price under the Heston model
complex<double> phi_hes_j (double r,double q, double omega,double lamda,double nu_bar,double u_0,double eta,double rho,double freq, double mu_j, double var_j);
double psi(double c,double d,double a,double b,double k);           //calculate the cosine series coefficients of g(y)=1
double chi(double c,double d,double a,double b,double k);           //calculate the cosine series coefficients of g(y)=exp(y)
double U(double k,int C_P,double a,double b);         //calculate the payoff function of the option in form of sin and cos
double a(double r,double q,double T,double lamda,double nu_bar,double u_0,double eta,double rho,double freq,double mu_j,double var_j);    //calculate the lower bound of the log-assest price
double b(double r,double q,double T,double lamda,double nu_bar,double u_0,double eta,double rho,double freq,double mu_j,double var_j);    //calculate the upper bound of teh log-assear price
//define the COS pricing function
double COS(double T,double lamda,double nu_bar,double u_0,double eta,double rho,double freq, double mu_j, double var_j, double r,double q,double S_0,double K,int C_P,int N);
#endif
