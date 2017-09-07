#ifndef qqmatrix2_hh
#define qqmatrix2_hh

#include <iostream>
#include <complex>



using namespace std;

double Qqmatrix(double p1, double p2, double g, double theta, double j3, double j4, double p3, double E3, double E4, double p4, double m)
{
	 
	 double Matrix =(0.1111111111111111049*pow(g,4)*(1.*pow(E3,3)*E4*j3*j4*p3 - 1.*E3*E4*j3*j4*pow(m,2)*p3 - 1.*pow(E3,2)*E4*j3*j4*p1*p3 - 
       1.*pow(E3,2)*E4*j3*j4*p2*p3 - 1.*E3*E4*j3*j4*pow(p3,3) + 1.*E4*j3*j4*p1*pow(p3,3) + 1.*E4*j3*j4*p2*pow(p3,3) - 
       2.*pow(E3,2)*pow(m,2)*p4 + 2.*pow(m,4)*p4 + 2.*E3*pow(m,2)*p1*p4 + 2.*E3*pow(m,2)*p2*p4 + 1.*pow(m,2)*pow(p3,2)*p4 + 
       1.*pow(E3,2)*j3*j4*p3*pow(p4,2) - 2.*j3*j4*pow(m,2)*p3*pow(p4,2) - 1.*j3*j4*pow(p3,3)*pow(p4,2) + 
       (pow(E3,3)*E4*j3*j4*(-2.*p1 + 2.*p2) + E3*E4*j3*j4*(2.*p1 - 2.*p2)*(pow(m,2) + pow(p3,2)) + 
          pow(E3,2)*j3*j4*(E4*(2.*pow(p1,2) - 2.*pow(p2,2)) + (-2.*p1 + 2.*p2)*pow(p4,2)) + 
          p3*(E4*j3*j4*(-2.*pow(p1,2) + 2.*pow(p2,2))*p3 + p4*(pow(m,2)*(-2.*p1 + 2.*p2) + j3*j4*(2.*p1 - 2.*p2)*p3*p4)))*cos(theta) + 
       p3*(1.*pow(E3,3)*E4*j3*j4 + E3*E4*j3*j4*(-1.*pow(m,2) - 1.*pow(p3,2)) + pow(E3,2)*j3*j4*(-1.*E4*p1 - 1.*E4*p2 + 1.*pow(p4,2)) + 
          p3*(E4*j3*j4*(1.*p1 + 1.*p2)*p3 + p4*(1.*pow(m,2) - 1.*j3*j4*p3*p4)))*cos(2.*theta)))/(pow(m,2)*p1*p2*p4);
    
            
	return Matrix;
}

#endif
