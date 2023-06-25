#include <cmath>

#ifndef GAUSS_INTEGRATION_H
#define GAUSS_INTEGRATION_H

#ifndef EPS_GAU
// EPS is the relative precision.
#define EPS_GAU 1.0e-14
#endif

void gauleg(int npoints, double xst, double xend, double *x, double *w){
	double z1,z,xm,xl,pp,p3,p2,p1;
	int m = (npoints + 1)/2;
	xm = 0.5*(xend + xst);
	xl = 0.5*(xend - xst);

	for(int i = 0; i < m; i++){
		z = std::cos(M_PI*(i + 0.75)/(npoints + 0.5));
		do{
			p1=1.0;
			p2=0.0;
			for(int j = 0; j < npoints; j++){
				p3=p2;
				p2=p1;
				p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
			}
			pp = npoints*(z*p1 - p2)/(z*z - 1.0);
			z1 = z;
			z = z1 - p1/pp;
		}
		while(std::fabs(z-z1) > EPS_GAU);
		x[i] = xm - xl*z;
		x[npoints-1-i] = xm + xl*z;
		w[i] = 2.0*xl/((1.0 - z*z)*pp*pp);
		w[npoints-1-i] = w[i];
	}
}

#endif