#include <iostream>
#include <cmath>
#include <iomanip>

#include <star_equations.h>

#ifndef STAR_SHOOTER_H
#define STAR_SHOOTER_H

#define SHOOT_TRIES 30

class Boson_Star_Solver{
public:
	Boson_Star_AdS bstar;

	double xst;
	double xend;

	double phi0;
	double delta0;

	Boson_Star_Solver(int d, double l, double omega, double phi0i, double xsti, double xendi): 
						bstar(d, l, omega), xst(xsti), xend(xendi),
						phi0(phi0i), delta0(0.0) {};

	void dphi_difference_shot(double &d_phi, double &d_delta){
		double *xint = new double[2];
		xint[0] = xst; xint[1] = xend;
		double *yint = new double[8];
		yint[0] = delta0;
 		yint[1] = phi0;

 		bstar.integrate_grid(1, xint, yint);

 		Boson_Star_AdS::Boundary_Asymptotic asym_bound(bstar);
 		d_phi = yint[7] - asym_bound.evaluate_dphi();
 		d_delta = asym_bound.del0 - bstar.y[1];

 		delete [] xint;
 		delete [] yint;

 		return;
	}

	void shoot_for_omega(double omega_st, double omega_end){
		double dphi_st, dphi_end, delta_diff;
		bstar.omega = omega_st;
		dphi_difference_shot(dphi_st, delta_diff);
		bstar.omega = omega_end;
		dphi_difference_shot(dphi_end, delta_diff);
		if(dphi_st*dphi_end > 0.0){
			std::cout << "Improper omega interval chosen" << std::endl;
			std::cout << "phi' difference for omega " << omega_st << " is " << dphi_st << std::endl;
			std::cout << "phi' difference for omega " << omega_end << " is " << dphi_end << std::endl;
			std::cout << "Differences should have opposite signs" << std::endl;
			throw("improper omega interval");
		}
		double dphi_try, omega_try;

		double *yasym = new double[4];
		
		//Shooting loop
		int at;
		for(at = 0; at < SHOOT_TRIES; ++at){
			omega_try = 0.5*(omega_st + omega_end);
			bstar.omega = omega_try;
			dphi_difference_shot(dphi_try, delta_diff);
			if(dphi_try*dphi_st > 0.0){
				dphi_st = dphi_try;
				omega_st = omega_try;
			}
			else{
				dphi_end = dphi_try;
				omega_end = omega_try;
			}

			//Dynamically adjusting delta0 while shooting for improved accuracy
			delta0 += delta_diff;
			omega_st *= std::exp(-delta_diff);
			omega_end *= std::exp(-delta_diff);

			if(std::fabs(dphi_try) < 1.0e-6)
				break;
		}

		delete [] yasym;

		if(at == SHOOT_TRIES){
			std::cout << std::setprecision(16);
			std::cout << "Final omega interval [" << omega_st << ", " << omega_end << "]" << std::endl;
			std::cout << "Final phi' differences [" << dphi_st << ", " << dphi_end << "]" << std::endl;
			std::cout << "Shooting didn't converge" << std::endl;
			throw("shooting didn't converge");
		}

		// delta0 = delta_diff;
		// bstar.omega *= std::exp(-delta_diff);
	
		return;
	}

	void integrate_bose_star(int npoints, double *xarr, double *yarr){
		double h = M_PI_2/(1.0*(npoints-1)), xtmp = 0.0;
		for(int i = 0; i < npoints; ++i, xtmp += h) xarr[i] = xtmp;
		for(int i = 0; i < 4*npoints; ++i) yarr[i] = 0.0;

		// Determining the start of integration domain
		int ist = 0;
		double *xint = xarr;
		double *yint = yarr;
		for(xtmp = 0.0; xtmp < xst; xtmp += h){
			++ist;
		}
		xtmp -= h;
		--ist;
		xint += ist;
		yint += 4*ist;

		Boson_Star_AdS::Origin_Asymptotic asym_orig(bstar, delta0, phi0);

		// Determining the end of integration domain
		int iend = 0;
		for(; xtmp < xend; xtmp += h){
			++iend;
		}
		xint[iend] = xend;

		// Integrating
		yint[0] = delta0;
 		yint[1] = phi0;
 		xint[0] = xst;

 		bstar.integrate_grid(iend, xint, yint);

 		// Filling values before xst with asymptotics
 		xtmp = 0.0;
 		yint = yarr;
 		for(int i = 0; i <= ist; ++i){
 			asym_orig.evaluate(xtmp, yint);
 			yint += 4;
 			xtmp += h;
 		}
 		xarr[ist] = xtmp-h;

 		// Filling values after xend with asymptotics
 		Boson_Star_AdS::Boundary_Asymptotic asym_bound(bstar);

 		xarr[ist+iend] = xtmp = h*(ist+iend);
 		yint = yarr + 4*(ist + iend);
 		for(int i = ist + iend; i < npoints-1; ++i){
 			asym_bound.evaluate(xtmp, yint);
 			xtmp += h;
 			yint += 4;
 		}
 		yint[0] = 1.0;
 		yint[1] = yint[2] = yint[3] = 0.0;
	}
};

#endif