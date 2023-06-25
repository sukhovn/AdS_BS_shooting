#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>

#include <star_equations.h>
#include <gauss_legendre.h>

#ifndef STAR_SHOOTER_H
#define STAR_SHOOTER_H

#define SHOOT_TRIES 40
#define INTEGRATION_POINTS 300

class Boson_Star_Solver{
	double xst;
	double xend;

public:
	Boson_Star_AdS bstar;

	double phi0;
	double delta0;
	double a0;

	double star_mass;
	double star_charge;

	//Mass integration part
	std::vector<double> x_origin;
	std::vector<double> x_integral;
	std::vector<double> x_bound;

	std::vector<double> weight_gauss;

	Boson_Star_Solver(int d, double l, double omega, double phi0i, double xsti, double xendi): 
						bstar(d, l, omega), xst(xsti), xend(xendi),
						phi0(phi0i), delta0(0.0), a0(0.0),
						weight_gauss(INTEGRATION_POINTS)
	{
		if(d < 2){
			throw("dimension less than two");
		}

		std::vector<double> x_gauss(INTEGRATION_POINTS);
		gauleg(INTEGRATION_POINTS, 0.0, M_PI_2, x_gauss.data(), weight_gauss.data());

		for(double point : x_gauss){
			if(point <= xst) x_origin.push_back(point);
			else if(point >= xend) x_bound.push_back(point);
			else x_integral.push_back(point);
		}
	};

	void dphi_difference_shot(double &d_phi, double &d_delta, double &d_a0){
		double *ycurr = new double[bstar.nfunc];

		double mass_integral = 0.0;
		std::vector<double>::iterator weight = weight_gauss.begin();

		int npoints = 1;

		bstar.dense = true;

		//Setting up asymptotic
		Boson_Star_AdS::Origin_Asymptotic asym(bstar, delta0, phi0, a0);

		for(double xcoord : x_origin){
			asym.evaluate(xcoord, ycurr);
			mass_integral += (*weight) * bstar.energy_density(xcoord, ycurr);
			++weight;
		}

		//Setting up integration
		bstar.x = xst;
		asym.evaluate(xst, bstar.y);
		bstar.derivatives(xst, bstar.y, bstar.dy);

		//Integrating
		std::vector<double>::iterator point = x_integral.begin();
		while(bstar.x + bstar.h < xend){
			bstar.step();
			while(point != x_integral.end() && bstar.x > *point){
				for(int i = 0; i < bstar.nfunc; ++i)
					ycurr[i] = bstar.dense_out(i, *point);
				mass_integral += (*weight) * bstar.energy_density(*point, ycurr);
				++weight;
				++point;
			}
		}

		bstar.fixed_step(xend - bstar.x);
		for(; point != x_integral.end(); ++point){
			for(int i = 0; i < bstar.nfunc; i++)
				ycurr[i] = bstar.dense_out(i, *point);
			mass_integral += (*weight) * bstar.energy_density(*point, ycurr);
			++weight;
		}

 		Boson_Star_AdS::Boundary_Asymptotic asym_bound(bstar);
		d_phi = bstar.y[3] - asym_bound.evaluate_dphi();
		d_delta = asym_bound.del0;
		d_a0 = asym_bound.a0;

		for(double xcoord : x_bound){
			mass_integral += (*weight) * asym_bound.evaluate_energy_density(xcoord);
			++weight;
		}

		mass_integral *= bstar.kappa*std::pow(bstar.l, bstar.d)/(bstar.d - 1.0);

		// Finding star mass and charge
		if(bstar.d == 2){
			if(bstar.kappa == 0.0)
				star_mass = 0.0;
			star_mass = -0.5*asym_bound.A0/std::pow(asym_bound.rho, bstar.d);
		}
		else
		star_mass = mass_integral;

		// if(bstar.kappa == 0.0)
		// 	star_mass = 0.0;
		// else{
		// 	double lpow = 1.0;
		// 	if(bstar.d > 2)
		// 		lpow = std::pow(bstar.l, bstar.d-2);
		// 	star_mass = -0.5*asym_bound.A0*lpow/std::pow(asym_bound.rho, bstar.d);
		// }

		if(bstar.el == 0)
			star_charge = 0.0;
		else
			star_charge = (bstar.d-2)*asym_bound.a1/bstar.l/std::pow(asym_bound.rho, bstar.d-2)*std::exp(asym_bound.del0)/bstar.el;

		delete [] ycurr;

 		return;
	}

	void shoot_for_omega(double omega_st, double omega_end){
		double dphi_st, dphi_end, delta_diff, a0_diff;
		bstar.omega = omega_st;
		dphi_difference_shot(dphi_st, delta_diff, a0_diff);
		bstar.omega = omega_end;
		dphi_difference_shot(dphi_end, delta_diff, a0_diff);
		if(dphi_st*dphi_end > 0.0){
			std::cout << "Improper omega interval chosen" << std::endl;
			std::cout << "phi' difference for omega " << omega_st << " is " << dphi_st << std::endl;
			std::cout << "phi' difference for omega " << omega_end << " is " << dphi_end << std::endl;
			std::cout << "Differences should have opposite signs" << std::endl;
			throw("improper omega interval");
		}
		double dphi_try, omega_try;
		
		// Shooting loop
		int at;
		for(at = 0; at < SHOOT_TRIES; ++at){
			omega_try = 0.5*(omega_st + omega_end);
			bstar.omega = omega_try;
			dphi_difference_shot(dphi_try, delta_diff, a0_diff);
			if(dphi_try*dphi_st > 0.0){
				dphi_st = dphi_try;
				omega_st = omega_try;
			}
			else{
				dphi_end = dphi_try;
				omega_end = omega_try;
			}

			if(std::fabs(dphi_try) < 1.0e-6)
				break;

			// Dynamically adjusting delta0 and a0
			// omega_st -= bstar.el*a0_diff;
			// omega_end -= bstar.el*a0_diff;
			// a0 -= a0_diff;
			// delta0 -= delta_diff;
			// omega_st *= std::exp(delta_diff);
			// omega_end *= std::exp(delta_diff);
			// a0 *= std::exp(delta_diff);
		}

		if(at == SHOOT_TRIES){
			std::cout << std::setprecision(16);
			std::cout << "Final omega interval [" << omega_st << ", " << omega_end << "]" << std::endl;
			std::cout << "Final phi' differences [" << dphi_st << ", " << dphi_end << "]" << std::endl;
			std::cout << "Shooting didn't converge" << std::endl;
			throw("shooting didn't converge");
		}

		// Adjusting delta0 and a0
		bstar.omega -= bstar.el*a0_diff;
		a0 -= a0_diff;
		delta0 -= delta_diff;
		bstar.omega *= std::exp(delta_diff);
		a0 *= std::exp(delta_diff);

		return;
	}

	void integrate_bose_star(int npoints, double *xarr, double *yarr){
		int nfunc_save = 7;
		double h = M_PI_2/(1.0*(npoints-1)), xtmp = 0.0;
		for(int i = 0; i < npoints; ++i, xtmp += h) xarr[i] = xtmp;
		for(int i = 0; i < nfunc_save*npoints; ++i) yarr[i] = 0.0;

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
		yint += nfunc_save*ist;

		Boson_Star_AdS::Origin_Asymptotic asym_orig(bstar, delta0, phi0, a0);

		// Determining the end of integration domain
		int iend = 0;
		for(; xtmp < xend; xtmp += h){
			++iend;
		}
		xint[iend] = xend;

		// Integrating
		bstar.dense = true;
		bstar.x = xst;
		asym_orig.evaluate(xst, bstar.y);
		bstar.derivatives(xst, bstar.y, bstar.dy);

		int st = 1;
		while(bstar.x + bstar.h < xend){
			bstar.step();
			while(bstar.x > xint[st]){
				yint += nfunc_save;
				for(int i = 0; i < bstar.nfunc; i++) yint[i] = bstar.dense_out(i, xint[st]);
				yint[6] = bstar.energy_density(xint[st], yint);
				st++;
			}
		}

		bstar.fixed_step(xint[iend] - bstar.x);
		for(; st <= iend; st++){
			yint += nfunc_save;
			for(int i = 0; i < bstar.nfunc; i++) yint[i] = bstar.dense_out(i, xint[st]);
			yint[6] = bstar.energy_density(xint[st], yint);
		}

 		// Filling values before xst with asymptotics
 		xtmp = 0.0;
 		yint = yarr;
 		for(int i = 0; i <= ist; ++i){
 			asym_orig.evaluate(xtmp, yint);
			yint[6] = bstar.energy_density(xtmp, yint);
			yint += nfunc_save;
 			xtmp += h;
 		}
 		xarr[ist] = xtmp-h;

 		// Filling values after xend with asymptotics
 		Boson_Star_AdS::Boundary_Asymptotic asym_bound(bstar);

 		xarr[ist+iend] = xtmp = h*(ist+iend);
		yint = yarr + nfunc_save*(ist + iend);
 		for(int i = ist + iend; i < npoints-1; ++i){
 			asym_bound.evaluate(xtmp, yint);
			yint[6] = asym_bound.evaluate_energy_density(xtmp);
 			xtmp += h;
			yint += nfunc_save;
 		}
		asym_bound.evaluate_boundary(yint);
		yint[6] = asym_bound.evaluate_energy_density(M_PI_2);
	}
};

#endif