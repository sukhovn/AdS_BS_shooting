#include <cstring>
#include <star_equations.h>
#include <star_shooter.h>

// params contain solution parameters
// params[0] = l
// params[1] = kappa
// params[2] = mass
// params[3] = el
// params[4] = phi0
// params[5] = omega
// params[6] = delta0
// params[7] = a0

extern "C" int bose_shooter(int d, double *params, double xst, double xend, double omega_st, double omega_end, double *params_out){
	if(d < 2){
		std::cout << "Only dimensions two or higher are allowed" << std::endl;
		return 1;
	}

	Boson_Star_Solver solver(d, params[0], 1.0, params[4], xst, xend);
	solver.bstar.kappa = params[1];
	solver.bstar.mass = params[2];
	solver.bstar.el = params[3];
	
	try{
		solver.shoot_for_omega(omega_st, omega_end);
	}
	catch(const char *error){
		if(strcmp(error, "stepsize underflow in RKDP_stepper") == 0){
			std::cout << "Runge-Kutta solver cannot integrate in this region" << std::endl;
			return 1;
		}
		else if(strcmp(error, "improper omega interval") == 0){
			return 1;
		}
		else if(strcmp(error, "shooting didn't converge") == 0){
			return 1;
		}
		else if(strcmp(error, "2D charged boson star not permitted") == 0){
			std::cout << "2D charged Boson Star has a singular electric potential" << std::endl;
			return 1;
		}
		else{
			std::cout << "Shooting exited with unknown exception" << std::endl;
			return -1;
		}
	}

	params_out[0] = solver.bstar.omega;
	params_out[1] = solver.delta0;
	params_out[2] = solver.a0;
	params_out[3] = solver.star_mass;
	params_out[4] = solver.star_charge;

	return 0;
}


extern "C" int get_difference_with_asymptotic(int d, double *params, double xst, double xend, double *dy_diff){
	if(d < 2){
		std::cout << "Only dimensions two or higher are allowed" << std::endl;
		return 1;
	}

	Boson_Star_Solver solver(d, params[0], params[5], params[4], xst, xend);
	solver.bstar.kappa = params[1];
	solver.bstar.mass = params[2];
	solver.bstar.el = params[3];
	solver.delta0 = params[6];
	solver.a0 = params[7];

	double d_delta, d_a0;

	try{
		solver.dphi_difference_shot(*dy_diff, d_delta, d_a0);
	}
	catch(const char *error){
		if(strcmp(error, "stepsize underflow in RKDP_stepper") == 0){
			std::cout << "Runge-Kutta solver cannot integrate in this region" << std::endl;
			return 1;
		}
		else if(strcmp(error, "2D charged boson star not permitted") == 0){
			std::cout << "2D charged Boson Star has a singular electric potential" << std::endl;
			return 1;
		}
		else{
			std::cout << "Shooting exited with unknown exception" << std::endl;
			return -1;
		}
	}

	return 0;
}


extern "C" int get_star_profile(int d, double *params, double xst, double xend, int npoints, double *xarr, double *yarr){
	if(d < 2){
		std::cout << "Only dimensions two or higher are allowed" << std::endl;
		return 1;
	}

	Boson_Star_Solver solver(d, params[0], params[5], params[4], xst, xend);
	solver.bstar.kappa = params[1];
	solver.bstar.mass = params[2];
	solver.bstar.el = params[3];
	solver.delta0 = params[6];
	solver.a0 = params[7];
	
	try{
		solver.integrate_bose_star(npoints, xarr, yarr);
	}
	catch(const char *error){
		if(strcmp(error, "stepsize underflow in RKDP_stepper") == 0){
			std::cout << "Runge-Kutta solver cannot integrate in this region" << std::endl;
			return 1;
		}
		else if(strcmp(error, "2D charged boson star not permitted") == 0){
			std::cout << "2D charged Boson Star has a singular electric potential" << std::endl;
			return 1;
		}
		else{
			std::cout << "Shooting exited with unknown exception" << std::endl;
			return -1;
		}
	}
	
	return 0;
}