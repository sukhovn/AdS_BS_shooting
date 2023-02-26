#include <cstring>
#include <star_equations.h>
#include <star_shooter.h>

// params contain solution parameters
// params[0] = l
// params[1] = kappa
// params[2] = mass
// params[3] = phi0

extern "C" int bose_shooter(int d, double *params, double xst, double xend, double omega_st, double omega_end, double *omega, double *delta0){
	Boson_Star_Solver solver(d, params[0], 1.0, params[3], xst, xend);
	solver.bstar.kappa = params[1];
	solver.bstar.mass = params[2];
	
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
		else{
			std::cout << "Shooting exited with unknown exception" << std::endl;
			return -1;
		}
	}

	*omega = solver.bstar.omega;
	*delta0 = solver.delta0;
	return 0;
}

extern "C" int get_star_profile(int d, double *params, double xst, double xend, double omega, double delta0, int npoints, double *xarr, double *yarr){
	Boson_Star_Solver solver(d, params[0], omega, params[3], xst, xend);
	solver.bstar.kappa = params[1];
	solver.bstar.mass = params[2];
	solver.delta0 = delta0;
	
	try{
		solver.integrate_bose_star(npoints, xarr, yarr);
	}
	catch(const char *error){
		if(strcmp(error, "stepsize underflow in RKDP_stepper") == 0){
			std::cout << "Runge-Kutta solver cannot integrate in this region" << std::endl;
			return 1;
		}
		else{
			std::cout << "Shooting exited with unknown exception" << std::endl;
			return -1;
		}
	}

	
	return 0;
}