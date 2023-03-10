#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <iomanip>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <ode_solver.h>
#include <star_equations.h>
#include <star_shooter.h>

void save_grid(const char *file_name, const int nfunc, const int npoints, const double *xarr, const double *yarr){
	std::ofstream out;
	out.open(file_name);

	out << "x";
	for(int i = 0; i < nfunc; i++) out << "," << "y" << i;
	out << std::endl;

	const double *ycurr = yarr;
	for(int i = 0; i < npoints; i++){
		out << xarr[i];
		for(int i = 0; i < nfunc; i++) out << "," << ycurr[i];
		out << std::endl;
		ycurr += nfunc;
	}

	out.close();
	return;
}

// This routine runs ordinary 4th order Runge-Kutta solver

int main(int argc, char const *argv[]){
	struct stat sttus = {0};
	if (stat("Data", &sttus) == -1) {
		mkdir("Data", 0700);
	}

	double phi0 = 0.1, omega = 3.3;
	double el = 0.1;

	for(int i = 1; i < argc; i++){
		if(strcmp(argv[i], "-phi0") == 0){
			phi0 = atof(argv[i+1]);
			i++;
		}
		else if(strcmp(argv[i], "-omega") == 0){
			omega = atof(argv[i+1]);
			i++;
		}
		else if(strcmp(argv[i], "-el") == 0){
			el = atof(argv[i+1]);
			i++;
		}
		else{
			std::cout << "Unknown key received" << std::endl;
			return 0;
		}
	}

	int npoints = 501;
	double *xarr = new double[npoints];
	double *yarr = new double[6*npoints];

	Boson_Star_Solver solver(3, 1.0, omega, phi0, 0.001, 1.3);
	solver.bstar.el = el;

	try{
		solver.shoot_for_omega(3.0, 4.0);
	}
	catch(const char *error){
		printf("Error code: %s\n", error);
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

	std::cout << "Omega found: " << solver.bstar.omega << std::endl;
	
	solver.integrate_bose_star(npoints, xarr, yarr);
	save_grid("Data/solution.dat", 6, npoints, xarr, yarr);

	return 0;
}