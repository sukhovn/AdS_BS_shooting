#include <eigen_solver.h>
#include <bose_star_perturbations.h>

// params contain solution parameters
// params[0] = l
// params[1] = kappa
// params[2] = mass
// params[3] = el
// params[4] = phi0
// params[5] = omega
// params[6] = delta0
// params[7] = a0

extern "C" int bose_mode_frequency(int d, double *params, int npoints, double xst, double xend, int nmodes, double *frequencies, double *precision){
	if(d < 3){
		std::cout << "Only dimensions three or higher are allowed" << std::endl;
		return 1;
	}

	if(std::fabs(params[3]) < 1.0e-5){
		std::cout << "The code requires the gauge coupling to be nonzero" << std::endl;
		return 1;
	}

	Boson_Star_Perturbations star_perturbations(d, params[0], params[5], params[4], params[6], params[7], xst, xend);
	star_perturbations.update_lagrangian_parameters(params[1], params[2], params[3]);

	try{
		star_perturbations.evaluate_unperturbed_fields(npoints);
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
			std::cout << "Mode  exited with unknown exception" << std::endl;
			return -1;
		}
	}

	SpectrumSolver solver(&star_perturbations, npoints);
	solver.compute_spectrum(true);

	std::vector<int> modes_select;
	int indx = 0;
	for(int i = 0; i < solver.real_freq.size(); ++i){
		if(std::fabs(solver.real_freq[i]) < 1.0e-5)
			continue;
		modes_select.push_back(i);
		frequencies[indx++] = solver.real_freq[i];
		if(indx == nmodes)
			break;
	}

	if(indx < nmodes){
		std::cout << "Warning! Not enough modes obtained, consider increasing resolution." << std::endl;
		return 1;
	}

	solver.eigenvectors_precision(modes_select, precision);

	return 0;
}


extern "C" int bose_mode_profiles(int d, double *params, int npoints, double xst, double xend, int nmodes, double *frequencies, double *profiles){
	if(d < 3){
		std::cout << "Only dimensions three or higher are allowed" << std::endl;
		return 1;
	}

	if(std::fabs(params[3]) < 1.0e-5){
		std::cout << "The code requires the gauge coupling to be nonzero" << std::endl;
		return 1;
	}

	Boson_Star_Perturbations star_perturbations(d, params[0], params[5], params[4], params[6], params[7], xst, xend);
	star_perturbations.update_lagrangian_parameters(params[1], params[2], params[3]);

	try{
		star_perturbations.evaluate_unperturbed_fields(npoints);
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
			std::cout << "Mode  exited with unknown exception" << std::endl;
			return -1;
		}
	}

	SpectrumSolver solver(&star_perturbations, npoints);
	solver.compute_spectrum(true);

	std::vector<int> modes_select;
	int indx = 0;
	for(int i = 0; i < solver.real_freq.size(); ++i){
		if(std::fabs(solver.real_freq[i]) < 1.0e-5)
			continue;
		modes_select.push_back(i);
		frequencies[indx++] = solver.real_freq[i];
		if(indx == nmodes)
			break;
	}

	if(indx < nmodes){
		std::cout << "Warning! Not enough modes obtained, consider increasing resolution." << std::endl;
		return 1;
	}

	indx = 0;
	for(int el : modes_select){
		for(int i = 0; i < solver.npoints; ++i){
			double cos_x = std::cos(solver.diff_matrix.x[i]);
			profiles[indx++] = std::pow(cos_x, 3.0)*solver.real_eigenvectors[el](3*i);
			profiles[indx++] = solver.real_eigenvectors[el](3*i+1);
			profiles[indx++] = std::pow(cos_x, 2.0)*solver.real_eigenvectors[el](3*i+2);
		}
	}

	return 0;
}