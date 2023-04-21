#include <iostream>
#include <iomanip>
#include <tests.h>
#include <eigen_solver.h>
#include <bose_star_perturbations.h>

void save_grid(const char *file_name, const int nfunc, const int npoints, const double *xarr, const double *yarr){
	std::ofstream out;
	out.open(file_name);

	out << std::setprecision(16);
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

int main(int argc, char const *argv[]){
	std::cout << "Running tests\n" << std::endl;
	eigen_solver_test(argc, argv);
	stag_diff_matrix_test(argc, argv);
	spectrum_solver_test(argc, argv);

	std::cout << "Finding Boson Star perturbation frequencies\n" << std::endl;
	int n = 20;
	for(int i = 1; i < argc; i++){
		if(strcmp(argv[i], "-n") == 0){
			n = atoi(argv[i+1]);
			i++;
		}
		else{
			std::cout << "Unknown key received" << std::endl;
			return 0;
		}
	}

	int d = 3;
	double l = 1.0;
	double kappa = 1.0;
	double mass = 0.0;
	double el = 1.0;
	double omega = 2.7057244879610716;
	double phi0 = 0.5;
	double xst = 0.001;
	double xend = 1.4;
	double delta0 = 0.375180752629896;
	double a0 = 0.1908970702516414;

	// int d = 4;
	// double l = 1.0;
	// double kappa = 1.0;
	// double mass = 0.0;
	// double el = 1.0;
	// double omega = 3.734428385858555;
	// double phi0 = 0.5;
	// double xst = 0.001;
	// double xend = 1.4;
	// double delta0 = 0.33011092330920855;
	// double a0 = 0.09274485030170596;

	Boson_Star_Perturbations star_perturbations(d, l, omega, phi0, delta0, a0, xst, xend);
	star_perturbations.update_lagrangian_parameters(kappa, mass, el);

	star_perturbations.evaluate_unperturbed_fields(n);
	// save_grid("test_save.csv", 6, n, star_perturbations.star_profile.xarr, star_perturbations.star_profile.yarr);

	SpectrumSolver solver(&star_perturbations, n);
	solver.compute_spectrum(true);
	solver.summary();
	// solver.save_real_eigenvalues_csv("eigenvectors.csv");

	std::vector<int> select = {1, 2, 3, 4, 5};
	std::vector<double> precision(5);
	solver.eigenvectors_precision(select, precision.data());

	for(int i = 0; i < 5; ++i){
		std::cout << "Estimated spectral precision for eigenvalue " << select[i] << " is " << precision[i] << std::endl;

		Eigen::VectorXd res = solver.matr_a * solver.real_eigenvectors[select[i]] -
								solver.real_freq[select[i]]*solver.matr_b * solver.real_eigenvectors[select[i]];

		std::cout << "Eigenvector error is: " << res.cwiseAbs().maxCoeff() << '\n' << std::endl;
	}

	return 0;
}