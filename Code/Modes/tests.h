#include <iostream>
#include <iomanip>
#include <stag_diff_matrix.h>
#include <eigen_solver.h>

// Routines to tst EigenSolver class
void eigen_solver_test(int argc, char const *argv[]){
	std::cout << "Running generalized eigenvalue solver tests\n" << std::endl; 
	
	EigenSolver solver(1, 3);

	Eigen::MatrixXd matr_a(3,3), matr_b(3,3);

	matr_a << 2, 1, 0,
			  1, 3, 1,
			  0, 1, 2;

	matr_b << 1, 0, 0,
			  0, 0, 0,
			  0, 0, 1;	

	solver.update_a(0, 0, matr_a);
	solver.update_b(0, 0, matr_b);

	solver.compute_eigenvalues();
	solver.summary();

	// for(int i = 0; i < 3; ++i){
	// 	for(int j = 0; j < 3; ++j){
	// 		if(j) std::cout << ", ";
	// 		std::cout << solver.matr_b(i, j);
	// 	}
	// 	std::cout << std::endl;
	// }

	// Eigen::VectorXd res = matr_a*solver.real_eigenvectors[0] - solver.real_freq[0]*matr_b*solver.real_eigenvectors[0];
	// std::cout << res.cwiseAbs().maxCoeff() << std::endl;

	for(int i = 0; i < 70; ++i)
		std::cout << '-';
	std::cout << '\n' << std::endl;
}


// Routines to test Chebyshev Differentiation Staggered_Differentiation_Matrix class
void stag_diff_matrix_test(int argc, char const *argv[]){
	std::cout << "Running Chebyshev differentiation matrix for staggered grid tests\n" << std::endl;

	int n = 20;

	for(int i = 1; i < argc; i++){
		if(strcmp(argv[i], "-n") == 0){
			n = atoi(argv[i+1]);
			i++;
		}
		// else{
		// 	std::cout << "Unknown key received" << std::endl;
		// 	return;
		// }
	}

	struct : Test_Function{
		double operator()(double x) const{
			return std::exp(x);
		}

		double deriv(double x) const{
			return std::exp(x);
		}

		double second_deriv(double x) const{
			return std::exp(x);
		}
	}	exp_test_func;

	struct : Test_Function{
		double c = 1.0;

		double operator()(double x) const{
			return 1.0/(1.0 + c*x*x);
		}

		double deriv(double x) const{
			return -2.0*c*x/std::pow(1.0 + c*x*x, 2.0);
		}

		double second_deriv(double x) const{
			return 2.0*c*(3.0*c*x*x - 1.0)/std::pow(1.0 + c*x*x, 3.0);
		}
	}	frac_test_func;

	// Staggered_Differentiation_Matrix stag_matr(n, -1.0, 1.0);

	// std::cout << std::setprecision(16);
	// for(int i = 0; i < n; ++i){
	// 	for(int j = 0; j < n; ++j)
	// 		std::cout << stag_matr.d1matr[i][j] << ' ';
	// 	std::cout << std::endl;
	// }

	Staggered_Differentiation_Matrix_Test stag_matr_test(n, -1.0, 1.0);
	std::cout << "Exponential function test" << std::endl;
	stag_matr_test.test_function(&exp_test_func);
	std::cout << "\nFractional function test" << std::endl;
	stag_matr_test.test_function(&frac_test_func);

	std::cout << '\n';
	for(int i = 0; i < 70; ++i)
		std::cout << '-';
	std::cout << '\n' << std::endl;
}


// Routines to test Spectrum Solver class
template<class FuncType>
Eigen::VectorXd fill_function(int nfunc, int n, double xst, double xend, const FuncType &func){
	double x;
	Eigen::VectorXd vec(nfunc*n);

	int indx = 0;
	for(int i = 0; i < n; ++i){
		x = xst + 0.5 * (xend - xst) * (std::cos(M_PI * (n - i - 0.5)/(1.0 * n)) + 1.0);
		for(int fn = 0; fn < nfunc; ++fn)
			vec(indx++) = func(x, fn);
	}

	return vec;
}

// Finding the spectrum of f'' + \omega^2 f = 0

class Simple_Oscillator{
public:
	const int nfunc = 1;
	double xst;
	double xend;

	Eigen::MatrixXd d2_block;
	Eigen::MatrixXd d1_block;
	Eigen::MatrixXd d0_block;
	Eigen::MatrixXd freq_block;

	Eigen::MatrixXd prefactor;
	Eigen::MatrixXd dprefactor;
	Eigen::MatrixXd ddprefactor;

	Simple_Oscillator(): xst(0.0), xend(M_PI), d2_block(nfunc, nfunc), d1_block(nfunc, nfunc),
							d0_block(nfunc, nfunc), freq_block(nfunc, nfunc),
							prefactor(nfunc, nfunc), dprefactor(nfunc, nfunc), ddprefactor(nfunc, nfunc) {};

	void apply_prefactor(double x){
		//Evaluating prefactor
		prefactor << std::sin(x);
		dprefactor << std::cos(x);
		ddprefactor << -std::sin(x);

		//Applying prefactor
		d0_block = d0_block*prefactor + d1_block*dprefactor + d2_block*ddprefactor;
		d1_block = d1_block*prefactor + 2.0*d2_block*dprefactor;
		d2_block = d2_block*prefactor;
		freq_block = freq_block*prefactor;
	}

	void evaluate_blocks(int i, double x){
		d2_block << 1.0;
		d1_block << 0.0;
		d0_block << 0.0;
		freq_block << 1.0;

		apply_prefactor(x);
	}
};

// Finding the spectrum of f' - g = 0
//						   g' + \omega^2 f = 0

class First_Order_Oscillator{
public:
	const int nfunc = 2;
	double xst;
	double xend;

	Eigen::MatrixXd d2_block;
	Eigen::MatrixXd d1_block;
	Eigen::MatrixXd d0_block;
	Eigen::MatrixXd freq_block;

	Eigen::MatrixXd prefactor;
	Eigen::MatrixXd dprefactor;
	Eigen::MatrixXd ddprefactor;

	First_Order_Oscillator(): xst(0.0), xend(M_PI), d2_block(nfunc, nfunc), d1_block(nfunc, nfunc),
							d0_block(nfunc, nfunc), freq_block(nfunc, nfunc),
							prefactor(nfunc, nfunc), dprefactor(nfunc, nfunc), ddprefactor(nfunc, nfunc) {};

	void apply_prefactor(double x){
		//Evaluating prefactor
		prefactor << std::sin(x), 0, 0, 1;
		dprefactor << std::cos(x), 0, 0, 0;
		ddprefactor << -std::sin(x), 0, 0, 0;

		//Applying prefactor
		d0_block = d0_block*prefactor + d1_block*dprefactor + d2_block*ddprefactor;
		d1_block = d1_block*prefactor + 2.0*d2_block*dprefactor;
		d2_block = d2_block*prefactor;
		freq_block = freq_block*prefactor;
	}

	void evaluate_blocks(int i, double x){
		d2_block << 0, 0,
					0, 0;
		d1_block << 1, 0,
					0, 1;
		d0_block << 0, -1,
					0, 0;
		freq_block << 0, 0,
					  1, 0;

		apply_prefactor(x);
	}
};

void spectrum_solver_test(int argc, char const *argv[]){
	std::cout << "Running generalized eigenvalue solver tests\n" << std::endl;
	std::cout << "Finding mode frequencies for f'' + omega^2 f = 0\n" << std::endl;

	int n = 20;
	
	for(int i = 1; i < argc; i++){
		if(strcmp(argv[i], "-n") == 0){
			n = atoi(argv[i+1]);
			i++;
		}
		// else{
		// 	std::cout << "Unknown key received" << std::endl;
		// 	return;
		// }
	}

	Simple_Oscillator osc;
	SpectrumSolver<Simple_Oscillator> solver(&osc, n);

	solver.compute_spectrum();
	solver.summary();

	// Eigen::VectorXd func_unit = fill_function(1, n, 0, M_PI, [](double x, int fn){return 1.0;});
	// Eigen::VectorXd func_sin = fill_function(1, n, 0, M_PI, [](double x, int fn){return std::sin(x);});
	// Eigen::VectorXd func_cos = fill_function(1, n, 0, M_PI, [](double x, int fn){return std::cos(x);});

	// Eigen::VectorXd res = solver.matr_a * func_unit + func_sin;
	// std::cout << res.cwiseAbs().maxCoeff() << std::endl;

	// Eigen::VectorXd res = solver.matr_a * solver.real_eigenvectors[0] - 
	// 						solver.real_freq[0]*solver.matr_b*solver.real_eigenvectors[0];
	// std::cout << res.cwiseAbs().maxCoeff() << std::endl;

	// solver.save_real_eigenvalues_csv("eigenvectors.csv");

	for(int i = 0; i < 70; ++i)
		std::cout << '-';
	std::cout << '\n' << std::endl;

	std::cout << "Running generalized eigenvalue solver tests\n" << std::endl;
	std::cout << "Finding mode frequencies for:\n";
	std::cout << "f' - g = 0\n";
	std::cout << "g' + omega^2 f = 0\n" << std::endl;

	First_Order_Oscillator osc_first;
	SpectrumSolver<First_Order_Oscillator> solver_first(&osc_first, n);

	solver_first.compute_spectrum(true);
	solver_first.summary();

	// Eigen::VectorXd func_sin = fill_function(2, n, 0, M_PI, [](double x, int fn){
	// 															if(fn == 0)
	// 																return 1.0;
	// 															else
	// 																return std::cos(x);
	// 															});

	// Eigen::VectorXd func_res = fill_function(2, n, 0, M_PI, [](double x, int fn){
	// 															if(fn == 0)
	// 																return 0.0;
	// 															else
	// 																return -std::sin(x);
	// 															});

	// Eigen::VectorXd res = solver_first.matr_a * func_sin - func_res;
	// for(int i = 0; i < n; ++i){
	// 	std::cout << res(2*i) << ", " << res(2*i+1) << std::endl;
	// }
	// std::cout << res.cwiseAbs().maxCoeff() << std::endl;


	int freq_num = 1;
	std::cout << "Verifying eigenvalue with the frequency: " << solver_first.real_freq[freq_num] << '\n' << std::endl;

	std::cout << "Corresponding eigenvector is:" << std::endl;
	for(int i = 0; i < n; ++i){
		std::cout << solver_first.real_eigenvectors[freq_num](2*i) << ", " << solver_first.real_eigenvectors[freq_num](2*i+1) << std::endl;
	}

	Eigen::VectorXd res = solver_first.matr_a * solver_first.real_eigenvectors[freq_num] - 
							solver_first.real_freq[freq_num]*solver_first.matr_b * solver_first.real_eigenvectors[freq_num];
	
	std::cout << "\nEigenvector error is: " << res.cwiseAbs().maxCoeff() << std::endl;

	// solver_first.save_real_eigenvalues_csv("eigenvectors.csv");

	std::cout << '\n';
	for(int i = 0; i < 70; ++i)
		std::cout << '-';
	std::cout << '\n' << std::endl;
}