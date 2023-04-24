#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <stag_diff_matrix.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <fftw3.h>

#ifndef EIGEN_SOLVER_H
#define EIGEN_SOLVER_H

#define SMALL 1.0e-10

class EigenSolver{
public:
	const int size;
	const int npoints;
	const int nfunc;

	Eigen::MatrixXd matr_a;
	Eigen::MatrixXd matr_b;
	Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges;

public:
	int infinite_freq;
	std::vector<double> real_freq;
	std::vector<std::complex<double>> complex_freq;

	std::vector<Eigen::VectorXd> real_eigenvectors;
	std::vector<Eigen::VectorXcd> complex_eigenvectors;

	EigenSolver(unsigned int npoints_, unsigned int nfunc_): npoints(npoints_), nfunc(nfunc_), size(nfunc_*npoints_), infinite_freq(0){
		matr_a.resize(size, size);
		matr_b.resize(size, size);
		matr_a.setZero();
		matr_b.setZero();
	}

	void update_a(unsigned int bl_i, unsigned int bl_j, const Eigen::MatrixXd &block){
		if(bl_i >= npoints || bl_j >= npoints){
			std::cout << "Modifying block (" << bl_i << ", " << bl_j << ") out of range" << std::endl;
			throw("block out of range");
		}

		matr_a.block(nfunc*bl_i, nfunc*bl_j, nfunc, nfunc) = block;
	}

	void update_b(unsigned int bl_i, unsigned int bl_j, const Eigen::MatrixXd &block){
		if(bl_i >= npoints || bl_j >= npoints){
			std::cout << "Modifying block (" << bl_i << ", " << bl_j << ") out of range" << std::endl;
			throw("block out of range");
		}

		matr_b.block(nfunc*bl_i, nfunc*bl_j, nfunc, nfunc) = block;
	}

	void compute_eigenvalues(bool compute_eigenvectors = false){
		infinite_freq = 0;
		real_freq.resize(0);
		complex_freq.resize(0);

		ges.compute(matr_a, matr_b, compute_eigenvectors);
		Eigen::VectorXd betas = ges.betas();
		Eigen::VectorXcd eigval = ges.eigenvalues();

		std::vector<int> real_eigenvalues;
		std::vector<int> complex_eigenvalues;

		for(int i = 0; i < betas.size(); ++i){
			if(std::fabs(betas(i)) < SMALL){
				++infinite_freq;
				continue;
			}
			if(std::fabs(eigval(i).imag()) < SMALL){
				real_eigenvalues.push_back(i);
				continue;
			}
			complex_eigenvalues.push_back(i);
		}

		auto compare_real_freq = [&eigval](int i, int j){return eigval(i).real() < eigval(j).real();};
		std::sort(real_eigenvalues.begin(), real_eigenvalues.end(), compare_real_freq);
		std::sort(complex_eigenvalues.begin(), complex_eigenvalues.end(), compare_real_freq);

		for(int i : real_eigenvalues){
			real_freq.push_back(eigval(i).real());
			if(compute_eigenvectors)
				real_eigenvectors.push_back(ges.eigenvectors().col(i).real());
		}

		for(int i : complex_eigenvalues){
			complex_freq.push_back(eigval(i));
			if(compute_eigenvectors)
				complex_eigenvectors.push_back(ges.eigenvectors().col(i));
		}
	}

	void summary(){
		std::ios oldState(nullptr);
		oldState.copyfmt(std::cout);

		std::cout << "Number of infinite eigenvalues: " << infinite_freq << '\n' << std::endl;

		// std::cout << std::scientific;
		std::cout << "Number of real eigenvalues: " << real_freq.size() << std::endl;
		for(double el : real_freq){
			std::cout << el << std::endl;
		}
		std::cout << std::endl;

		std::cout << "Number of complex eigenvalues: " << complex_freq.size() << std::endl;
		for(std::complex<double> el : complex_freq){
			std::cout << '(' << el.real() << ", " << el.imag() << ')' << std::endl;
		}
		std::cout << std::endl;

		std::cout.copyfmt(oldState);
	}

	template <class>
	friend class SpectrumSolver;
};

template<class Equations_Class>
class SpectrumSolver: public EigenSolver{
	void fill_matrices(){
		matr_a.setZero();
		matr_b.setZero();

		for(int i = 0; i < npoints; ++i){
			eq_system->evaluate_blocks(i, diff_matrix.x[i]);
			update_a(i, i, eq_system->d2_block*diff_matrix.d2matr[i][i]+
						   eq_system->d1_block*diff_matrix.d1matr[i][i]+
						   eq_system->d0_block);

			update_b(i, i, -eq_system->freq_block);

			for(int j = 0; j < npoints; ++j){
				if(i == j) continue;
				update_a(i, j, eq_system->d2_block*diff_matrix.d2matr[i][j]+
							   eq_system->d1_block*diff_matrix.d1matr[i][j]);
			}
		}
	}

public:
	Staggered_Differentiation_Matrix diff_matrix;

	Equations_Class *eq_system;

	SpectrumSolver(Equations_Class *eq_sys, unsigned int npoints_):
		EigenSolver(npoints_, eq_sys->nfunc), diff_matrix(npoints_, eq_sys->xst, eq_sys->xend), eq_system(eq_sys) {};

	void compute_spectrum(bool compute_eigenvectors = false){
		fill_matrices();
		compute_eigenvalues(compute_eigenvectors);
	}

	void save_real_eigenvalues_csv(const std::string &file_name){
		std::ofstream out(file_name);

		out << 'x';

		for(int i = 0; i < real_eigenvectors.size(); ++i){
			for(int j = 0; j < nfunc; ++j){
				out << ",eig_" << i << "_f" << j;
			}
		}

		out << std::endl;

		for(int i = 0; i < npoints; ++i){
			out << diff_matrix.x[i];
			for(int j = 0; j < real_eigenvectors.size(); ++j){
				for(int fn = 0; fn < nfunc; ++fn)
					out << ',' << real_eigenvectors[j](nfunc*i+fn);
			}
			out << std::endl;
		}
	}

	//Checks decay of spectral coefficients
	void eigenvectors_precision(int nmodes, double *precision){
		double *spectra = new double[size];
		fftw_r2r_kind kind = FFTW_REDFT10;
		fftw_plan planCOS = fftw_plan_many_r2r(1, &npoints, nfunc, spectra, NULL, nfunc, 1, spectra, NULL, nfunc, 1, &kind, FFTW_ESTIMATE);

		int indx = 0;

		for(int mode = 0; mode < nmodes; ++mode){
			for(int i = 0; i < size; ++i){
				spectra[i] = real_eigenvectors[mode](i);
			}

			fftw_execute_r2r(planCOS, spectra, spectra);

			double tmp_front = 0.0, tmp_back = 0.0;
			for(int i = 0; i < 3*nfunc; ++i){
				tmp_front += spectra[i];
				tmp_back += spectra[size-1-i];
			}

			precision[indx++] = std::fabs(tmp_back/tmp_front);
		}

		fftw_destroy_plan(planCOS);
		delete [] spectra;
	}
};

#endif