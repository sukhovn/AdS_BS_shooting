#include <cmath>
#include <fftw3.h>

#ifndef STAG_DIFF_MATRIX_H
#define STAG_DIFF_MATRIX_H

// Quick and dirty way of computing Chebyshev differentiation matrix on a staggered grid

class Staggered_Differentiation_Matrix{
	unsigned int n;
	double xst;
	double xend;

	double *dmatr;
	double *ddmatr;

	void fill_x(){
		for(int i = 0; i < n; ++i)
			x[i] = xst + 0.5 * (xend - xst) * (std::cos(M_PI * (n - i - 0.5)/(1.0 * n)) + 1.0);
	}

	void cheb_diff(double *specti, double *specto){
		double tmp;
		double sp1 = specti[n-2];
		specto[n-2] = - 4.0*(n-1)*specti[n-1]/(xend-xst);
		specto[n-1] = 0.0;
		for(int i = n-3; i >= 0; i--){
			tmp = specti[i];
			specto[i] = specto[i+2] - 4.0*(i+1)*sp1/(xend-xst);
			sp1 = tmp;
		}
		return;
	}

	void fill_matrix(){
		double *point_spect = new double[n];
		double *column_tmp = new double[n];
		double tmp_pos;

		fftw_plan planCOS = fftw_plan_r2r_1d(n, point_spect, column_tmp, FFTW_REDFT01, FFTW_ESTIMATE);
		
		for(int i = 0; i < n; ++i){
			tmp_pos = M_PI*(i + 0.5)/(1.0*n);
			for(int j = 0; j < n; ++j)
				point_spect[j] = std::cos(tmp_pos*j)/(1.0*n);

			cheb_diff(point_spect, point_spect);
			fftw_execute_r2r(planCOS, point_spect, column_tmp);

			for(int j = 0; j < n; ++j)
				d1matr[j][i] = column_tmp[j];

			cheb_diff(point_spect, point_spect);
			fftw_execute_r2r(planCOS, point_spect, column_tmp);

			for(int j = 0; j < n; ++j)
				d2matr[j][i] = column_tmp[j];
		}

		fftw_destroy_plan(planCOS);

		delete [] point_spect;
		delete [] column_tmp;
	}

public:
	double *x;
	double **d1matr;
	double **d2matr;

	Staggered_Differentiation_Matrix(unsigned int n_, double xst_, double xend_): n(n_), xst(xst_), xend(xend_){
		x = new double[n];
		dmatr = new double[n*n];
		ddmatr = new double[n*n];
		d1matr = new double*[n];
		d2matr = new double*[n];

		double *tmpd = dmatr, *tmpdd = ddmatr;
		for(int i = 0; i < n; ++i){
			d1matr[i] = tmpd;
			d2matr[i] = tmpdd;
			tmpd += n;
			tmpdd += n;
		}

		fill_x();
		fill_matrix();
	}

	~Staggered_Differentiation_Matrix(){
		delete [] x;
		delete [] d1matr;
		delete [] d2matr;
		delete [] dmatr;
		delete [] ddmatr;
	}

	friend class Staggered_Differentiation_Matrix_Test;
};

// Testing routines

struct Test_Function{
	virtual double operator()(double x) const = 0;
	virtual double deriv(double x) const = 0;
	virtual double second_deriv(double x) const = 0;
};

class Staggered_Differentiation_Matrix_Test: public Staggered_Differentiation_Matrix{
public:
	Staggered_Differentiation_Matrix_Test(unsigned int n_, double xst_, double xend_):
		Staggered_Differentiation_Matrix(n_, xst_, xend_) {};

	void test_function(const Test_Function *test_func){
		double *func_vals = new double[n];
		double *deriv_vals = new double[n];
		double *second_deriv_vals = new double[n];

		for(int i = 0; i < n; ++i)
			func_vals[i] = (*test_func)(x[i]);

		for(int i = 0; i < n; ++i){
			deriv_vals[i] = 0.0;
			for(int j = 0; j < n; ++j)
				deriv_vals[i] += d1matr[i][j]*func_vals[j];

			second_deriv_vals[i] = 0.0;
			for(int j = 0; j < n; ++j)
				second_deriv_vals[i] += d2matr[i][j]*func_vals[j];
		}

		double tmp;
		double d1_max_diff = 0.0;
		double d2_max_diff = 0.0;

		for(int i = 0; i < n; ++i){
			tmp = std::fabs(deriv_vals[i] - test_func->deriv(x[i]));
			// std::cout << deriv_vals[i] - test_func->deriv(x[i]) << std::endl;
			if(tmp > d1_max_diff)
				d1_max_diff = tmp;
		
			tmp = std::fabs(second_deriv_vals[i] - test_func->second_deriv(x[i]));
			if(tmp > d2_max_diff)
				d2_max_diff = tmp;
		}

		std::ios oldState(nullptr);
		oldState.copyfmt(std::cout);

		std::cout << std::scientific;
		std::cout << "Testing Chebyshev differentiantion matrix for staggered grid calculator" << std::endl;
		std::cout << "Maximal difference in first derivative calculation: " << d1_max_diff << std::endl;
		std::cout << "Maximal difference in second derivative calculation: " << d2_max_diff << std::endl;

		std::cout.copyfmt(oldState);

		delete [] func_vals;
		delete [] deriv_vals;
		delete [] second_deriv_vals;
	}
};

#endif