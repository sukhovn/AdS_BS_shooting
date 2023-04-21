#include <star_equations.h>
#include <eigen3/Eigen/Dense>

class Boson_Star_Profile{
	int npoints;
	double *xstorage;
	double *ystorage;

public:
	double xst;
	double xend;

	double phi0;
	double delta0;
	double a0;

	Boson_Star_AdS bstar;

	double *xarr;
	double *yarr;

	Boson_Star_Profile(int d, double l, double omega, double phi0_, double delta0_, double a0_, double xst_, double xend_):
						bstar(d, l, omega), xst(xst_), xend(xend_), phi0(phi0_), delta0(delta0_), a0(a0_),
						xstorage(nullptr), ystorage(nullptr), xarr(nullptr), yarr(nullptr), npoints(0) {};

	~Boson_Star_Profile(){
		if(xstorage != nullptr)
			delete [] xstorage;
		if(ystorage != nullptr)
			delete [] ystorage;
	}

	void gauge_out_omega(){
		a0 -= bstar.omega/bstar.el;
		bstar.omega = 0;
	}

	template<class Fill_Func>
	void fill_solution(int npoints_, const Fill_Func &func){
		npoints = npoints_;
		xstorage = new double[npoints+1];
		ystorage = new double[6*(npoints+1)];
		xarr = xstorage+1;
		yarr = ystorage+6;

		xstorage[0] = 0.0;
		for(int i = 0; i < npoints; ++i)
			xarr[i] = func(npoints, i);

		for(int i = 0; i < 6*(npoints+1); ++i)
			ystorage[i] = 0.0;

		int ist = 0, iend;
		double *xint = xarr;
		double *yint = yarr;
		double xst_tmp, xend_tmp;

		while(xarr[ist] < xst)
			++ist;
		--ist;

		xint += ist;
		yint += 6*ist;
		xst_tmp = xint[0];

		Boson_Star_AdS::Origin_Asymptotic asym_orig(bstar, delta0, phi0, a0);

		for(iend = 0; iend < npoints-1; ++iend){
			if(xarr[iend] > xend)
				break;
		}

		xend_tmp = xarr[iend];
		if(xend_tmp > xend){
			xarr[iend] = xend;
		}

		// Integrating
		yint[0] = delta0;
 		yint[1] = phi0;
		yint[2] = a0;
 		xint[0] = xst;

 		bstar.integrate_grid(iend-ist, xint, yint);

 		// Filling values before xst with asymptotics
 		xint[0] = xst_tmp;
 		yint = yarr;
 		for(int i = 0; i <= ist; ++i){
 			asym_orig.evaluate(xarr[i], yint);
			yint += 6;
 		}

 		if(xend_tmp <= xend)
 			return;

		// Filling values after xend with asymptotics
 		Boson_Star_AdS::Boundary_Asymptotic asym_bound(bstar);

 		xarr[iend] = xend_tmp;
		yint = yarr + 6*iend;
		// yint[0] = -10;
 		for(int i = iend; i < npoints; ++i){
 			if(M_PI_2 - xarr[i] < 1.0e-10)
 				asym_bound.evaluate_boundary(yint);
 			else
	 			asym_bound.evaluate(xarr[i], yint);
			yint += 6;
 		}
	}
};


class Boson_Star_Perturbations{
public:
	const int nfunc = 3;
	double xst;
	double xend;

	Boson_Star_Profile star_profile;

	Eigen::MatrixXd d2_block;
	Eigen::MatrixXd d1_block;
	Eigen::MatrixXd d0_block;
	Eigen::MatrixXd freq_block;

	Eigen::MatrixXd prefactor;
	Eigen::MatrixXd dprefactor;
	Eigen::MatrixXd ddprefactor;

	Boson_Star_Perturbations(int d, double l, double omega, double phi0_, double delta0_, double a0_, double xst_, double xend_): 
							xst(0.0), xend(M_PI_2), d2_block(nfunc, nfunc), d1_block(nfunc, nfunc),
							d0_block(nfunc, nfunc), freq_block(nfunc, nfunc),
							prefactor(nfunc, nfunc), dprefactor(nfunc, nfunc), ddprefactor(nfunc, nfunc),
							star_profile(d, l, omega, phi0_, delta0_, a0_, xst_, xend_) {};

	void update_lagrangian_parameters(double kappa, double mass, double el){
		star_profile.bstar.kappa = kappa;
		star_profile.bstar.mass = mass;
		star_profile.bstar.el = el;
	}

	void evaluate_unperturbed_fields(int npoints){
		star_profile.gauge_out_omega();
		star_profile.fill_solution(npoints, [](int n, int i){
			return M_PI_4 * (1.0 - std::cos(M_PI * (i + 0.5)/(1.0 * n)));
		});
	}

	// Fields go as {\phi_1, \phi_2, a_1}
	void apply_prefactor(double x){
		double sin_x = std::sin(x);
		double cos_x = std::cos(x);

		//Evaluating prefactor
		prefactor << std::pow(cos_x, 3.0),	0,	0,
					 0,	1,	0,
					 0,	0, std::pow(cos_x, 2.0);
		dprefactor << -3.0*cos_x*cos_x*sin_x,	0,	0,
					  0,	0,	0,
					  0,	0, -2.0*cos_x*sin_x;
		ddprefactor << -3.0*std::pow(cos_x, 3.0) + 6.0*cos_x*sin_x*sin_x, 0, 0,
						0,	0,	0,
						0,	0,	-2.0*(cos_x*cos_x - sin_x*sin_x);

		//Applying prefactor
		d0_block = d0_block*prefactor + d1_block*dprefactor + d2_block*ddprefactor;
		d1_block = d1_block*prefactor + 2.0*d2_block*dprefactor;
		d2_block = d2_block*prefactor;
		freq_block = freq_block*prefactor;
	}

	void evaluate_blocks(int indx, double x){
		if(std::fabs(x - star_profile.xarr[indx]) > 1.0e-12){
			std::cout << "Mismatch between unperturbed profile and collocation greed detected" << std::endl;
			throw("grid mismatch");
		}
		double d = star_profile.bstar.d;
		double &l = star_profile.bstar.l;
		double &kappa = star_profile.bstar.kappa;
		double &m = star_profile.bstar.mass;
		double &el = star_profile.bstar.el;
		double *y = star_profile.yarr + 6*indx;

		d2_block << 1.0, -2*el*y[4]*y[2], 0.0,
					0.0, 1.0, 0.0,
					0.0, (2*el*std::pow(l, 2))/std::pow(std::cos(x), 2)*std::exp(-2*y[1])*y[0]*std::pow(y[2], 2), 1.0;
		
		d1_block << (d - 1 - std::cos(2*x)*(1 - y[0]))/(y[0]*std::sin(x)*std::cos(x)) - (std::exp(2*y[1])*kappa*std::pow(std::cos(x), 3)*std::sin(x)*std::pow(y[5], 2))/((d - 1)*std::pow(l, 2)*y[0]) - (2*std::pow(l, 2)*std::pow(m, 2)*kappa*std::tan(x)*std::pow(y[2], 2))/((d - 1)*y[0]),
					(2*el*kappa*std::sin(2*x)*std::pow(y[2], 2)*y[5]*y[3])/(d - 1) - (2*el*y[2]*y[4])/y[0]*((d - 1 - std::cos(2*x)*(1 - y[0]))/(std::sin(x)*std::cos(x)) - (std::exp(2*y[1])*kappa*std::pow(std::cos(x), 3)*std::sin(x)*std::pow(y[5], 2))/((d - 1)*std::pow(l, 2)) - (4*std::pow(l, 2)*std::pow(m, 2)*kappa*std::tan(x)*std::pow(y[2], 2))/(d - 1)) - 4*el*y[4]*y[3] - (4*el*kappa*y[4]*std::pow(y[2], 2))/(std::pow(d - 1, 2)*std::pow(l, 2)*y[0])*((d - 1)*std::pow(l, 2)*(d - 1 - std::cos(2*x)) - std::exp(2*y[1])*kappa*std::pow(std::cos(x), 4)*std::pow(std::sin(x), 2)*std::pow(y[5], 2) - 2*std::pow(l, 4)*std::pow(m, 2)*kappa*std::pow(std::sin(x), 2)*std::pow(y[2], 2))*y[3],
					0.0,//
					0.0,
					(d - 1 - std::cos(2*x)*(1 - y[0]))/(y[0]*std::sin(x)*std::cos(x)) - (2*std::pow(el, 2)*std::exp(2*y[1])*kappa*std::pow(y[4], 2)*std::sin(2*x)*std::pow(y[2], 2))/((d - 1)*std::pow(y[0], 2)) - (std::pow(l, 2)*std::pow(y[2], 2))/(y[0]*std::pow(std::cos(x), 2))*((std::pow(m, 2)*kappa*std::sin(2*x))/(d - 1) + (2*std::pow(el, 2)*y[4])/y[5]) - (std::exp(2*y[1])*kappa*std::pow(std::cos(x), 3)*std::sin(x)*std::pow(y[5], 2))/(std::pow(l, 2)*y[0]*(d - 1)) + (2*y[3])/y[2],
					-(el*std::exp(2*y[1])*y[4])/(std::pow(y[0], 2)*y[5]),//
					-(2*kappa*std::sin(2*x)*y[5]*y[3])/(d - 1),
					2*el*y[2]*((std::pow(l, 2)*std::exp(-2*y[1])*(d - 1 - (1 - y[0])*std::cos(2*x)))/(std::sin(x)*std::pow(std::cos(x), 3))*y[2] - (kappa*std::cos(x)*std::sin(x)*y[2]*std::pow(y[5], 2))/(d - 1) + 2*std::pow(el, 2)*y[4]*std::pow(y[2], 3)*((2*std::pow(l, 2)*kappa*y[4]*std::tan(x))/((d - 1)*y[0]) - std::pow(std::pow(l, 2)/std::pow(std::cos(x), 2), 2)*std::exp(-2*y[1])/y[5] + (std::exp(2*y[1])*std::pow(kappa, 2)*std::pow(y[4], 2)*std::pow(std::sin(2*x), 2)*y[5])/(std::pow(d - 1, 2)*std::pow(y[0], 2))) + (2*std::pow(l, 2))/std::pow(std::cos(x), 2)*std::exp(-2*y[1])*y[0]*y[3] - (2*std::exp(-2*y[1])*std::pow(l, 4)*std::pow(m, 2)*kappa*std::sin(x)*std::pow(y[2], 3))/((d - 1)*std::pow(std::cos(x), 3))),
					(d - 2 + std::cos(2*x))/(std::sin(x)*std::cos(x)) - (2*std::pow(el, 2)*std::pow(l, 2)*y[4]*y[0]*std::pow(y[2], 2))/(y[5]*std::pow(y[0], 2)*std::pow(std::cos(x), 2)) + (kappa*std::sin(2*x)*(std::pow(el, 2)*std::exp(2*y[1])*std::pow(y[4], 2)*std::pow(y[2], 2) - std::pow(y[0], 2)*std::pow(y[3], 2)))/(std::pow(y[0], 2)*(d - 1));

		d0_block << -(3*std::pow(el, 2)*std::exp(2*y[1])*std::pow(y[4], 2))/std::pow(y[0], 2) - (std::pow(l, 2)*std::pow(m, 2))/(y[0]*std::pow(std::cos(x), 2)) + (4*kappa*y[3])/(y[0]*(d - 1))*(-2*std::pow(l, 2)*std::pow(m, 2)*std::tan(x)*y[2] + ((d - 1 - std::cos(2*x)) - (2*std::pow(l, 2)*std::pow(m, 2)*kappa*std::pow(std::sin(x), 2)*std::pow(y[2], 2))/(d - 1) - (std::exp(2*y[1])*kappa*std::pow(std::cos(x), 4)*std::pow(std::sin(x), 2)*std::pow(y[5], 2))/((d - 1)*std::pow(l, 2)))*y[3]),
					0.0,
					0.0,//
					(2*el*std::exp(2*y[1])*y[4]*(d - 1 + kappa*std::sin(2*x)*y[2]*y[3]))/((d - 1)*std::pow(y[0], 2)*y[2]),//
					0.0,
					(el*std::exp(2*y[1]))/std::pow(y[0], 2),//
					-(2*std::pow(el, 2)*std::exp(2*y[1])*kappa*std::pow(y[4], 2)*std::sin(2*x)*y[2]*y[5]*(d - 1 + 2*kappa*std::sin(2*x)*y[2]*y[3]))/(std::pow(d - 1, 2)*std::pow(y[0], 2)),
					0.0,
					-(2*std::pow(el, 2)*std::exp(2*y[1])*kappa*y[4]*std::sin(2*x)*std::pow(y[2], 2)*y[5])/((d - 1)*std::pow(y[0], 2));
					
		freq_block << std::exp(2*y[1])/std::pow(y[0], 2), 0.0, 0.0,
					  0.0, std::exp(2*y[1])/std::pow(y[0], 2), 0.0,
					  0.0, -(2*el*std::exp(2*y[1])*kappa*y[4]*std::sin(2*x)*std::pow(y[2], 2)*y[5])/((d - 1)*std::pow(y[0], 2)), 0.0;

		apply_prefactor(x);
	}
};