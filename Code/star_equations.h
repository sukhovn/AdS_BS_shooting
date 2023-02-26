#include <ode_solver.h>
#include <cmath>

#ifndef STAR_EQUATIONS_H
#define STAR_EQUATIONS_H

#ifndef SQR
#define SQR(a) a*a
#endif

//Solves qubic equation a*x^3 + b*x = c
double solve_qubic(double a, double b, double c){
	if(std::fabs(a) > 1.0e-4){
		double det = std::sqrt(SQR(c/(2.0*a)) + std::pow(b/(3.0*a),3.0));
		return std::pow(det + c/(2.0*a), 1.0/3.0) - std::pow(det - c/(2.0*a), 1.0/3.0);
	}
	return c/b*(1.0 - a*c*c/std::pow(b, 3.0) + 3.0*a*a*std::pow(c, 4.0)/std::pow(b, 6.0));
}

class Boson_Star_AdS: public RKDP_stepper{
public:
	int d;
	double l;
	double omega;
	double kappa;
	double mass;

	Boson_Star_AdS() = delete;
	Boson_Star_AdS(int di, double li, double omi): RKDP_stepper(4), 
		d(di), l(li), omega(omi), kappa(1.0), mass(1.0) {};

	// Potential and potential derivatives
	double vpot(double phi2) const{
		return mass*mass*phi2;
	}
	
	double dvpot(double phi2) const{
		return mass*mass;
	}

	double ddvpot(double phi2) const{
		return 0.0;
	}

	// Variables
	// y[0] <- A(x)
	// y[1] <- delta(x)
	// y[2] <- phi(x)
	// y[3] <- phi'(x)

	void derivatives(const double x, const double *y, double *dy) const{
		double cosx = std::cos(x), sinx = std::sin(x);
		double expdel = std::exp(2.0*y[1]);
		double grav_factor = (d - 2.0*SQR(cosx))/(sinx*cosx)*(y[0] - 1.0);
		double kinetic = 2.0*sinx*cosx*(expdel/y[0]*SQR(omega*y[2]) + y[0]*y[3]*y[3])/(1.0*(d-1));
		double potential = 2.0*l*l*sinx/((d-1)*cosx)*vpot(y[2]*y[2]);
		double freq_part = (expdel/y[0]*SQR(omega) - SQR(l/cosx)*dvpot(y[2]*y[2]))/y[0];
		
		dy[0] = -grav_factor - kappa*(kinetic + potential);
		dy[1] = -kappa*kinetic/y[0];
		dy[2] = y[3];
		dy[3] = -(dy[0]/y[0] - dy[1] + (d-1)/(sinx*cosx))*y[3] - freq_part*y[2];
	}

	struct Origin_Asymptotic{
		double del0;
		double phi0;

		double A2;
		double del2;
		double phi2;

		double A4;
		double del4;
		double phi4;

		Origin_Asymptotic(const Boson_Star_AdS &bstar, double del0i, double phi0i): del0(del0i), phi0(phi0i){
			double expdel = std::exp(2.0*del0);
			double vpot0 = bstar.vpot(phi0*phi0);
			double dvpot0 = bstar.dvpot(phi0*phi0);
			double ddvpot0 = bstar.ddvpot(phi0*phi0);

			A2 = -2.0*bstar.kappa*(expdel*SQR(phi0*bstar.omega) + SQR(bstar.l)*vpot0)/(1.0*bstar.d*(bstar.d-1));
			del2 = -bstar.kappa*expdel*SQR(phi0*bstar.omega)/(1.0*(bstar.d-1));
			phi2 = -(expdel*SQR(bstar.omega) - SQR(bstar.l)*dvpot0)*phi0/(2.0*bstar.d);

			A4 = -2.0*(A2*(SQR(bstar.d) - 1.0) + bstar.kappa*(12.0*phi2*phi2 - (2.0 + 3.0*A2 - 6.0*del2)*expdel*SQR(phi0*bstar.omega)
							+ SQR(bstar.l)*vpot0 + 6.0*phi0*phi2*(expdel*SQR(bstar.omega) + SQR(bstar.l)*dvpot0)))/(3.0*(bstar.d+2)*(bstar.d-1));
			del4 = bstar.kappa*(expdel*phi0*(phi0*(1.0 + 3.0*(A2 - del2)) - 3.0*phi2)*SQR(bstar.omega) - 6.0*SQR(phi2))/(3.0*(bstar.d-1));
			phi4 = -(4.0*phi2*(3.0*(A2-del2) + bstar.d-1)/3.0 - expdel*phi0*(A2 - 2.0*del2)*SQR(bstar.omega) + 
						(phi2 - A2*phi0)*(expdel*SQR(bstar.omega) - SQR(bstar.l)*dvpot0) - SQR(bstar.l)*phi0*(dvpot0 + 2.0*phi0*phi2*ddvpot0))/
						(4.0*(bstar.d+2));
		};

		void evaluate(double xa, double *ya){
			ya[0] = 1.0 + A2*xa*xa + A4*xa*xa*xa*xa;
			ya[1] = del0 + del2*xa*xa + del4*xa*xa*xa*xa;
			ya[2] = phi0 + phi2*xa*xa + phi4*xa*xa*xa*xa;
			ya[3] = 2.0*phi2*xa + 4.0*phi4*xa*xa*xa;
		}
	};

	void set_initial_conditions(const double xi, double *yi){
		x = xi;

		Origin_Asymptotic asym(*this, yi[0], yi[1]);

		asym.evaluate(xi, y);
		derivatives(xi, y, dy);

		// for(int i = 0; i < 4; ++i) yi[i] = y[i];
		// std::cout << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << std::endl;
		// throw("Stop");
	}

	struct Boundary_Asymptotic{
		int da;

		double rho;
		double np;

		double phi0;
		double A0;

		double a1;
		double a2;
		double b1;
		double b2;

		double del0;

		Boundary_Asymptotic(const Boson_Star_AdS &bstar){
			da = bstar.d;
			rho = M_PI_2 - bstar.x;
			double m2 = bstar.dvpot(0.0); 
			np = 0.5*(da + std::sqrt(SQR(da) + 4.0*m2*SQR(bstar.l)));

			a1 = 1.0 + (da - 3.0)/3.0*rho*rho;
			a2 = 4.0*bstar.kappa*np*(np - 1.0)/((2.0*np - da)*(da - 1.0));
			b1 = (np*(np+da-2) - 3.0*SQR(bstar.omega))/(2.0*(np+1.0) - da)/6.0;
			if(da == 2)
				b2 = 3.0*np*np/(2.0*(np+1.0) - da)/6.0;
			else 
				b2 = 0.0;

			if(da == 2)
				phi0 = solve_qubic(a2*b2/a1, 1.0 + b1*rho*rho - b2/a1*(bstar.y[0] - 1.0), bstar.y[2]);
			else
				phi0 = bstar.y[2]/(1.0 + b1*rho*rho);

			A0 = (bstar.y[0] - 1.0 - a2*phi0*phi0)/a1;
			del0 = bstar.kappa*np*phi0*phi0/(1.0*(da - 1));
		}

		double evaluate_dphi(){
			return -(np/rho + (np + 2.0)*(b1*rho - b2*A0/rho))*phi0;
		}

		void evaluate(double xa, double *ya){
			double rhoa = M_PI_2 - xa;
			double rho_frac = rhoa/rho;
			ya[0] = 1.0 + A0*(1.0 + (da - 3.0)/3.0*rhoa*rhoa)*std::pow(rho_frac, da) + a2*std::pow(rho_frac, 2.0*np)*phi0*phi0;
			ya[1] = del0*std::pow(rho_frac, 2.0*np);
			ya[2] = (1.0 + b1*rhoa*rhoa - b2*A0*SQR(rho_frac))*phi0*std::pow(rho_frac, np);
			ya[3] = -(np/rho + (np+2.0)*(b1*rhoa*rho_frac - b2*SQR(rho_frac)*A0/rho))*std::pow(rho_frac, np-1)*phi0;
		}
	};
};

#endif