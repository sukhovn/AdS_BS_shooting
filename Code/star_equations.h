#include <ode_solver.h>
#include <cmath>

#ifndef STAR_EQUATIONS_H
#define STAR_EQUATIONS_H

#ifndef SQR
#define SQR(a) a*a
#endif

// Solves qubic equation a*x^3 + b*x = c
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
	double el;

	Boson_Star_AdS() = delete;
	Boson_Star_AdS(int di, double li, double omi): RKDP_stepper(6),
		d(di), l(li), omega(omi), kappa(1.0), mass(1.0), el(0.1) {};

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
	// y[4] <- a(x)
	// y[5] <- a'(x)

	void derivatives(const double x, const double *y, double *dy) const{
		double cosx = std::cos(x), sinx = std::sin(x);
		double expdel = std::exp(2.0*y[1]);
		double grav_factor = (d - 2.0*SQR(cosx))/(sinx*cosx)*(y[0] - 1.0);
		double em_energy = sinx*cosx*SQR(y[5]*cosx)*expdel/(l*l*(d-1));
		double kinetic = 2.0*sinx*cosx*(expdel/y[0]*SQR((omega - el*y[4])*y[2]) + y[0]*y[3]*y[3])/(1.0*(d-1));
		double potential = 2.0*l*l*sinx/((d-1)*cosx)*vpot(y[2]*y[2]);
		double freq_part = (expdel/y[0]*SQR((omega - el*y[4])) - SQR(l/cosx)*dvpot(y[2]*y[2]))/y[0];
		
		dy[0] = -grav_factor - kappa*(kinetic + potential + em_energy);
		dy[1] = -kappa*kinetic/y[0];
		dy[2] = y[3];
		dy[3] = -(dy[0]/y[0] - dy[1] + (d-1)/(sinx*cosx))*y[3] - freq_part*y[2];
		dy[4] = y[5];
		dy[5] = - (dy[1] + (2*cosx*cosx + d - 3.0)/(sinx*cosx))*y[5] - 2.0*el*(omega - el*y[4])*SQR(l*y[2]/cosx)/y[0];
	}

	struct Origin_Asymptotic{
		double del0;
		double phi0;
		double a0;

		double A2;
		double del2;
		double phi2;
		double a2;

		double A4;
		double del4;
		double phi4;
		double a4;

		Origin_Asymptotic(const Boson_Star_AdS &bstar, double del0i, double phi0i, double a0i): del0(del0i), phi0(phi0i), a0(a0i){
			double expdel = std::exp(2.0*del0);
			double vpot0 = bstar.vpot(phi0*phi0);
			double dvpot0 = bstar.dvpot(phi0*phi0);
			double ddvpot0 = bstar.ddvpot(phi0*phi0);
			double omegahat = bstar.omega - bstar.el*a0;

			A2 = -2.0*bstar.kappa*(expdel*SQR(phi0*omegahat) + SQR(bstar.l)*vpot0)/(1.0*bstar.d*(bstar.d-1));
			del2 = -bstar.kappa*expdel*SQR(phi0*omegahat)/(1.0*(bstar.d-1));
			phi2 = -(expdel*SQR(omegahat) - SQR(bstar.l)*dvpot0)*phi0/(2.0*bstar.d);
			a2 = - bstar.el*omegahat*SQR(bstar.l*phi0)/(1.0*bstar.d);

			A4 = -2.0*(A2*(SQR(bstar.d) - 1.0) + 3.0*bstar.kappa*(2.0*expdel*SQR(a2/bstar.l) + 4.0*phi2*phi2 -
					expdel*phi0*omegahat*((A2*phi0 - 2.0*(phi2 + phi0*del2))*omegahat + 2.0*a2*bstar.el*phi0) - 2.0*expdel*SQR(omegahat*phi0)/3.0) +
					SQR(bstar.l)*bstar.kappa*(vpot0 + 6.0*phi0*phi2*dvpot0))/(3.0*(bstar.d+2)*(bstar.d-1));
			del4 = - bstar.kappa*(4.0*phi2*phi2 - 2.0*expdel*phi0*omegahat*(a2*bstar.el*phi0 + omegahat*(A2*phi0 - phi2 - phi0*del2)) -
					2.0*expdel*SQR(omegahat*phi0)/3.0)/(2.0*(bstar.d - 1));
			phi4 = -(4.0*phi2*(3.0*(A2 - del2) + bstar.d - 1)/3.0 - expdel*omegahat*phi0*(2.0*bstar.el*a2 + omegahat*(A2 - 2*del2)) +
					(phi2 - A2*phi0)*(expdel*SQR(omegahat) - SQR(bstar.l)*dvpot0) - SQR(bstar.l)*phi0*(dvpot0 + 2.0*phi0*phi2*ddvpot0))/(4.0*(bstar.d+2));
			a4 = -(2.0*a2*(2.0*(bstar.d - 4) + 6.0*del2 - 3.0*SQR(bstar.el*bstar.l*phi0))/3.0 -
					2.0*bstar.el*SQR(bstar.l)*phi0*omegahat*((A2 - 1.0)*phi0 - 2.0*phi2))/(4.0*(bstar.d+2));
		};

		void evaluate(double xa, double *ya){
			ya[0] = 1.0 + A2*xa*xa + A4*xa*xa*xa*xa;
			ya[1] = del0 + del2*xa*xa + del4*xa*xa*xa*xa;
			ya[2] = phi0 + phi2*xa*xa + phi4*xa*xa*xa*xa;
			ya[3] = 2.0*phi2*xa + 4.0*phi4*xa*xa*xa;
			ya[4] = a0 + a2*xa*xa + a4*xa*xa*xa*xa;
			ya[5] = 2.0*a2*xa + 4.0*a4*xa*xa*xa;
		}
	};

	void set_initial_conditions(const double xi, double *yi){
		x = xi;

		Origin_Asymptotic asym(*this, yi[0], yi[1], yi[2]);

		asym.evaluate(xi, y);
		derivatives(xi, y, dy);

		// for(int i = 0; i < 4; ++i) yi[i] = y[i];
		// std::cout << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << std::endl;
		// throw("Stop");
	}

	struct Boundary_Asymptotic{
		double da;
		double np;

		double rho;

		double phi0;
		double A0;
		double a0;
		double a1;

		double b1;
		double b2;
		double b3;
		double b4;
		double b5;

		double del0;

		// Variables
		// y[0] <- A(x)
		// y[1] <- delta(x)
		// y[2] <- phi(x)
		// y[3] <- phi'(x)
		// y[4] <- a(x)
		// y[5] <- a'(x)

		Boundary_Asymptotic(const Boson_Star_AdS &bstar){
			da = bstar.d;
			rho = M_PI_2 - bstar.x;
			double tmp = bstar.dvpot(0.0), tmpem;
			np = 0.5*(da + std::sqrt(SQR(da) + 4.0*tmp*SQR(bstar.l)));

			if(bstar.d < 3 && bstar.el != 0.0){
				throw("2D charged boson star not permitted");
			}

			// Determining electromagnetic potential asymptotic coefficients
			if(bstar.el == 0){
				a1 = 0.0;
				a0 = 0.0;
				b4 = 0.0;
			}
			else{
				a1 = - rho*bstar.y[5]/(da - 2.0)/(1.0 + da/3.0*rho*rho);
				a0 = bstar.y[4] - a1*(1.0 + (da - 2.0)/3.0*rho*rho);
				b4 = bstar.kappa*(da - 2.0)/(SQR(bstar.l)*(da - 1.0));
			}

			// Determining scalar field asymptotic coefficients
			b1 = np*(np + da - 2.0)/(6.0*(2.0*(np + 1.0) - da));
			b2 = SQR((bstar.omega - bstar.el*a0))/(4.0*(np + 1.0) - 2.0*da);
			b3 = bstar.kappa*np/(da - 1.0);
			b5 = 2.0*bstar.kappa*np*np/(2.0*np - da)/(da - 1.0);

			tmp = std::exp(2.0*bstar.y[1]);
			tmpem = b4*(1.0 + (2.0*da - 5.0)*rho*rho/3.0)*SQR(a1*rho);
			phi0 = solve_qubic(2.0*b2*b3*tmp*rho*rho + 0.5*b5/(1.0 + (da - 3.0)*rho*rho/3.0),
								1.0 + (b1 - b2*tmp)*rho*rho - 0.5*(bstar.y[0] - 1.0 - tmp*tmpem)/
								(1.0 + (da - 3.0)*rho*rho/3.0), bstar.y[2]);

			// Determining gravitational asymptotic coefficients
			del0 = bstar.y[1] - b3*phi0*phi0;
			b2 *= std::exp(2.0*del0);
			b4 *= std::exp(2.0*del0);

			A0 = (bstar.y[0] - 1.0 - b5*phi0*phi0 - std::exp(2.0*del0)*tmpem)/
					(1.0 + (da - 3.0)*rho*rho/3.0);
		}

		double evaluate_dphi(){
			return -(np/rho + (np + 2.0)*(b1 - b2)*rho)*phi0;
		}

		void evaluate(double xa, double *ya){
			double rhoa = M_PI_2 - xa;
			double rho_frac = rhoa/rho;
			ya[0] = 1.0 + A0*(1.0 + (da - 3.0)/3.0*rhoa*rhoa)*std::pow(rho_frac, da) + b5*std::pow(rho_frac, 2.0*np)*phi0*phi0 +
					b4*(1.0 + (2.0*da - 5.0)/3.0*rhoa*rhoa)*SQR(a1*rho)*std::pow(rho_frac, 2.0*da-2.0);
			ya[1] = del0 + b3*phi0*phi0*std::pow(rho_frac, 2.0*np);
			ya[2] = (1.0 + (b1 - b2)*rhoa*rhoa - 0.5*A0)*phi0*std::pow(rho_frac, np);
			ya[3] = -(np/rho + (np+2.0)*(b1 - b2)*rhoa*rhoa/rho - 0.5*(np + da)*A0/rho)*std::pow(rho_frac, np-1.0)*phi0;
			if(da > 2){
				ya[4] = a0 + a1*(1.0 + (da - 2.0)/3.0*rhoa*rhoa)*std::pow(rho_frac, da-2.0);
				double tmp = (da == 3.0 ? 1.0 : std::pow(rho_frac, da-3.0));
				ya[5] = -(da - 2.0)*a1*(tmp/rho)*(1.0 + da/3.0*rhoa*rhoa);
			}
			else{
				ya[4] = 0.0;
				ya[5] = 0.0;
			}
		}

		void evaluate_boundary(double *ya){
			ya[0] = 1.0;
			ya[1] = del0;
			ya[2] = 0.0;
			ya[3] = 0.0;
			if(da > 2){
				ya[4] = a0;
				if(da == 3)
					ya[5] = -(da - 2.0)*a1/rho;
				else
					ya[5] = 0.0;
			}
			else{
				ya[4] = 0.0;
				ya[5] = 0.0;
			}
		}
	};
};

#endif