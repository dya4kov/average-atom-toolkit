#pragma once
#include <cmath>
#include <functional>
#include <numeric-tools/ODE/types.h>
#include <numeric-tools/specfunc/fermi-dirac/complete.h>

namespace AATools {
namespace TF {
namespace ODE {

using ::numtools::ODE::Array;
using ::numtools::ODE::Dimension;

using ::numtools::specfunc::FermiDirac;
using ::numtools::specfunc::FD::Half;

using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;

class RotatePoints {
protected:
	struct RHSRP2 {
		static const Dimension dim = 2;

		double     xUp()   { return _xUp;   }
		double     xDown() { return _xDown; }
		Array<dim> yUp()   { return _yUp;   }
		Array<dim> yDown() { return _yDown; }

		RHSRP2() : V(1.0), mu(1.0), T(1.0), energyArg(1.0), lambdaArg(1.0) {
                _yOld.fill(-1.0); p1 = -1.0;
                _xOld =     1.0;  p2 = -1.0;
                rhs = std::bind(&RHSRP2::rhsT, this, _1, _2, _3);
                eval_a();
            }

        void reset  () { _xUp = 0.0; _xDown = 0.0; }
        void set_mu (const double& _mu) { mu = _mu; }
        void set_l  (const double& _l)  { lambdaArg = _l; }
        void set_e  (const double& _e)  { energyArg = _e; }
		void set_V  (const double& _V)  { V = _V; eval_a(); }
        void set_T  (const double& _T)  { T = _T;
            rhs = T > 1e-10 
                    ? std::bind(&RHSRP2::rhsT,  this, _1, _2, _3) 
                    : std::bind(&RHSRP2::rhsT0, this, _1, _2, _3);
        }

		void operator() (const double& x, Array<dim> &y, Array<dim> &dydx) { rhs(x, y, dydx); }
	private:
		void eval_a() { 
            a = std::pow(2.0, 7.0/6.0)
              * std::pow(3.0, 2.0/3.0)
              * std::pow(M_PI, -5.0/3.0)
              * std::pow(V, 2.0/3.0);
        }
        // actual rhs 
        std::function<void(const double&, Array<dim>&, Array<dim>&)> rhs;
        // rhs at T > 0
        void rhsT(const double& x, Array<dim>& y, Array<dim>& dydx) {
            double phi = y[0] + x*x*mu;
            dydx[0] = 2.0*x*y[1];
            if (x > 0) {
                dydx[1] = 2.0*a*x*x*x*std::sqrt(T)*T*FDhalf(phi/(T*x*x));
            }
            else {
                dydx[1] = 4.0/3.0*a*std::sqrt(phi)*phi;
            }
            p2 = energyArg*x*x + y[0] - lambdaArg / (x*x);
            
            if (p1 < 0 && p2 >= 0 && x < _xOld) {
                _yUp = _yOld; _yDown = y;
                _xUp = _xOld; _xDown = x;
            }

            p1 = p2;
            _yOld = y;
            _xOld = x;
        }
        // rhs at T = 0
        void rhsT0(const double& x, Array<dim>& y, Array<dim>& dydx) {
            double phi = y[0] + x*x*mu;
            dydx[0] = 2.0*x*y[1];
            dydx[1] = 4.0/3.0*a*std::sqrt(phi)*phi;

            p2 = energyArg*x*x + y[0] - lambdaArg / (x*x);
            
            if (p1 < 0 && p2 >= 0 && x < _xOld) {
                _yUp = _yOld; _yDown = y;
                _xUp = _xOld; _xDown = x;
            }

            p1 = p2;
            _yOld = y;
            _xOld = x;
        }

		double _xOld;
        double _xUp;
		double _xDown;

		Array<dim> _yOld;
		Array<dim> _yUp;
		Array<dim> _yDown;

		double a;
		double energyArg;
		double lambdaArg;
        double T;
        double V;
        double mu;
		double p1;
		double p2;
		
		FermiDirac<Half> FDhalf;
	};

	struct RHSRP1 {

		static const Dimension dim = 2;

        double     xUp()   { return _xUp;   }
		double     xDown() { return _xDown; }
		Array<dim> yUp()   { return _yUp;   }
		Array<dim> yDown() { return _yDown; }
		
		RHSRP1() : V(1.0), mu(1.0), T(1.0), energyArg(1.0), lambdaArg(1.0) {
                _yOld.fill(-1.0); 
                p1   = -1.0;
				p2   = -1.0;
                rhs = std::bind(&RHSRP1::rhsT, this, _1, _2, _3);
                eval_a();
            }

        void reset  () { _xUp = 0.0; _xDown = 0.0; }
        void set_mu (const double& _mu) { mu = _mu; }
        void set_l  (const double& _l)  { lambdaArg = _l; }
        void set_e  (const double& _e)  { energyArg = _e; }
        void set_V  (const double& _V)  { V = _V; eval_a(); }
        void set_T  (const double& _T)  { T = _T; 
            rhs = T > 1e-10 
                    ? std::bind(&RHSRP1::rhsT,  this, _1, _2, _3) 
                    : std::bind(&RHSRP1::rhsT0, this, _1, _2, _3);
        }

		void operator() (const double& x, Array<dim> &y, Array<dim> &dydx) { rhs(x, y, dydx); }
	private:
        void eval_a() { 
            a = std::pow(2.0, 7.0/6.0)
              * std::pow(3.0, 2.0/3.0)
              * std::pow(M_PI, -5.0/3.0)
              * std::pow(V, 2.0/3.0);
        }
        // actual rhs 
        std::function<void(const double&, Array<dim>&, Array<dim>&)> rhs;
        // rhs at T > 0
        void rhsT(const double& x, Array<dim>& y, Array<dim>& dydx) {
            double phi = y[0] + x*x*mu;
            dydx[0] = 2.0*x*y[1];
            if (x > 0) {
                dydx[1] = 2.0*a*x*x*x*std::sqrt(T)*T*FDhalf(phi/(T*x*x));
            }
            else {
                dydx[1] = 4.0/3.0*a*std::sqrt(phi)*phi;
            }
            p2 = energyArg*x*x + y[0] - lambdaArg / (x*x);

            if (p1 >= 0 && p2 < 0 && x < _xOld) {
                _yUp = _yOld; _yDown = y;
                _xUp = _xOld; _xDown = x;
            }

            p1 = p2;
            _yOld = y;
            _xOld = x;
        }
        // rhs at T = 0
        void rhsT0(const double& x, Array<dim>& y, Array<dim>& dydx) {
            double phi = y[0] + x*x*mu;
            dydx[0] = 2.0*x*y[1];
            dydx[1] = 4.0/3.0*a*std::sqrt(phi)*phi;
            p2 = energyArg*x*x + y[0] - lambdaArg / (x*x);

            if (p1 >= 0 && p2 < 0 && x < _xOld) {
                _yUp = _yOld; _yDown = y;
                _xUp = _xOld; _xDown = x;
            }

            p1 = p2;
            _yOld = y;
            _xOld = x;
        }

		double _xOld;
        double _xUp;
		double _xDown;

		Array<dim> _yOld;
		Array<dim> _yUp;
		Array<dim> _yDown;

		double a;
		double energyArg;
		double lambdaArg;
        double T;
        double V;
        double mu;
		double p1;
		double p2;

		FermiDirac<Half> FDhalf;
	};
};

}
}
}