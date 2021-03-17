#include <cmath>
#include <iostream>

#include <average-atom-toolkit/atom/thomas-fermi.h>
#include <average-atom-toolkit/atom/semiclassic.h>

namespace aatk {
namespace atom {

using TFAtom = ::aatk::atom::ThomasFermiAtom;

// Atom::Atom(ConstructorArgs args) : 
// 	acc(nullptr),
//     phiSpline(nullptr),
//     dphiSpline(nullptr),
//     densSpline(nullptr),
//     eLevelStart({-1e+4, -1e+3, -1e+2, -1e+1, -1.0, 0.0, 1.0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5})
// {
// 	acc = gsl_interp_accel_alloc();
// 	reset(args);
// }

// old style
SemiclassicAtom::SemiclassicAtom(
	double _V,
	double _T,
	double _Z,
	double _tolerance,
	int    _meshSize,
	int    _nmax
) : acc(nullptr),
    phiSpline(nullptr),
    dphiSpline(nullptr),
    densSpline(nullptr),
    eLevelStart({-1e+4, -1e+3, -1e+2, -1e+1, -1.0, 0.0, 1.0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5})
{
	acc = gsl_interp_accel_alloc();
	reset(_V, _T, _Z, _tolerance, _meshSize, _nmax);
}

// void Atom::reset(ConstructorArgs args) {
// 	reset(
// 		args.V, 
// 		args.T, 
// 		args.Z, 
// 		args.tolerance, 
// 		args.meshSize, 
// 		args.nmax
// 	);
// }

void SemiclassicAtom::reset(
	double _V,
	double _T,
	double _Z,
	double _tolerance,
	int    _meshSize,
	int    _nmax
) {
	V = _V;
	T = _T;
	Z = _Z;
	tolerance = _tolerance;
	meshSize = _meshSize;
	nmax = _nmax;
	r0 = std::pow(3.0*V / 4.0 / M_PI, 1.0 / 3.0);
    // Thomas-Fermi by default
	TFAtom tfAtom(_V, _T, _Z, _tolerance, _meshSize);
    M = tfAtom.chemicalPotential();
    // setup mesh
    mesh.resize(meshSize);
    std::vector<double> xmesh(meshSize);
    double umin = 0.0;
	double umax = 1.0;
	for (int i = 0; i < meshSize; ++i) {
		mesh[i] = umin + i*(umax - umin)/(meshSize - 1);
		xmesh[i] = mesh[i]*mesh[i];
	}
	xmesh[0] = 1e-10;
	dens.resize(xmesh.size());
	tfAtom.electronDensity(xmesh.data(), dens.data(), xmesh.size());
	for (int i = 0; i < dens.size(); ++i) {
		double x = mesh[i]*mesh[i];
		dens[i] *= 4.0*M_PI*r0*r0*x*x;
	}
	if (densSpline != nullptr) gsl_spline_free(densSpline);
    densSpline = gsl_spline_alloc(gsl_interp_cspline, meshSize);
    gsl_spline_init(densSpline, mesh.data(), dens.data(), meshSize);
	pot.resize(meshSize, 0.0);
	dpot.resize(meshSize, 0.0);
	evaluate_potential();
	// prepare arrays for energy levels
    eLevel.clear();
    eLevelReady.clear();
    eLevel.resize(nmax + 1); 
    eLevelReady.resize(nmax + 1);
    eLevel[0].resize(0);
    eLevelReady[0].resize(0);
    for (std::size_t n = 1; n <= nmax; ++n) {
        eLevel[n].resize(n);
        if (n <= nmax) eLevelReady[n].resize(n);
        else eLevelReady[n].resize(n, false);
    }
    nUpdate = 0;
}

SemiclassicAtom::~SemiclassicAtom() {
	gsl_spline_free(phiSpline);
	gsl_spline_free(dphiSpline);
	gsl_spline_free(densSpline);
	gsl_interp_accel_free(acc);
}

void SemiclassicAtom::U(const double *x, double *y, std::size_t n) {
	for (std::size_t i = 0; i < n; ++i) {
		y[i]  = gsl_spline_eval(phiSpline, std::sqrt(x[i]), acc)/x[i];
	}
	return;
}

void SemiclassicAtom::xU(const double *x, double *y, std::size_t n) {
	for (std::size_t i = 0; i < n; ++i) {
		y[i] = gsl_spline_eval(phiSpline, std::sqrt(x[i]), acc);
	}
	return;
}

void SemiclassicAtom::x2dU(const double *x, double *y, std::size_t n) {
	for (std::size_t i = 0; i < n; ++i) {
		double u = std::sqrt(x[i]);
		double y0 = gsl_spline_eval(phiSpline, u, acc);
		double y1 = gsl_spline_eval(dphiSpline, u, acc);
		y[i]  = 0.5*y1*u + y0;
	}
	return;
}

void SemiclassicAtom::electronDensity(const double* x, double* dens, std::size_t n, double eb) {
	for (std::size_t i = 0; i < n; ++i) {
		dens[i] = gsl_spline_eval(densSpline, std::sqrt(x[i]), acc);
	}
	return;
}

void SemiclassicAtom::update(double mixing) {
	// 1. evaluate energy levels
    // if(useContinuous){

    //     int n;
    //     double E_curr = 0;

    //     for ( n = 1; n <= Zcharge && E_curr <= 0; ++n) {
    //         if (n > nmax){
    //             nmax++;
    //             std::vector<double> l_vector(nmax);
    //             std::vector<bool> l_bool_vector(nmax);
    //             eLevel.push_back(l_vector);
    //             eLevelReady.push_back(l_bool_vector);
    //         }

    //         for (int l = 0; l < n; ++l) {
    //             evaluateEnergyLevel(n, l);
    //             E_curr = eLevel[n][l];
    //         }
    //     }
    //     nmax = --n;
    // }
    // else{
        for (int n = 1; n <= nmax; ++n) {
            for (int l = 0; l < n; ++l) {
                evaluate_energy_level(n, l);
            }
        }
    // }
	// 2. evaluate new chemical potential
	chemPotReady = false;
	evaluate_chemical_potential();
	// 3. evaluate new electron density
// #ifdef ENABLE_MULTITHREADING
// 	if (pool.size() > 0) {
// 		std::size_t num_tasks = nmax*(nmax + 1)/2;
// 		std::vector<std::future<std::vector<double>>> results(num_tasks);

// 		std::size_t i = 0;
//     	for (int n = 1; n <= nmax; ++n) {
//     		for (int l = 0; l < n; ++l) {
//     			double enl = energyLevel(n, l);
//     			double lambda = l + 0.5;
//     			results[i] = pool.enqueue(&Atom::waveFunctionVec, this, enl, lambda, std::cref(x));
//     			++i;	
//     		}
//     	}

//     	i = 0;
//     	for (int n = 1; n <= nmax; ++n) {
//     		for (int l = 0; l < n; ++l) {
//     			auto Rnl = results[i].get();
//     			auto Nnl = electronStatesDiscrete(n, l);
//     			for (std::size_t k = 0; k < x.size(); ++k) density[k] += Nnl*Rnl[k]*Rnl[k];
//     			++i;
//     		}
// 		}
// 	}
// 	else 
// #endif
	// {
		auto density = std::vector<double>(mesh.size(), 0.0);
		std::vector<double> Rnl(mesh.size());
		for (int n = 1; n <= nmax; ++n) {
    		for (int l = 0; l < n; ++l) {
    			auto enl = energyLevel(n, l);
                double lambda = l + 0.5;
    			// auto Rnl = waveFunction(enl, lambda, x);
    			evaluate_wave_function(mesh.data(), Rnl.data(), mesh.size(), enl, lambda);
    			auto Nnl = electronStatesDiscrete(n, l);
    			for (std::size_t k = 0; k < mesh.size(); ++k) density[k] += Nnl*Rnl[k]*Rnl[k];
    		}
    	}
        // if (useContinuous){
        //     for (std::size_t k = 0; k < x.size(); ++k) {
        //         density[k] += electronDensityContinuous(x[k])*(4.0 * M_PI * pow(x[k]*r0,2.0));
        //     }
        // }
	// }
    if (nUpdate > 0) {
        // 4. mixing with previous
        for (std::size_t k = 0; k < mesh.size(); ++k) 
        	density[k] = mixing*gsl_spline_eval(densSpline, mesh[k], acc) + (1.0 - mixing)*density[k];
    }
    // densityInterpolation = new Spline(u, density);
	gsl_spline_free(densSpline);
	densSpline = gsl_spline_alloc(gsl_interp_cspline, meshSize);
	dens = density;
    gsl_spline_init(densSpline, mesh.data(), dens.data(), meshSize);
	evaluate_potential();
    nUpdate++;
}

}
}