#include <cmath>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <average-atom-toolkit/semiclassic/atom.h>
#include <average-atom-toolkit/thomas-fermi/atom/electron-density.h>
#include <average-atom-toolkit/thomas-fermi/eos/chemical-potential.h>

namespace aatk {
namespace semiclassic {

using ::aatk::TF::ChemicalPotential;
using ::aatk::TF::ElectronDensity;

static const int DTFsize = 601;

Atom::Atom(double _V, double _T, double _Z, int _nmax, bool _useContinuous,  double _tolerance
#ifdef ENABLE_MULTITHREADING
    ,ThreadPool& threads
#endif
	) :
	useContinuous(_useContinuous),
    tolerance(_tolerance),
    eLevelStart({-1e+3, -1e+2, -1e+1, -1.0, 0.0, 1.0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5})
#ifdef ENABLE_MULTITHREADING
   ,pool(threads)
#endif
{
	reset(_V, _T, _Z, _nmax);
}

Atom::~Atom() {
	delete densityInterpolation;
}

double Atom::V() {
	return volume;
}
double Atom::T() {
	return temperature;
}
double Atom::Z() {
	return Zcharge;
}
double Atom::M() {
	return chemPot;
}

double Atom::N_max(){
    return nmax;
}

std::vector<double> Atom::sorted_mesh(const double* mesh, std::size_t size) {
	std::vector<std::size_t> idx(size);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
       [&mesh](std::size_t i1, std::size_t i2) {
        return mesh[i1] < mesh[i2];
       }
    );
    std::vector<double> x(size);
    std::size_t k = 0;
    for (auto i : idx) {
		x[k] = mesh[i]; ++k;
    }
    return x;
}

void Atom::update(double mixing) {
	std::vector<double> x(DTFsize);
	double umin = 1e-3;
	double umax = 1.0;

	for (int i = 0; i < DTFsize; ++i) {
		double u = umin + i*(umax - umin)/(DTFsize - 1);
		x[i] = u*u;
	}
	update(x.data(), x.size(), mixing);
}

void Atom::update(const std::vector<double>& mesh, double mixing) {
    auto x = sorted_mesh(mesh.data(), mesh.size());
	auto d = electronDensity(mesh);
	auto u = x; for (auto& u_i : u) u_i = std::sqrt(u_i);
	densityInterpolation = new Spline(u, d);
}

void Atom::update(const double* mesh, std::size_t size, double mixing) {
    auto x = sorted_mesh(mesh, size);
	auto density = std::vector<double>(x.size(), 0.0);
	auto u = x; for (auto& u_i : u) u_i = std::sqrt(u_i);
	// 1. evaluate energy levels
#ifdef ENABLE_MULTITHREADING
	if (pool.size() > 0) {
		std::size_t num_tasks = nmax*(nmax + 1)/2;
		std::vector<std::future<int>> results(num_tasks);

		std::size_t i = 0;
    	for (int n = 1; n <= nmax; ++n) {
    		for (int l = 0; l < n; ++l) {
    			results[i] = pool.enqueue(&Atom::evaluateEnergyLevel, this, n, l);
    			++i;	
    		}
    	}

    	for (auto &result: results) result.get();
	}
	else 
#endif
	{
        if(useContinuous){

            int n;
            double E_curr = 0;

            for ( n = 1; n <= Zcharge && E_curr <= 0; ++n) {
                if (n > nmax){
                    nmax++;
                    std::vector<double> l_vector(nmax);
                    std::vector<bool> l_bool_vector(nmax);
                    eLevel.push_back(l_vector);
                    eLevelReady.push_back(l_bool_vector);
                }

                for (int l = 0; l < n; ++l) {
                    evaluateEnergyLevel(n, l);
                    E_curr = eLevel[n][l];
                }
            }
            nmax = --n;
        }
        else{
            for (int n = 1; n <= nmax; ++n) {
                for (int l = 0; l < n; ++l) {
                    evaluateEnergyLevel(n, l);
                }
            }
        }
	}
	// 2. evaluate new chemical potential
	// chemPotReady = false;
	evaluateChemicalPotential();
	// 3. evaluate new electron density
#ifdef ENABLE_MULTITHREADING
	if (pool.size() > 0) {
		std::size_t num_tasks = nmax*(nmax + 1)/2;
		std::vector<std::future<std::vector<double>>> results(num_tasks);

		std::size_t i = 0;
    	for (int n = 1; n <= nmax; ++n) {
    		for (int l = 0; l < n; ++l) {
    			double enl = energyLevel(n, l);
    			double lambda = l + 0.5;
    			results[i] = pool.enqueue(&Atom::waveFunctionVec, this, enl, lambda, std::cref(x));
    			++i;	
    		}
    	}

    	i = 0;
    	for (int n = 1; n <= nmax; ++n) {
    		for (int l = 0; l < n; ++l) {
    			auto Rnl = results[i].get();
    			auto Nnl = electronStates(n, l);
    			for (std::size_t k = 0; k < x.size(); ++k) density[k] += Nnl*Rnl[k]*Rnl[k];
    			++i;
    		}
		}
	}
	else 
#endif
	{
		for (int n = 1; n <= nmax; ++n) {
    		for (int l = 0; l < n; ++l) {
    			auto enl = energyLevel(n, l);
                double lambda = l + 0.5;
    			auto Rnl = waveFunction(enl, lambda, x);
    			auto Nnl = electronStates(n, l);
    			for (std::size_t k = 0; k < x.size(); ++k) density[k] += Nnl*Rnl[k]*Rnl[k];
    		}
    	}

        if (useContinuous){
            for (std::size_t k = 0; k < x.size(); ++k) {
                density[k] += electronDensityContinuous(x[k])*(4.0 * M_PI * pow(x[k]*r0,2.0));
            }
        }
	}
    if (nUpdate > 0) {
        // 4. mixing with previous
        auto& oldDensity = *densityInterpolation;
        for (std::size_t k = 0; k < x.size(); ++k) density[k] = mixing*oldDensity(u[k]) + (1.0 - mixing)*density[k];
    }
    densityInterpolation = new Spline(u, density);
    nUpdate++;
}

void Atom::reset(double _V, double _T, double _Z, int _nmax) {
	if (_V > 0)    volume = _V;
	if (_T >= 0)   temperature = _T;
	if (_Z > 0)    Zcharge = _Z;
	if (_nmax > 0) nmax = _nmax;

	r0 = std::pow(3.0*volume/(4.0*M_PI), 1.0/3.0);

	ChemicalPotential MTF;
    MTF.setZ(Zcharge);
    chemPot = MTF(volume, temperature);
    bool chemPotReady = true;

    ElectronDensity DTF;
    DTF.setVTZ(volume, temperature, Zcharge);
	std::vector<double> x(DTFsize);
	std::vector<double> u(DTFsize);
	double umin = 1e-3;
	double umax = 1.0;

	for (int i = 0; i < DTFsize; ++i) {
		u[i] = umin + i*(umax - umin)/(DTFsize - 1);
		x[i] = u[i]*u[i];
	}

	auto d = DTF(x);
	for (int i = 0; i < DTFsize; ++i) {
		d[i] *= 4.0*M_PI*r0*r0*x[i]*x[i];
	}

    densityInterpolation = new Spline(u, d);

    // prepare levels
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

}
}