#include <cmath>
#include <iostream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>

#include <numeric-toolkit/specfunc/fermi-dirac/complete.h>
#include <numeric-toolkit/specfunc/fermi-dirac/incomplete.h>


#include <average-atom-toolkit/atom/semiclassic.h>

namespace aatk {
namespace atom {

double energyDensityContinuousFunc(double x, void * params){
    SemiclassicAtom * tempAtom = (SemiclassicAtom *)params;

    numtk::specfunc::FermiDirac<numtk::specfunc::FD::ThreeHalf>        FD_ThreeHalf;
    numtk::specfunc::FermiDiracInc<numtk::specfunc::FDI::ThreeHalfInc> FD_ThreeHalf_Inc;
    numtk::specfunc::FermiDirac<numtk::specfunc::FD::Half> FD_Half;
    numtk::specfunc::FermiDiracInc<numtk::specfunc::FDI::HalfInc> FD_Half_Inc;
    const double E0 = tempAtom->boundaryEnergyValue(); 
    const double V_r = tempAtom->U(x);
    const double mu = tempAtom->chemicalPotential(); 
    const double T = tempAtom->temperature();
    const double factor = std::sqrt(2) * std::pow( T,3.0 / 2.0) / (std::pow(M_PI, 2) );
    const double y0 = (V_r + E0)/T;
    double result;

    if (y0 <= 0){
        result = T * FD_ThreeHalf((V_r + mu)/T) -  V_r * FD_Half((V_r + mu)/T) ;
    }
    else{
        result = T * (FD_ThreeHalf((V_r + mu)/T) - FD_ThreeHalf_Inc((V_r + mu)/T,y0)) -
                V_r * (FD_Half((V_r + mu)/T) - FD_Half_Inc((V_r + mu)/T,y0));
    }

    //return result * factor;

    return  x * x * result * factor;
}
double SemiclassicAtom::energyContinuous() {
    double result, error;
    SemiclassicAtom * tempAtom = (SemiclassicAtom *)this;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_function Func;
    Func.function = &energyDensityContinuousFunc;
    Func.params = tempAtom;

    gsl_integration_qags (&Func,1e-6,1,tolerance ,tolerance,1000,w,&result, &error); // epsabs = 1e-6
    gsl_integration_workspace_free (w);

    return 4 * M_PI * result * std::pow(r0,3);
}

double SemiclassicAtom::energyFull(){
    double result,integral, error;
    double Enl;
    double Nnl;
    result = 0.0;
    if (nUpdate == 0) {
        for (int n = 1; n <= nmax; ++n) {
            for (int l = 0; l < n; ++l) {
                evaluate_energy_level(n, l);
            }
        }
        evaluate_boundary_energy(); // remove me ?
    }


    for (int n = 1; n <= nmax; n++){
        for (int l = 0; l < n; l++ ){
            Enl = energyLevel(n, l);
            Nnl = electronStatesDiscrete(n, l);
            result += Enl * Nnl;        
        }
    }

    if (useContinuous){
        result += energyContinuous();
    }

    return result;
}

double internalEnergyFunc (double x, void * params){
    SemiclassicAtom * tempAtom = (SemiclassicAtom *)params;
    double result;
    double Z  = tempAtom->Znucleus();
    double r0 = tempAtom->radius();
    double U  = tempAtom->U(x);
    double n  = tempAtom->electronDensity(x) ;// / 4.0 * M_PI * std::pow(x*r0,2.0);
    //temp_cell->electronDensityDiscrete(x) + temp_cell->electronDensityContinuous(x);
    result = n * (  U -  Z / (r0*x) );

    return  result;
}

double SemiclassicAtom::internalEnergy(){
    double result = 0.0;
    double integral = 0.0;
    double error;
    SemiclassicAtom * tempAtom = (SemiclassicAtom *)this;

    result += energyFull();

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_function Func;
    Func.function = &internalEnergyFunc;
    Func.params = tempAtom;

    gsl_integration_qags (&Func,1e-6,1,
                          tolerance ,tolerance * 10,1000,w,&integral, &error);
    gsl_integration_workspace_free (w);
    result += 0.5 * r0 * integral;//* M_PI * std::pow(r0,3)


    return result;
}

double entropyFunc(double x, void * params){
    SemiclassicAtom * tempAtom = (SemiclassicAtom *)params;
    numtk::specfunc::FermiDirac<numtk::specfunc::FD::ThreeHalf>        FD_ThreeHalf;
    numtk::specfunc::FermiDiracInc<numtk::specfunc::FDI::ThreeHalfInc> FD_ThreeHalf_Inc;

    double mu     = tempAtom->chemicalPotential();
    double E0     = tempAtom->boundaryEnergyValue();
    double T      = tempAtom->temperature();
    double U      = tempAtom->U(x);
    double n_c    = tempAtom->electronDensityContinuous(x);
    double result = 0.0;
    double y0     = (E0 + U) / T;


    result -= 5.0 * std::sqrt(2.0) * std::pow(T, 3.0 / 2.0) / (3 * M_PI * M_PI) * FD_ThreeHalf((U + mu) / T);

    if ( y0 > 0){
        result += 2.0 * std::sqrt(2.0) / (3 * M_PI * M_PI) * std::log(1.0 + std::exp( (mu - E0 ) / T)) *
                  std::pow(E0 + U, 3.0/2.0);

        result += 5.0 * std::sqrt(2.0) *
                std::pow(T, 3.0 / 2.0) / (3 * M_PI * M_PI) * FD_ThreeHalf_Inc((U + mu) / T, y0);
    }

    result += (1 / T) * U * n_c;

    return  result * x * x;
}

double SemiclassicAtom::entropy(){
    double result = 0.0;
    double integral = 0.0;
    double error, Enl, Nnl;
    //double T = temperature;
    SemiclassicAtom * tempAtom = (SemiclassicAtom *)this;
    result = -M * Z / T; // Z or N ?

    if (nUpdate == 0) {
        for (int n = 1; n <= nmax; ++n) {
            for (int l = 0; l < n; ++l) {
                evaluate_energy_level(n, l);
            }
        }
    evaluate_boundary_energy(); // remove me ?
    }

    for (int n = 1; n <= nmax; n++){
        for (int l = 0; l < n; l++ ){
            Enl = energyLevel(n, l);
            Nnl = electronStatesDiscrete(n,l);
            double exp_ind = (M - Enl) / T;
            if ((useContinuous && Enl < boundaryEnergy) || !useContinuous ){
                if (exp_ind > 50.0){
                result += 2.0 * (2 * l + 1) * ( (M - Enl ) / T) ;
            }
            else{
                result += 2.0 * (2 * l + 1) * std::log(1.0 + std::exp( (M - Enl ) / T)) ;

            }
            }

            result += Nnl * Enl / T;
        }
    }

    if (useContinuous){
        gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
        gsl_function Func;
        Func.function = &entropyFunc;
        Func.params = tempAtom;

        gsl_integration_qags (&Func,1e-6,1,
                              tolerance ,tolerance * 10,1000,w,&integral, &error);
        gsl_integration_workspace_free (w);
        result -= 4 * M_PI * std::pow(r0,3) * integral;
    }


    return result;
}



double SemiclassicAtom::pressure() {
    numtk::specfunc::FermiDirac<numtk::specfunc::FD::ThreeHalf>        FD_ThreeHalf;
    numtk::specfunc::FermiDiracInc<numtk::specfunc::FDI::ThreeHalfInc> FD_ThreeHalf_Inc;
    double result = 0.0;
    double Nnl,Enl,Rnl_div_1,Rnl_div_2,Rnl = 0.0;
    double h = 1.0 / (meshSize - 1);
    int index_end = meshSize - 1;
//    std::vector<double> x(meshSize); //{ 1.0 - 2*h,1.0 - h,1.0};
    std::vector<double> R_vec(mesh.size());
    double factor_d = -1.0 / (8.0 * M_PI * std::pow(r0,3.0));
    double factor_c = std::pow(2 * T, 5.0/2.0) / (6 * M_PI * M_PI);
    double lambda;

    if (nUpdate == 0) {
        for (int n = 1; n <= nmax; ++n) {
            for (int l = 0; l < n; ++l) {
                evaluate_energy_level(n, l);
            }
        }
        evaluate_boundary_energy(); // remove me ?
    }

    double y0 = (boundaryEnergy + M)/T;
    

    for (int n = 1; n <= nmax; n++){
        for (int l = 0; l < n; l++ ){
            lambda = l + 0.5;
            Enl = energyLevel(n, l);
            Nnl = electronStatesDiscrete(n, l);

            if ((useContinuous && Enl < boundaryEnergy)|| !useContinuous ){
                waveFunction(mesh.data(), R_vec.data(), mesh.size(), Enl, lambda);
                Rnl = R_vec[index_end];

                // f(n-2) - 4 f(n-1) + 3 f(n) /2h
                Rnl_div_1 = (R_vec[index_end-2] - 4 * R_vec[index_end-1] + 3 * R_vec[index_end] ) / (2 * h);

                Rnl_div_2 = - 2*(Enl - l * (l + 1) / (2 * r0 * r0))* Rnl;
                result += factor_d * Nnl *
                        (Rnl * Rnl_div_1  + r0 * Rnl * Rnl_div_2 - r0 * Rnl_div_1 * Rnl_div_1 );

//                 if (std::abs(result) > 0.0 || Rnl != 0.0){
//                     std::cout << "n = " << n << " l = " << l << std::endl;
//                     std::cout << "Rnl = " << Rnl << " Rnl_div_1 = " << Rnl_div_1 << std::endl;
//                 }
            }
        }
    }
  //  std::cout << "Pressure_discrete = "<< result <<std::endl;

    if (y0 <= 0){
        result += factor_c * FD_ThreeHalf(M/T);
    }else{
        result += factor_c * (FD_ThreeHalf(M/T) - FD_ThreeHalf_Inc(M/T, y0));
    }


    return result;
}




}
}

