#include <cmath>
#include <iostream>
#include <average-atom-toolkit/atom/semiclassic.h>
#include <iomanip>
#include <fstream>


// linspace c++ analogue
//template <typename T>
//std::vector<T> linspace(T a, T b, size_t N) {
//    T h = (b - a) / static_cast<T>(N-1);
//    std::vector<T> xs(N);
//    typename std::vector<T>::iterator x;
//    T val;
//    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
//        *x = val;
//    return xs;
//}

int main(){
    std::cout << "Test" << '\n';
    const double k = 36.75, N_a = 6.022140857E+23, avol = std::pow(5.2917720859e-9,3);
    aatk::atom::SemiclassicAtom cell;
    double Te,Tet, rho, A, Z, r0;
    const int n_max = 20;
    bool Continuous = true;
    const double E_harte = 27.21;
    double sigma = 0.0 / E_harte;
    int    meshSize = 1600;
    double tolerance = 1.e-6;
    // Initial parameters of Temperature Te (ev) and density rho (g/cm^3)
    // Tet is reduced temperature
    Te = 90;//114.0;//77.0;//27 wf;//100.0;
    Tet = k * Te * 1E-3;
    rho = 2.698900;//0.1;//
    // Setting material properties: A - atomic mass ,Z  - number of element in periodic table
    A = 58.69;//55.845;//A = 56.0 Fe//107E0 Ag//197E0 Au// Al 26.98 // Sn 118.7 // W 183.84 // Ni 58.69
    Z = 28.00;//26.0;//Z = 26.0; Fe//47E0 Ag//79.0E0 Au// Al 13.0 // Sn 50 // W 74 // Ni 28
    r0 = 1.3882E0 * std::pow(A / rho,1.0 / 3.0);
    double V = A/(N_a*rho*avol);
    // setting Cell initial parameters
    cell.reset(V,Tet,Z,tolerance,meshSize ,n_max, Continuous);

//    //x grid
//    const int N = 1000;
//    std::vector<double> x;
//    double x_mesh[N];
//    x = linspace(1.0/N, 1.0, N);

//    for(int i = 0; i < x.size(); i++){
//        x_mesh[i] = x[i] * x[i];
//    }

    int iterations = 1;
    const int iterations_max = 50;
    double E_iner_0 = 0.76874512 * std::pow(Z,7.0/3.0);
    // Test
    std::cout << "self - consistent cycle:" << std::endl;
    std::cout << "Z = " << Z << " T = " << Te << " eV" << std::endl;
    const double dE_max = 1e-5;
    double dE = 2 * dE_max , E_cur = 1;//cell.energyFull();


    while((abs(dE) > dE_max && iterations < iterations_max) || iterations <= 2){
        cell.update( 0.75);

        std::cout << "Iteration № " << iterations << std::endl;
        std::cout << "N = "<< std::setprecision(8) << cell.electronStatesDiscrete() + cell.electronStatesContinuous() << " " << "Mu = "<< cell.chemicalPotential() << std::endl;
        std::cout << "n_max = " << cell.discreteLevelsNumber() << std::endl;
        std::cout << "Inner energy  =" << cell.internalEnergy() + E_iner_0  << std::endl;// E_iner_0 + 1.5 * Tet
        std::cout << "Entropy =" << cell.entropy() << std::endl;// E_iner_0 + 1.5 * Tet
        std::cout << "Pressure =" << cell.pressure() << std::endl;// E_iner_0 + 1.5 * Tet


        double E_new = cell.energyFull();
        dE = abs(E_cur - E_new) / -(E_cur + E_new);
        E_cur = E_new;
        std::cout << "dE= " <<  dE << std::endl;
        iterations++;
    }

//    std::ofstream out;
//    out.open("sc.txt"); // Path to file
//    if (out.good()){
//        std::cout << "stream is Ok" << std::endl;
//    }
//
//    for(auto x_cur : x){
//        out << x_cur*x_cur << " " << cell.electronDensity(x_cur*x_cur)<< std::endl;
//    }

    std::cout << "test convergence completed" <<  std::endl;
    std::cout << "E0 = "<< cell.boundaryEnergyValue() << std::endl;
    std::cout << "Nu = "<< -cell.chemicalPotential() / Tet << std::endl;
    double Z_calculated = cell.electronStatesDiscrete() + cell.electronStatesContinuous() ;
    std::cout << "Z = " << Z <<" " <<Z_calculated <<std::endl;
    std::cout << "dZ = "<< std::setprecision(8)<< Z - Z_calculated<< std::endl;

    return 0;
}
