#include <iostream>
#include <average-atom-toolkit/semiclassic/atom.h>
#include <iomanip>
#include <fstream>


// linspace c++ analogue
template <typename T>
std::vector<T> linspace(T a, T b, size_t N) {
    T h = (b - a) / static_cast<T>(N-1);
    std::vector<T> xs(N);
    typename std::vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}

int main(){
    std::cout << "Test" << '\n';
    const double k = 36.75, N_a = 6.022140857E+23, avol = pow(5.2917720859e-9,3);
    aatk::semiclassic::Atom cell;
    double Te,Tet, rho, A, Z, r0;
    const int n_max = 8;

    // Initial parameters of Temperature Te (ev) and density rho (g/cm^3)
    // Tet is reduced temperature
    Te = 100;
    Tet = k * Te * 1E-3;
    rho = 0.1;

    // Setting material properties: A - atomic mass ,Z  - number of element in periodic table
    A = 197.0;//55.845;//A = 56.0 Fe//107E0 Ag//197E0 Au
    Z = 79.0;//26.0;//Z = 26.0; Fe//47E0 Ag//79.0E0 Au
    r0 = 1.3882E0 * pow(A / rho,1.0 / 3.0);
    double V = A/(N_a*rho*avol);

    // setting Cell initial parameters
    cell.reset(V,Tet,Z,n_max);

    //x grid
    const int N = 1000;
    std::vector<double> x;
    double x_mesh[N];
    x = linspace(1.0/N, 1.0, N);

    for(int i = 0; i < x.size(); i++){
        x_mesh[i] = x[i] * x[i];
    }

    int iterations = 1;
    const int iterations_max = 73;

    // Test
    std::cout << "self - consistent cycle:" << std::endl;

    while(iterations <= iterations_max){
        cell.update(x_mesh,N, 0.75);

        std::cout << "Iteration № " << iterations << std::endl;
        std::cout << "N = "<< std::setprecision(8) << cell.electronStates() << " " << "Mu = "<< cell.M() << std::endl;
        std::cout << "n_max = " << cell.N_max() << std::endl;
        iterations++;
    }

    std::ofstream out;
    out.open("/Users/SashaP/soft/test_files/conv_cont_atom.txt"); //
    if (out.good()){
        std::cout << "stream is Ok" << std::endl;
    }

    for(auto x_cur : x){
        out << x_cur*x_cur << " " << cell.electronDensity(x_cur*x_cur)<< std::endl;
//      cout << x_cur*x_cur  << " " << cell.electronDensity(x_cur*x_cur )<< endl;
    }

    std::cout << "test convergence completed" <<  std::endl;
    std::cout << "Nu = "<< -cell.M() / Tet << std::endl;
    double Z_calculated = cell.electronStates(cell.M()) + cell.electronStatesContinuous(cell.M()) ;
    std::cout << "Z = " << Z <<" " <<Z_calculated <<std::endl;
    std::cout << "dZ = "<< std::setprecision(8)<< Z - Z_calculated<< std::endl;
}
