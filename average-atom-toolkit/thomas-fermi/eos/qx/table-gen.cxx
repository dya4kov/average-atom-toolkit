#include <fstream>
#include <string>
#include <cstddef>

void gen_table_cxx(std::string opath, double* data, std::size_t vSize, std::size_t tSize) {
    opath += "/table.cxx";
    std::ofstream out(opath.c_str());
    out << "#include <average-atom-toolkit/thomas-fermi/eos/qx/chemical-potential.h>" << std::endl;
    out << "std::vector<double> aatk::TF::qx::ChemicalPotential::table({" << std::endl;

    out.setf(std::ios::scientific | std::ios::left );
    out.width(20);
    out.precision(12);

    for (int t = 0; t < tSize - 1; ++t) {
    for (int v = 0; v < vSize; ++v) {

    out << data[t*vSize + v] << "," << std::endl;

    }
    }

    for (int v = 0; v < vSize - 1; ++v) {

    out << data[(tSize - 1)*vSize + v] << "," << std::endl;

    }

    out << data[tSize*vSize - 1] << std::endl;

    out << "});" << std::endl;

    out << std::endl;

    out.close();
}

int main(int argc, char *argv[]) {
    std::string ipath = argv[1];
    std::string opath = argv[2];
    ipath += "/table";
    std::ifstream in(ipath.c_str());

    int vSize; in.read((char*) &vSize, sizeof(int));
    int tSize; in.read((char*) &tSize, sizeof(int));
    double* data = new double[vSize*tSize];

    in.read((char*) data, sizeof(double)*vSize*tSize);

    gen_table_cxx(opath, data, vSize, tSize);

    in.close();
}