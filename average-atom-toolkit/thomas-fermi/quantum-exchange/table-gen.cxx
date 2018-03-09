#include <fstream>
#include <string>
#include <cstddef>

void gen_table_cxx(std::string opath, double* data, std::size_t vSize, std::size_t tSize) {
    opath += "/table.cxx";
	std::ofstream out(opath.c_str());
	out << "#include <cmath>" << std::endl;
	out << "#include <average-atom-toolkit/thomas-fermi/quantum-exchange/table.h>" << std::endl;
	out << "aatk::TF::QE::Table::Table() : ";
	out << "vSize(" << vSize << "), "
	    << "tSize(" << tSize << "), "
	    << "lgV0(-10.0), "
	    << "lgT0(-10.0), "
	    << "lgVstep(0.1), "
	    << "lgTstep(0.1) {} " << std::endl;

	out << "double aatk::TF::QE::Table::operator()(const int& v, const int& t) const {  " << std::endl;
	out << "	return table[v + t*vSize]; " << std::endl;
	out << "} " << std::endl;

	out << "double aatk::TF::QE::Table::operator()(const double& lgV, const double& lgT) const {" << std::endl;
    out << "    int v, t;" << std::endl;
    out << "    v = (int) std::floor((lgV - lgV0)/lgVstep);" << std::endl;
    out << "    t = (int) std::floor((lgT - lgT0)/lgTstep);" << std::endl;
    out << "    if (t < 0) return zeroFit(lgV);" << std::endl;
    out << "    double result = 0.0;" << std::endl;
    out << "    if ((v >= 0 && v < vSize) || (t < tSize)) {" << std::endl;
    out << "        double f00, f01, f10, f11;" << std::endl;
    out << "        f00 = -table[v     +       t*vSize];" << std::endl;
    out << "        f01 = -table[v     + (t + 1)*vSize];" << std::endl;
    out << "        f10 = -table[v + 1 +       t*vSize];" << std::endl;
    out << "        f11 = -table[v + 1 + (t + 1)*vSize];" << std::endl;
    out << "        f00 = std::log10(std::abs(f00));" << std::endl;
    out << "        f01 = std::log10(std::abs(f01));" << std::endl;
    out << "        f10 = std::log10(std::abs(f10));" << std::endl;
    out << "        f11 = std::log10(std::abs(f11));" << std::endl;
    out << "        double V0 = (lgV - lgV0)/lgVstep - 1.0*v;" << std::endl;
    out << "        double T0 = (lgT - lgT0)/lgTstep - 1.0*t;" << std::endl;
    out << "        double a[2][2];" << std::endl;
    out << "        a[0][0] = f00;" << std::endl;
    out << "        a[1][0] = f10 - f00;" << std::endl;
    out << "        a[0][1] = f01 - f00;" << std::endl;
    out << "        a[1][1] = f11 + f00 - (f10 + f01);" << std::endl;
    out << "        for (int i = 0; i < 2; ++i) {" << std::endl;
    out << "            for (int j = 0; j < 2; ++j) {" << std::endl;
    out << "                result += a[i][j]*std::pow(V0, i)*std::pow(T0, j);" << std::endl;
    out << "            }" << std::endl;
    out << "        }" << std::endl;
    out << "        result = -std::pow(10.0, result);" << std::endl;
    out << "    }" << std::endl;
    out << "    return result;" << std::endl;
    out << "}" << std::endl;

    out << std::endl;

    out << "double aatk::TF::QE::Table::zeroFit(const double& lgV) const {" << std::endl;
	out << "	const double B0 =  1.20752290577393;" << std::endl;
    out << "    const double B1 = -0.338439954152781;" << std::endl;
    out << "    const double B2 = -0.00667391360555085;" << std::endl;
    out << "    const double B3 = -0.00202025760222042;" << std::endl;
    out << "    const double B4 = -1.77130362887407E-4;" << std::endl;
    out << "    const double B5 =  8.67013549505259E-6;" << std::endl;
    out << "    const double B6 =  1.89781271548473E-6;" << std::endl;
    out << "    const double B7 =  6.62246701639077E-8;" << std::endl;

    out << "    double lgV1 = lgV;" << std::endl;
    out << "    double lgV2 = lgV1*lgV1;" << std::endl;
    out << "    double lgV3 = lgV1*lgV2;" << std::endl;
    out << "    double lgV4 = lgV2*lgV2;" << std::endl;
    out << "    double lgV5 = lgV3*lgV2;" << std::endl;
    out << "    double lgV6 = lgV3*lgV3;" << std::endl;
    out << "    double lgV7 = lgV4*lgV3;" << std::endl;

    out << "    double lgPsi = B0 + B1*lgV1 + B2*lgV2 + B3*lgV3 + B4*lgV4 + B5*lgV5 + B6*lgV6 + B7*lgV7;" << std::endl;
    out << "    return -std::pow(10.0, lgPsi);" << std::endl;
	out << "}" << std::endl;

    out << "std::vector<double> aatk::TF::QE::Table::table({" << std::endl;

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