#include <iostream>
#include <fstream>
#include <string>

#define FTTF_TABLE_VSIZE 181
#define FTTF_TABLE_TSIZE 201

int main(int argc, char *argv[]) {
	std::string ipath = argv[1];
	std::string opath = argv[2];
	std::ifstream in(ipath + "/table.dat");
	std::ofstream out(opath + "/table.cxx");
	out << "#include <average-atom-tools/thomas-fermi/atom/table.h>" << std::endl;
	out << "AATools::TF::PotTable::PotTable() : ";
	out << "vSize(" << FTTF_TABLE_VSIZE << "), "
	    << "tSize(" << FTTF_TABLE_TSIZE << "), "
	    << "lgV0(-10.0), "
	    << "lgT0(-10.0), "
	    << "lgVstep(0.1), "
	    << "lgTstep(0.1) { " << std::endl;
	double value;
	for (int t = 0; t < FTTF_TABLE_TSIZE; ++t) {
		for (int v = 0; v < FTTF_TABLE_VSIZE; ++v) {
			in >> value;
			out << "table[" << t*FTTF_TABLE_VSIZE + v << "] = ";
			out.setf(std::ios::scientific | std::ios::left );
			out.width(20);
    		out.precision(12);
    		out << value << ";" << std::endl;
		}
	}
	out << "}";
	in.close();
	out.close();
}