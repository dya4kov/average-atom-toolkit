#pragma once
#include <vector>

namespace AATools {
namespace TF {

class Table {
public:
    Table();
    double operator()(const int& v, const int& t) const;
    double operator()(const double& lgV, const double& lgT) const;
private:
    const int    vSize;
    const int    tSize;
    const double lgV0;
    const double lgT0;
    const double lgVstep;
    const double lgTstep;
    double zeroFit(const double& lgV) const;
    static std::vector<double> table;
};

}
}