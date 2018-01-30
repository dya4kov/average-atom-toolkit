#pragma once
#include <cmath>

#define FTTF_TABLE_VSIZE 181
#define FTTF_TABLE_TSIZE 201

namespace AATools {
namespace TF {    

struct PotTable {
    PotTable();
    const unsigned vSize;
    const unsigned tSize;
    const double lgV0;
    const double lgT0;
    const double lgVstep;
    const double lgTstep;
    double operator()(const double& lgV, const double& lgT) {
        int v, t;
        v = (int) std::floor((lgV - lgV0)/lgVstep);
        t = (int) std::floor((lgT - lgT0)/lgTstep);
        if (t < 0) return zeroFit(lgV);
        double result = 0.0;
        if ((v >= 0 && v < vSize) || (t < tSize)) {
            double f00, f01, f10, f11;
            f00 = table[v     +       t*FTTF_TABLE_VSIZE];
            f01 = table[v     + (t + 1)*FTTF_TABLE_VSIZE];
            f10 = table[v + 1 +       t*FTTF_TABLE_VSIZE];
            f11 = table[v + 1 + (t + 1)*FTTF_TABLE_VSIZE];
            int sign = 0;
            sign = (sign == 0 && f00 > 0.0 && f01 > 0.0 && f10 > 0.0 && f11 > 0.0) ?  1 : 0;
            sign = (sign == 0 && f00 < 0.0 && f01 < 0.0 && f10 < 0.0 && f11 < 0.0) ? -1 : 0;
            if (sign != 0) {
            	f00 = std::log10(std::abs(f00));
                f01 = std::log10(std::abs(f01));
                f10 = std::log10(std::abs(f10));
                f11 = std::log10(std::abs(f11));
            }
            double V0 = (lgV - lgV0)/lgVstep - 1.0*v;
            double T0 = (lgT - lgT0)/lgTstep - 1.0*t;
            double a[2][2];
            a[0][0] = f00;
            a[1][0] = f10 - f00;
            a[0][1] = f01 - f00;
            a[1][1] = f11 + f00 - (f10 + f01);
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    result += a[i][j]*std::pow(V0, i)*std::pow(T0, j);
                }
            }
            if (sign != 0) result = sign*std::pow(10.0, result);
        }
        return result;
    }
    double operator()(const unsigned& v, const unsigned& t) { return table[v + t*FTTF_TABLE_VSIZE]; }
private:
	double zeroFit(const double& lgV) {
		const double B0 =  0.648742997083556;
        const double B1 = -0.704984628856768;
        const double B2 = -0.0224226496439102;
        const double B3 = -0.00419385235723519;
        const double B4 = -3.75915351702641E-5;
        const double B5 =  3.94764845762704E-5;
        const double B6 =  5.4828018180471E-7;
        const double B7 = -1.49964096611993E-7;

        double lgV1 = lgV;
        double lgV2 = lgV1*lgV1;
        double lgV3 = lgV1*lgV2;
        double lgV4 = lgV2*lgV2;
        double lgV5 = lgV3*lgV2;
        double lgV6 = lgV3*lgV3;
        double lgV7 = lgV4*lgV3;

        double lgMu = B0 + B1*lgV1 + B2*lgV2 + B3*lgV3 + B4*lgV4 + B5*lgV5 + B6*lgV6 + B7*lgV7;
        return std::pow(10.0, lgMu);
	}
    double table[FTTF_TABLE_VSIZE*FTTF_TABLE_TSIZE];
};

}
}