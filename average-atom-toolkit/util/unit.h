#pragma once
#include <string>
#include <iomanip>
#include <sstream>

namespace aatk {
namespace util {
namespace unit {

double  m_value = 1.0;
double kg_value = 1.0;
double  s_value = 1.0;

class Unit {
public:
    Unit() : l(0), m(0), t(0), si_value(1.0) {}
    Unit(int _l, int _m, int _t) : 
           l(_l),  m(_m),  t(_t), si_value(1.0) {}
    Unit(int _l, int _m, int _t, double val) : 
           l(_l),  m(_m),  t(_t) {
        double factor = 1.0;
        if (l < 0) for (int il = l; il < 0; ++il) factor *=  m_value;
        if (m < 0) for (int im = m; im < 0; ++im) factor *= kg_value;
        if (t < 0) for (int it = t; it < 0; ++it) factor *=  s_value;
        if (l > 0) for (int il = l; il > 0; --il) factor /=  m_value;
        if (m > 0) for (int im = m; im > 0; --im) factor /= kg_value;
        if (t > 0) for (int it = t; it > 0; --it) factor /=  s_value;
        si_value = factor*val;
    }
    Unit& operator=(const Unit& u) {
        l = u.l;
        m = u.m;
        t = u.t;
        si_value = u.si_value;
        return *this;
    }
    double value() const {
        double factor = 1.0;
        if (l < 0) for (int il = l; il < 0; ++il) factor /=  m_value;
        if (m < 0) for (int im = m; im < 0; ++im) factor /= kg_value;
        if (t < 0) for (int it = t; it < 0; ++it) factor /=  s_value;
        if (l > 0) for (int il = l; il > 0; --il) factor *=  m_value;
        if (m > 0) for (int im = m; im > 0; --im) factor *= kg_value;
        if (t > 0) for (int it = t; it > 0; --it) factor *=  s_value;
        return factor*si_value;
    }
    std::string toString() const {
        std::stringstream ss;
        ss << "<";
        ss << std::scientific << value();
        ss << ", ";
        ss << "L" << l << "M" << m << "T" << t << ">";
        return ss.str();
    }
private:
    friend Unit operator*(const Unit& u1, const Unit& u2);
    friend Unit operator/(const Unit& u1, const Unit& u2);
    friend Unit operator*(const double& val, const Unit& u);
    friend Unit operator/(const double& val, const Unit& u);
    friend void setUnitLength(const Unit& length);
    friend void setUnitMass(const Unit& mass);
    friend void setUnitTime(const Unit& time);
    friend bool similar(const Unit& u1, const Unit& u2);

    int l;
    int m;
    int t;
    double si_value;
};

Unit operator*(const Unit& u1, const Unit& u2) {
    return Unit(
        u1.l + u2.l,
        u1.m + u2.m,
        u1.t + u2.t,
        u1.value()*u2.value()
    );
}

Unit operator/(const Unit& u1, const Unit& u2) {
    return Unit(
        u1.l - u2.l,
        u1.m - u2.m,
        u1.t - u2.t,
        u1.value()/u2.value()
    );
}

Unit operator*(const double& val, const Unit& u) {
    return Unit(
        u.l,
        u.m,
        u.t,
        val*u.value()
    );
}

Unit operator/(const double& val, const Unit& u) {
    return Unit(
        -u.l,
        -u.m,
        -u.t,
        val/u.value()
    );
}

void setUnitLength(const Unit& length) {
    bool validUnit = true;
    validUnit = validUnit && length.l == 1;
    validUnit = validUnit && length.m == 0;
    validUnit = validUnit && length.t == 0;
    if (validUnit)
       m_value = m_value/length.value();
}

void setUnitMass(const Unit& mass) {
    bool validUnit = true;
    validUnit = validUnit && mass.l == 0;
    validUnit = validUnit && mass.m == 1;
    validUnit = validUnit && mass.t == 0;
    if (validUnit)
        kg_value = kg_value/mass.value();
}

void setUnitTime(const Unit& time) {
    bool validUnit = true;
    validUnit = validUnit && time.l == 0;
    validUnit = validUnit && time.m == 0;
    validUnit = validUnit && time.t == 1;
    if (validUnit)
        s_value = s_value/time.value();
}

bool similar(const Unit& u1, const Unit& u2) {
    bool result = true;
    result = result && u1.l == u2.l;
    result = result && u1.m == u2.m;
    result = result && u1.t == u2.t;
    return result;
}

// basic units
Unit _1_ = Unit(0, 0, 0);
Unit  m  = Unit(1, 0, 0);
Unit kg  = Unit(0, 1, 0);
Unit  s  = Unit(0, 0, 1);

Unit Avogadro  = 6.022140857e+23*_1_;
Unit qElectron = 1.60217657e-19*_1_;
Unit kBoltzman = 1.38064852e-23*_1_;

// length units
                   Unit  m2 =  m*m ; Unit  m3 =  m2*m ;
Unit km = 1.e+3*m; Unit km2 = km*km; Unit km3 = km2*km;
Unit cm = 1.e-2*m; Unit cm2 = cm*cm; Unit cm3 = cm2*cm;
Unit mm = 1.e-3*m; Unit mm2 = mm*mm; Unit mm3 = mm2*mm;
Unit um = 1.e-6*m; Unit um2 = um*um; Unit um3 = um2*um;
Unit nm = 1.e-9*m; Unit nm2 = nm*nm; Unit nm3 = nm2*nm;

Unit rBohr = 5.2917720859e-11*m;
Unit aVol  = rBohr*rBohr*rBohr;

// mass units
Unit  g = 1.e-3*kg;
Unit mg = 1.e-6*kg;
Unit ug = 1.e-9*kg;

Unit eMass = 9.10938291e-31*kg;
Unit aMass = 1.66054021e-27*kg;

// time units
Unit ms = 1.e-3*s;
Unit us = 1.e-6*s;
Unit ns = 1.e-9*s;

// energy
Unit  J = kg*m2/(s*s);
Unit GJ = 1.e+9*J;
Unit MJ = 1.e+6*J;
Unit kJ = 1.e+3*J;
Unit mJ = 1.e-3*J;
Unit uJ = 1.e-6*J;
Unit nJ = 1.e-9*J;

Unit eV = qElectron*J;
Unit Ry = 13.605693009*eV;
Unit hartree = 2.0*Ry;
Unit K  = kBoltzman*J;

// pressure
Unit Pa  = J/m3;
Unit GPa = 1.e+9*Pa;
Unit MPa = 1.e+6*Pa;
Unit kPa = 1.e+3*Pa;
Unit mPa = 1.e-3*Pa;
Unit uPa = 1.e-6*Pa;
Unit nPa = 1.e-9*Pa;

Unit  bar = 1.e+5*Pa;
Unit kbar = 1.e+3*bar;
Unit Mbar = 1.e+6*bar;

// misc
Unit c    = 2.99792458e+8*m/s;
Unit hbar = 1.0545718e-34*J*s;

Unit aTime = hbar/hartree;

}
}
}