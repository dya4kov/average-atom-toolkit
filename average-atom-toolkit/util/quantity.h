#pragma once
#include <string>
#include <iostream>
#include <cstddef>
#include <iomanip>
#include <sstream>

#include <average-atom-toolkit/util/unit.h>

namespace aatk {
namespace util {

class Quantity {
public:
    Quantity() : un(0,0,0), val(0.0) {}
    Quantity(const double& v, const unit::Unit& u) {
        val = v;
        un = u;
    }
    Quantity(const unit::Unit& q) {
        val = q.value();
        un = 1.0/val*q;
    }
    Quantity& operator=(const Quantity& q) {
        val = q.val;
        un = q.un;
        return *this;
    }
    Quantity& operator=(const unit::Unit& q) {
        val = q.value();
        un = 1.0/val*q;
        return *this;
    }
    Quantity& operator*=(const double& value) {
        val *= value;
    }
    Quantity& operator/=(const double& value) {
        val /= value;
    }
    Quantity& operator+=(const Quantity& q) {
        if (similar(un, q.un))
            val += q(un);
        else {
            std::cerr << "runtime error: cannot sum quantities " 
                      << toString() << " and " << q.toString()
                      << std::endl;
            exit(0);
        }
        return *this;
    }
    Quantity& operator-=(const Quantity& q) {
        if (similar(un, q.un))
            val -= q(un);
        else {
            std::cerr << "runtime error: cannot sum quantities " 
                      << toString() << " and " << q.toString()
                      << std::endl;
            exit(0);
        }
        return *this;
    }
    double operator()() const {
        return val;
    }
    double operator()(const unit::Unit& u) const {
    	double result;
    	if (similar(u, un))
    		result = val*un.value()/u.value();
    	else {
    		std::cerr << "runtime error: cannot express unit " 
                  << un.toString() << " through " << u.toString()
                  << std::endl;
        	exit(0);
    	}
    	return result;
    }

    unit::Unit unit() const { return un; }

    std::string toString() const {
        std::stringstream ss;
        ss << "<";
        ss << std::scientific << operator()();
        ss << ", ";
        ss << un.toString() << ">";
        return ss.str();
    }
private:
    friend Quantity operator+(const Quantity& q1, const Quantity& q2);
    friend Quantity operator-(const Quantity& q1, const Quantity& q2);
    friend Quantity operator*(const double& val, const Quantity& q);
    friend Quantity operator*(const Quantity& q1, const Quantity& q2);
    friend Quantity operator/(const Quantity& q1, const Quantity& q2);

    friend bool operator< (const Quantity& q1, const Quantity& q2);
    friend bool operator<=(const Quantity& q1, const Quantity& q2);
    friend bool operator==(const Quantity& q1, const Quantity& q2);
    friend bool operator>=(const Quantity& q1, const Quantity& q2);
    friend bool operator> (const Quantity& q1, const Quantity& q2);

    friend bool similar(const Quantity& q1, const Quantity& q2);

    unit::Unit un;
    double    val;
};

Quantity operator+(const Quantity& q1, const Quantity& q2) {
    Quantity q;
    if (similar(q1.un, q2.un))
        q = Quantity(q1(q1.un) + q2(q1.un), q1.un);
    else {
        std::cerr << "runtime error: cannot sum quantities " 
                  << q1.toString() << " and " << q2.toString()
                  << std::endl;
        exit(0);
    }
    return q;
}

Quantity operator-(const Quantity& q1, const Quantity& q2) {
    Quantity q;
    if (similar(q1.un, q2.un))
        q = Quantity(q1(q1.un) - q2(q1.un), q1.un);
    else {
        std::cerr << "runtime error: cannot sum quantities " 
                  << q1.toString() << " and " << q2.toString()
                  << std::endl;
        exit(0);
    }
    return q;
}

Quantity operator*(const double& val, const Quantity& q) {
	return Quantity(val*q.val, q.un);
}

Quantity operator*(const Quantity& q1, const Quantity& q2) {
    return Quantity(q1.val*q2.val, q1.un*q2.un);
}

Quantity operator/(const Quantity& q1, const Quantity& q2) {
    return Quantity(q1.val/q2.val, q1.un/q2.un);
}

bool operator< (const Quantity& q1, const Quantity& q2) {
    return q1(q1.un) <  q2(q1.un);
}
bool operator<=(const Quantity& q1, const Quantity& q2) {
    return q1(q1.un) <= q2(q1.un);
}
bool operator==(const Quantity& q1, const Quantity& q2) {
    return q1(q1.un) == q2(q1.un);
}
bool operator>=(const Quantity& q1, const Quantity& q2) {
    return q1(q1.un) >= q2(q1.un);
}
bool operator> (const Quantity& q1, const Quantity& q2) {
    return q1(q1.un) >  q2(q1.un);
}
bool similar(const Quantity& q1, const Quantity& q2) {
    return unit::similar(q1.un, q2.un);
}

}
}