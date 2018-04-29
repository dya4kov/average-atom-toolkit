#pragma once
#include <vector>
#include <iostream>
#include <iomanip>

#include <average-atom-toolkit/util/unit.h>
#include <average-atom-toolkit/util/quantity.h>

namespace aatk {
namespace util {

class QuantityArray {
public:
    // initialize with specified size
    QuantityArray(const unit::Unit& u1, const unit::Unit& u2, const std::size_t& size) {
        if (!unit::similar(u1, u2)) {
            std::cerr << "runtime error: cannot create QuantityArray for units " 
                  << u1.toString() << " and " << u2.toString()
                  << std::endl;
            exit(0);
        }
        double val1 = u1.value(); 
        double val2 = u2.value();
        unit::Unit un1 = 1.0/val1*u1;
        unit::Unit un2 = 1.0/val2*u2; 
        Quantity q1(val1, un1);
        Quantity q2(val2, un2);
        val1 = q1(un1);
        val2 = q2(un1);
        un = un1;
        if (val1 >= val2) {
            std::cerr << "runtime error: cannot create QuantityArray for range " 
                  << val1 << ", " << val2
                  << std::endl;
            exit(0);
        }
        array.resize(size);
        if (size == 0) return;
        if (size == 1) { array[0] = val1; return; }
        double valStep = (val2 - val1)/(size - 1);
        for (std::size_t i = 0; i < size; ++i) {
            array[i] = val1 + i*valStep;
        }
    }
    // initialize with specified size
    QuantityArray(const Quantity& q1, const Quantity& q2, const std::size_t& size) {
        if (!similar(q1, q2)) {
            std::cerr << "runtime error: cannot create QuantityArray for quantities " 
                  << q1.toString() << " and " << q2.toString()
                  << std::endl;
            exit(0);
        }
        un = q1.unit();
        double val1 = q1(un);
        double val2 = q2(un);
        if (val1 >= val2) {
            std::cerr << "runtime error: cannot create QuantityArray for range " 
                  << val1 << ", " << val2
                  << std::endl;
            exit(0);
        }
        array.resize(size);
        if (size == 0) return;
        if (size == 1) array[0] = val1;
        double valStep = (val2 - val1)/(size - 1);
        for (std::size_t i = 0; i < size; ++i) {
            array[i] = val1 + i*valStep;
        }
    }
    // initialize from array
    QuantityArray(const std::vector<double>& values, const unit::Unit& u) {
        array = values;
        un = u;
    }
    QuantityArray(const QuantityArray& qa) {
        array = qa.array;
        un = qa.un;
    }
    QuantityArray& operator=(const QuantityArray& qa) {
        array = qa.array;
        un = qa.un;
        return *this;
    }
    Quantity operator[](const std::size_t& i) {
        return Quantity(array[i], un);
    }
    unit::Unit unit() const { return un; }

    QuantityArray& operator()(const unit::Unit& u) {
        if (similar(u, un)) {
            for (auto&& val : array)
                val *= un.value()/u.value();
            un = u;
        }
        else {
            std::cerr << "runtime error: cannot express unit " 
                  << un.toString() << " through " << u.toString()
                  << std::endl;
            exit(0);
        }
        return *this;
    }
private:
    friend std::ostream& operator<<(std::ostream& stream, const QuantityArray& qa);

    std::vector<double> array;
    unit::Unit          un;
};

std::ostream& operator<<(std::ostream& stream, const QuantityArray& qa) {
    stream << qa.un.toString() << std::endl;
    stream << std::scientific << "[";
    for (std::size_t i = 0; i < qa.array.size() - 1; ++i)
        stream << qa.array[i] << " ";
    stream << qa.array[qa.array.size() - 1] << "]";
    return stream;
}

}
}