#pragma once

namespace aatk {

class EoS {
public:

    EoS();
    EoS(double Z, double A);

    setZ(double Z);
    setAmass(double A);

    EoS& operator()(std::string model);
    EoS& operator[](std::string params);

private:

    double Z;
    Quantity A;

};

}