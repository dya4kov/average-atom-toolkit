#include <vector>

#include <numeric-tools/ODE/solver.h>
#include <numeric-tools/ODE/stepper/PD853.h>

#include <average-atom-tools/thomas-fermi/atom/table.h>
#include <average-atom-tools/thomas-fermi/atom/ODEs/potential.h>
#include <average-atom-tools/thomas-fermi/atom/ODEs/action.h>
#include <average-atom-tools/thomas-fermi/atom/ODEs/tau.h>
#include <average-atom-tools/thomas-fermi/atom/ODEs/rotatePoints.h>

namespace AATools {
namespace TF {

using ::numtools::ODE::Solver;
using ::numtools::ODE::stepper::PD853;

using ::numtools::specfunc::FermiDirac;
using ::numtools::specfunc::FD::Half;

class Atom : ODE::Potential, ODE::Action, ODE::Tau, ODE::RotatePoints {
public:
    Atom();

    double xU    (const double& x);
    double xUDX  (const double& x);
    double eDens (const double& x);

    std::vector<double>& xU    (const std::vector<double>& x);
    std::vector<double>& xUDX  (const std::vector<double>& x);
    std::vector<double>& eDens (const std::vector<double>& x);

    double* xU    (double* x, size_t n);
    double* xUDX  (double* x, size_t n);
    double* eDens (double* x, size_t n);

    double e     (const int& n, const int& l);
    double rpIn  (const int& n, const int& l);
    double rpOut (const int& n, const int& l);
    double tau   (const int& n, const int& l);

    void setV(const double& _V);
    void setT(const double& _T);
    void setZ(const double& _Z);
    void setVTZ(const double& _V, const double& _T, const double& _Z);

    void setTolerance (const double& _eps);

private:

    void adjust(); bool needAdjust;
    void updateMu();
    int idxLevel(const int& n);

    double action (const double& e, const double& l);
    double rPoint1(const double& e, const double& l);
    double rPoint2(const double& e, const double& l);
    double tauEval(const double& e, const double& l);

    double V1, T1;
    double VZ, TZ;
    double r0();
    double eps;
    double mu;
    
    PotTable table;
    FermiDirac<Half> FDhalf;

    std::vector<double> rPointI;
    std::vector<bool>   rpIready;
    std::vector<double> rPointO;
    std::vector<bool>   rpOready;

    std::vector<bool>   eReady;
    std::vector<double> eLevel;
    std::vector<double> eStart;

    std::vector<bool>   tReady;
    std::vector<double> tLevel;
    
    RHSPotential rhsPot; Solver<PD853<RHSPotential>> potSolver;
    RHSAction    rhsAct; Solver<PD853<RHSAction>>    actSolver;
    RHSTau       rhsTau; Solver<PD853<RHSTau>>       tauSolver;
    RHSRP1       rhsRP1; Solver<PD853<RHSRP1>>       RP1Solver;
    RHSRP2       rhsRP2; Solver<PD853<RHSRP2>>       RP2Solver;

};

}
}