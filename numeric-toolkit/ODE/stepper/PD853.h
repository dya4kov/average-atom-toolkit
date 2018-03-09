#pragma once

#include <numeric-toolkit/ODE/types.h>

namespace numtk {
namespace ODE {
namespace stepper {

using ::numtk::ODE::Array;
using ::numtk::ODE::Dimension;

struct PD853_constants {
	const double 
	c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c14, c15, c16,
	b1, b6, b7, b8, b9, b10, b11, b12, bhh1, bhh2, bhh3,
	er1, er6, er7, er8, er9, er10, er11, er12,
	a21, a31, a32, a41, a43, a51, a53, a54, a61, a64, a65, a71, a74, a75, a76,
	a81, a84, a85, a86, a87, a91, a94, a95, a96, a97, a98, a101, a104, a105,	
	a106, a107, a108, a109, a111, a114, a115, a116, a117, a118, a119, a1110,
	a121, a124, a125, a126, a127, a128, a129, a1210, a1211, a141, a147, a148,
	a149, a1410, a1411, a1412, a1413, a151, a156, a157, a158, a1511, a1512,
	a1513, a1514, a161, a166, a167, a168, a169, a1613, a1614, a1615;
	PD853_constants() :
	   c2( 0.526001519587677318785587544488e-01), 
	   c3( 0.789002279381515978178381316732e-01), 
	   c4( 0.118350341907227396726757197510e+00), 
	   c5( 0.281649658092772603273242802490e+00), 
	   c6( 0.333333333333333333333333333333e+00), 
	   c7( 0.250000000000000000000000000000e+00), 
	   c8( 0.307692307692307692307692307692e+00), 
	   c9( 0.651282051282051282051282051282e+00), 
	  c10( 0.600000000000000000000000000000e+00), 
	  c11( 0.857142857142857142857142857142e+00), 
	  c14( 0.100000000000000000000000000000e+00), 
	  c15( 0.200000000000000000000000000000e+00), 
	  c16( 0.777777777777777777777777777778e+00),
	   b1( 5.42937341165687622380535766363e-02), 
	   b6( 4.45031289275240888144113950566e+00), 
	   b7( 1.89151789931450038304281599044e+00), 
	   b8(-5.80120396001058478146721142270e+00), 
	   b9( 3.11164366957819894408916062370e-01), 
	  b10(-1.52160949662516078556178806805e-01), 
	  b11( 2.01365400804030348374776537501e-01), 
	  b12( 4.47106157277725905176885569043e-02), 
	 bhh1(0.244094488188976377952755905512e+00), 
	 bhh2(0.733846688281611857341361741547e+00), 
	 bhh3(0.220588235294117647058823529412e-01),
	  er1( 0.13120044994194880732501029960e-01), 
	  er6(-0.12251564463762044407205697530e+01), 
	  er7(-0.49575894965725019152140799520e+00), 
	  er8( 0.16643771824549865369615304150e+01), 
	  er9(-0.35032884874997368168864872900e+00), 
	 er10( 0.33417911871301747902973188410e+00), 
	 er11( 0.81923206485115712465707426130e-01), 
	 er12(-0.22355307863886295258844278450e-01),
	  a21( 5.26001519587677318785587544488e-02), 
	  a31( 1.97250569845378994544595329183e-02), 
	  a32( 5.91751709536136983633785987549e-02), 
	  a41( 2.95875854768068491816892993775e-02), 
	  a43( 8.87627564304205475450678981324e-02), 
	  a51( 2.41365134159266685502369798665e-01),
	  a53(-8.84549479328286085344864962717e-01), 
	  a54( 9.24834003261792003115737966543e-01), 
	  a61( 3.70370370370370370370370370370e-02), 
	  a64( 1.70828608729473871279604482173e-01), 
	  a65( 1.25467687566822425016691814123e-01), 
	  a71( 3.71093750000000000000000000000e-02), 
	  a74( 1.70252211019544039314978060272e-01), 
	  a75( 6.02165389804559606850219397283e-02), 
	  a76(-1.75781250000000000000000000000e-02),
	  a81( 3.70920001185047927108779319836e-02), 
	  a84( 1.70383925712239993810214054705e-01), 
	  a85( 1.07262030446373284651809199168e-01), 
	  a86(-1.53194377486244017527936158236e-02), 
	  a87( 8.27378916381402288758473766002e-03), 
	  a91( 6.24110958716075717114429577812e-01), 
	  a94(-3.36089262944694129406857109825e+00), 
	  a95(-8.68219346841726006818189891453e-01), 
	  a96( 2.75920996994467083049415600797e+01), 
	  a97( 2.01540675504778934086186788979e+01), 
	  a98(-4.34898841810699588477366255144e+01), 
	 a101( 4.77662536438264365890433908527e-01), 
	 a104(-2.48811461997166764192642586468e+00), 
	 a105(-5.90290826836842996371446475743e-01),	
	 a106( 2.12300514481811942347288949897e+01), 
	 a107( 1.52792336328824235832596922938e+01), 
	 a108(-3.32882109689848629194453265587e+01), 
	 a109(-2.03312017085086261358222928593e-02), 
	 a111(-9.37142430085987325717040216580e-01), 
	 a114( 5.18637242884406370830023853209e+00), 
	 a115( 1.09143734899672957818500254654e+00), 
	 a116(-8.14978701074692612513997267357e+00), 
	 a117(-1.85200656599969598641566180701e+01), 
	 a118( 2.27394870993505042818970056734e+01), 
	 a119( 2.49360555267965238987089396762e+00), 
	a1110(-3.04676447189821950038236690220e+00),
	 a121( 2.27331014751653820792359768449e+00), 
	 a124(-1.05344954667372501984066689879e+01), 
	 a125(-2.00087205822486249909675718444e+00), 
	 a126(-1.79589318631187989172765950534e+01), 
	 a127( 2.79488845294199600508499808837e+01), 
	 a128(-2.85899827713502369474065508674e+00),
	 a129(-8.87285693353062954433549289258e+00), 
	a1210( 1.23605671757943030647266201528e+01), 
	a1211( 6.43392746015763530355970484046e-01), 
	 a141( 5.61675022830479523392909219681e-02), 
	 a147( 2.53500210216624811088794765333e-01), 
	 a148(-2.46239037470802489917441475441e-01),
	 a149(-1.24191423263816360469010140626e-01), 
	a1410( 1.53291798278765697312063226850e-01), 
	a1411( 8.20105229563468988491666602057e-03), 
	a1412( 7.56789766054569976138603589584e-03), 
	a1413(-8.29800000000000000000000000000e-03), 
	 a151( 3.18346481635021405060768473261e-02), 
	 a156( 2.83009096723667755288322961402e-02), 
	 a157( 5.35419883074385676223797384372e-02), 
	 a158(-5.49237485713909884646569340306e-02), 
	a1511(-1.08347328697249322858509316994e-04), 
	a1512( 3.82571090835658412954920192323e-04),
	a1513(-3.40465008687404560802977114492e-04), 
	a1514( 1.41312443674632500278074618366e-01), 
	 a161(-4.28896301583791923408573538692e-01), 
	 a166(-4.69762141536116384314449447206e+00), 
	 a167( 7.68342119606259904184240953878e+00), 
	 a168( 4.06898981839711007970213554331e+00), 
	 a169( 3.56727187455281109270669543021e-01), 
	a1613(-1.39902416515901462129418009734e-03), 
	a1614( 2.94751478915277233895562721490e+00), 
	a1615(-9.15095847217987001081870187138e+00) {}
};

template<class RHStype>
class PD853 : public PD853_constants {
	// Dormand-Prince fifth-order Runge-Kutta step with monitoring of 
	// local truncation error to ensure accuracy and adjust stepsize.
public:
	typedef RHStype RHS;

	PD853(
		Array<RHStype::dim> &_y, 
		Array<RHStype::dim> &_dydx, 
		double               &_x, 
		double               &_tolAbs, 
		double               &_tolRel, 
		double               &_EPS) : 
    	y(_y), dydx(_dydx), x(_x), 
    	tolAbs(_tolAbs), tolRel(_tolRel), EPS(_EPS),
    	PD853_constants() {}

	void dy(const double h, RHStype &rhs);

	void   makeStep(const double hTry, RHStype &rhs);
	double lastStep() { return hDid;  }
	double nextStep() { return hNext; }

	double error(const double h);

	struct Controller {
		Controller();
		bool success(const double err, double &h);
		double hNext,errOld;
		bool reject;
	};
private:
	double              &x;
	Array<RHStype::dim> &y, &dydx;
	Array<RHStype::dim> k2,k3,k4,k5,k6,k7,k8,k9,k10;
	Array<RHStype::dim> dydxnew;
	Array<RHStype::dim> yOut, yErr;
	Array<RHStype::dim> yErr2;
	Controller          controller;

	double &tolAbs, &tolRel, &EPS;
	double hDid, hNext;
};

template<class RHStype>
void PD853<RHStype>::makeStep(const double hTry, RHStype &rhs) {
	double h = hTry; // Set stepsize to the initial trial value
	for (;;) {
		dy(h, rhs); // Take a step
		double err = error(h); // Evaluate accuracy
		if (controller.success(err, h)) break; 
		// Else step rejected. Try again with
        // reduced h set by controller
		if (std::abs(h) <= EPS) break;
	}
	rhs(x + h, yOut, dydxnew);
	dydx = dydxnew; // Reuse last derivative evaluation for next step
	y = yOut;
	hDid = h;
	x += hDid;
	hNext = controller.hNext;
}

template<class RHStype>
void PD853<RHStype>::dy(const double h, RHStype &rhs) {
	Array<RHStype::dim> yTemp; yTemp.fill(0.0);
	Dimension i;
	for (i = 0; i < RHStype::dim; ++i) {
		yTemp[i] = y[i] + h*a21*dydx[i];
	}
	rhs(x + c2*h, yTemp, k2);
	for (i = 0; i < RHStype::dim; ++i) {
		yTemp[i] = y[i] + h*(a31*dydx[i] + a32*k2[i]);
	}
	rhs(x + c3*h, yTemp, k3);
	for (i = 0; i < RHStype::dim; ++i) {
		yTemp[i] = y[i] + h*(a41*dydx[i] + a43*k3[i]);
	}
	rhs(x + c4*h, yTemp, k4);
	for (i = 0; i < RHStype::dim; ++i) {
		yTemp[i] = y[i] + h*(a51*dydx[i] + a53*k3[i] + a54*k4[i]);
	}
	rhs(x + c5*h, yTemp, k5);
	for (i = 0; i < RHStype::dim; ++i) {
		yTemp[i] = y[i] + h*(a61*dydx[i] + a64*k4[i] + a65*k5[i]);
	}
	rhs(x + c6*h, yTemp, k6);
	for (i = 0; i < RHStype::dim; ++i) {
		yTemp[i] = y[i] + h*(a71*dydx[i] + a74*k4[i] + a75*k5[i] + a76*k6[i]);
	}
	rhs(x + c7*h, yTemp, k7);
	for (i = 0; i < RHStype::dim; ++i) {
		yTemp[i] = y[i] + h*(a81*dydx[i] + a84*k4[i] + a85*k5[i] + a86*k6[i] +
							 a87*k7[i]);
	}
	rhs(x + c8*h, yTemp, k8);
	for (i = 0; i < RHStype::dim; ++i) {
		yTemp[i] = y[i] + h*(a91*dydx[i] + a94*k4[i] + a95*k5[i] + a96*k6[i] + 
							 a97*k7[i] + a98*k8[i]);
	}
	rhs(x + c9*h, yTemp, k9);
	for (i = 0; i < RHStype::dim; ++i) {
		yTemp[i] = y[i] + h*(a101*dydx[i] + a104*k4[i] + a105*k5[i] + a106*k6[i] +
							 a107*k7[i] + a108*k8[i] + a109*k9[i]);
	}
	rhs(x + c10*h, yTemp, k10);
	for (i = 0; i < RHStype::dim; ++i) {
		yTemp[i] = y[i] + h*(a111*dydx[i] + a114*k4[i] + a115*k5[i] + a116*k6[i] +
							 a117*k7[i] + a118*k8[i] + a119*k9[i] + a1110*k10[i]);
	}
	rhs(x + c11*h, yTemp, k2);
	double xph = x + h;
	for (i = 0; i < RHStype::dim; ++i) {
		yTemp[i] = y[i] + h*(a121*dydx[i] + a124*k4[i] + a125*k5[i] + a126*k6[i] +
							 a127*k7[i] + a128*k8[i] + a129*k9[i] + a1210*k10[i] + 
							 a1211*k2[i]);
	}
	rhs(xph, yTemp, k3);
	for (i = 0; i < RHStype::dim; ++i) {
		k4[i] = b1*dydx[i] + b6*k6[i]  + b7*k7[i] + b8*k8[i] + b9*k9[i] +
				b10*k10[i] + b11*k2[i] + b12*k3[i];
		yOut[i] = y[i] + h*k4[i];
	}
	for (i = 0; i < RHStype::dim; ++i) {
		yErr[i] = k4[i] - bhh1*dydx[i] - bhh2*k9[i] - bhh3*k3[i];
		yErr2[i] = er1*dydx[i] + er6*k6[i]  + er7*k7[i] + er8*k8[i] + er9*k9[i] +
				   er10*k10[i] + er11*k2[i] + er12*k3[i];
	}
}

template<class RHStype>
double PD853<RHStype>::error(const double h) {
	double err = 0.0;
	double err2 = 0.0;
	double sk, deno, erri;
	for (Dimension i = 0; i < RHStype::dim; ++i) {
		sk = tolAbs + tolRel*std::max(std::abs(y[i]), std::abs(yOut[i]));
		erri = yErr[i]/sk; erri *= erri;
		err2 += erri;
		erri = yErr2[i]/sk; erri *= erri;
		err += erri;
	}
	deno = err + 0.01*err2;
	if (deno <= 0.0) {
		deno = 1.0;
	}
	return std::abs(h)*err*std::sqrt(1.0/(RHStype::dim*deno));
}

template<class RHStype>
PD853<RHStype>::Controller::Controller() : reject(false), errOld(1.0e-4) {}

template <class RHStype>
bool PD853<RHStype>::Controller::success(const double err, double &h) {
	// Returns true if err <= 1, false otherwise. If step was successful, 
	// sets hNext to the estimated optimal stepsize for the next step.
	// If the step failed, reduces h appropriately for another try.
	const double 
	beta = 0.04, 
	alpha = 1.0/8.0 - beta*0.2, 
	safe = 0.9, 
	minscale = 0.333,
	maxscale = 6.0;
	// Set beta to a nonzero value for PI control. beta = 0.04-0.08 is a good default.
	double scale;
	if (err <= 1.0) { // Step succeeded. Compute hNext.
		if (err == 0.0) {
			scale = maxscale;
		}
		else { // PI control if beta != 0.
			scale = safe*std::pow(err, -alpha)*std::pow(errOld, beta);
			// Ensure minscale <= hNext/h <= maxscale.
			if (scale < minscale) scale = minscale;  
			if (scale > maxscale) scale = maxscale;
		}
		// Don’t let step increase if last one was rejected
		if (reject) {
			hNext = h*std::min(scale,1.0);
		}
		else {
			hNext = h*scale;
		}
		errOld = std::max(err, 1.0e-4); // Bookkeeping for next call
		reject = false;
		return true;
	}
	else {
		scale = std::max(safe*std::pow(err, -alpha), minscale);
		h *= scale;
		reject = true;
		return false;
	}
}

} // end of namespace stepper
} // end of namespace ODE
} // end of namespace numtk