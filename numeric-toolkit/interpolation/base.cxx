#include <numeric-toolkit/interpolation/base.h>
#include <cmath>
#include <iostream>

using namespace numtk::interpolation;

Base::Base(
	const double*         x, 
	const double*         y, 
	std::size_t      xysize, 
	std::size_t _interpSize
) : xx(x, x + xysize), 
    yy(y, y + xysize), 
    interpSize(_interpSize),
    correlated(false),
    isaved(0) {
  	corDistance = std::max(1, (int) std::pow((double) xysize, 0.25));
}

Base::Base(
	const std::vector<double>&  x, 
	const std::vector<double>&  y, 
	std::size_t _interpSize
) : xx(x), yy(y), 
    interpSize(_interpSize),
    correlated(false),
    isaved(0) {
  	corDistance = std::max(1, (int) std::pow((double) x.size(), 0.25));
}

Base::Base(const Base& interpolator) {
	xx          = interpolator.xx;
	yy          = interpolator.yy;
	correlated  = interpolator.correlated;
	corDistance = interpolator.corDistance;
	interpSize  = interpolator.interpSize;
	isaved      = interpolator.isaved;
}

Base& Base::operator=(const Base& interpolator) {
	xx          = interpolator.xx;
	yy          = interpolator.yy;
	correlated  = interpolator.correlated;
	corDistance = interpolator.corDistance;
	interpSize  = interpolator.interpSize;
	isaved      = interpolator.isaved;
	return *this;
}

double Base::operator()(const double x) {
	return interpolate(correlated ? hunt(x) : locate(x), x);
}

double Base::interpolate(std::size_t location, double x) {
	return 0.0;
}

std::size_t Base::locate(const double x) {
	std::size_t iup, imid, ilow;
	std::size_t size = xx.size();
	if (size < 2 || interpSize < 2 || interpSize > size) 
		throw("locate size error");
	bool ascending = (xx[size - 1] >= xx[0]);
	ilow = 0;
	iup  = size - 1;
	while (iup - ilow > 1) {
		imid = (iup + ilow)/2;
		if (x >= xx[imid] == ascending)
			ilow = imid;
		else
			iup  = imid;
	}
	std::size_t distance = ilow > isaved ? ilow - isaved : isaved - ilow;
	correlated = distance < corDistance;
	isaved = ilow;
	std::size_t result1 = size - interpSize;
	std::size_t result2 = ilow > (interpSize - 2)/2 ? ilow - (interpSize - 2)/2 : 0; 
	return std::min(result1, result2);
}

std::size_t Base::hunt(const double x) {
	std::size_t iup, imid, ilow;
	std::size_t size = xx.size();
	std::size_t increment = 1;
	if (size < 2 || interpSize < 2 || interpSize > size) 
		throw("locate size error");
	bool ascending = (xx[size - 1] >= xx[0]);
	if (x >= xx[isaved] == ascending) { // hunt up
		ilow = isaved;
		for(;;) {
			iup = (ilow + increment > size - 1) ? (ilow + increment) : (size - 1);
			if (iup == size - 1) break;
			else if (x < xx[iup] == ascending) break;
			else { ilow = iup; increment += increment; }
		}
	}
	else { // hunt down
		iup = isaved;
		for(;;) {
			ilow = iup > increment ? iup - increment : 0;
			if (ilow == 0) break;
			else if (x >= xx[ilow] == ascending) break;
			else { iup = ilow; increment += increment; }
		}
	}
	// bisection phase
	while (iup - ilow > 1) {
		imid = (iup + ilow)/2;
		if (x >= xx[imid] == ascending)
			ilow = imid;
		else
			iup  = imid;
	}
	std::size_t distance = ilow > isaved ? ilow - isaved : isaved - ilow;
	correlated = distance < corDistance;
	isaved = ilow;
	std::size_t result1 = size - interpSize;
	std::size_t result2 = ilow > (interpSize - 2)/2 ? ilow - (interpSize - 2)/2 : 0; 
	return std::min(result1, result2);
}
