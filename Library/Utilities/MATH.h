#ifndef __MATH__
#define __MATH__

#include <cmath>

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
template<class T> inline T sqr(T x){return x*x;}
template<class T> inline T cube(T x){return x*x*x;}

template<class T> inline T sinc(const T x) // sin(x)/x
{if(std::abs(x)<1e-8) return 1;return sin(x)/x;}

template<class T> inline T finite(const T x)
{return (std::abs(x)<=__DBL_MAX__) && (x==x);}

template<class T> inline T minmag(T a,T b)
{if(std::abs(a)<std::abs(b)){return a;} return b;}

template<class T> T clamp(T a,T minimum,T maximum)
{if(a<minimum){return minimum;}if(a>maximum){return maximum;}return a;}

#endif
