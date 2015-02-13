#ifndef __MATH__
#define __MATH__

template<class T> inline T sqr(T x){return x*x;}
template<class T> inline T cube(T x){return x*x*x;}
template<class T>
inline T sinc(const T x) // sin(x)/x
{if(abs(x)<1e-8) return 1;return sin(x)/x;}



#endif
