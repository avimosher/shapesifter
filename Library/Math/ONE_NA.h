#ifndef __ONE_NA__
#define __ONE_NA__

#include <Math/Function.h>

namespace Mechanics{
// 1/|a|
template<class TV>
struct ONE_NA
{
    typedef typename TV::Scalar T;
    static T Evaluate(const TV& a,const T na){return 1/na;}
    static TV First_Derivative(const TV& a,const T na){return -a/cube(na);}
};
}
#endif
