#ifndef __SPRING_FORCE__
#define __SPRING_FORCE__

#include <Force/FORCE_TYPE.h>

namespace Mechanics{

template<class TV>
class SPRING_FORCE:public FORCE_TYPE<TV>
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
    enum DEFINITIONS{d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    typedef Matrix<T,t+d,1> FORCE_VECTOR;
    typedef Matrix<T,Dynamic,1> VECTOR;
public:
    using FORCE_TYPE<TV>::errors;
    struct SPRING{
        std::array<int,2> index;
        std::array<TV,2> offset;
        T target_distance;
        T stiffness;
    };
    std::vector<SPRING> springs;

    SPRING_FORCE(){}
    ~SPRING_FORCE(){}

    void Archive(cereal::BinaryOutputArchive& archive){}
    void Archive(cereal::BinaryInputArchive& archive){}

    void Archive(cereal::JSONOutputArchive& archive){}
    void Archive(cereal::JSONInputArchive& archive){}

    void Identify_Interactions_And_Compute_Errors(DATA<TV>& data,FORCE<TV>& force,const T dt,const T target_time,MATRIX_BUNDLE<TV>& system,bool stochastic);
    void Compute_Derivatives(DATA<TV>& data,FORCE<TV>& force,MATRIX_BUNDLE<TV>& system);
    DEFINE_TYPE_NAME("SPRING_FORCE")
};
}
#endif
