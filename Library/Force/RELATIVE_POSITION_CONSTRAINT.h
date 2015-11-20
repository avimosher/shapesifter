#ifndef __RELATIVE_POSITION_CONSTRAINT__
#define __RELATIVE_POSITION_CONSTRAINT__

#include <Force/FORCE_TYPE.h>

namespace Mechanics{

template<class TV>
class RELATIVE_POSITION_CONSTRAINT:public FORCE_TYPE<TV>
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
    enum DEFINITIONS{d=TV::RowsAtCompileTime,t=T_SPIN::RowsAtCompileTime};
    typedef Matrix<T,1,t+d> CONSTRAINT_VECTOR;
    typedef Matrix<T,t+d,1> FORCE_VECTOR;
    typedef Matrix<T,Dynamic,1> VECTOR;
public:
    using FORCE_TYPE<TV>::stored_forces;using FORCE_TYPE<TV>::errors;
    struct CONSTRAINT{
        int s1;
        TV v1;
        int s2;
        TV v2;
        T target_distance;

        template<class Archive>
        void serialize(Archive& archive)
        {archive(s1,v1,s2,v2,target_distance);}
    };
    std::vector<CONSTRAINT> constraints;

    RELATIVE_POSITION_CONSTRAINT(){}
    ~RELATIVE_POSITION_CONSTRAINT(){}

    void Archive(cereal::BinaryOutputArchive& archive){archive(CEREAL_NVP(errors));}
    void Archive(cereal::BinaryInputArchive& archive){archive(CEREAL_NVP(errors));}

    void Archive(cereal::JSONOutputArchive& archive){archive(CEREAL_NVP(errors));}
    void Archive(cereal::JSONInputArchive& archive){archive(CEREAL_NVP(errors));}

    void Linearize(DATA<TV>& data,FORCE<TV>& force,const T dt,const T target_time,MATRIX_BUNDLE<TV>& system,bool stochastic);
    void Viewer(const DATA<TV>& data,osg::Node* node);
    DEFINE_TYPE_NAME("RELATIVE_POSITION_CONSTRAINT")
};
}
#endif
