#ifndef __RIGID_STRUCTURE__
#define __RIGID_STRUCTURE__

#include <Data/FRAME.h>
#include <Data/MOMENT.h>
#include <Data/SUBSTRUCTURE.h>
#include <Data/TWIST.h>
#include <Utilities/CEREAL_HELPERS.h>
#include <Utilities/TYPE_UTILITIES.h>

namespace Mechanics{
template<class TV> class DATA;

template<class TV>
class RIGID_STRUCTURE
{
    typedef typename TV::Scalar T;
    typedef typename ROTATION<TV>::SPIN T_SPIN;
    bool initialized;
public:
    enum DEFINITIONS{PositionSize=FRAME<TV>::STATIC_SIZE,VelocitySize=TWIST<TV>::STATIC_SIZE};
    std::string name;
    FRAME<TV> frame;
    TWIST<TV> error;
    TWIST<TV> twist;
    MOMENT<TV> moi;
    bool kinematic;
    std::vector<SUBSTRUCTURE<TV>> substructures;

    RIGID_STRUCTURE():initialized(false),kinematic(false){}
    ~RIGID_STRUCTURE(){}

    DiagonalMatrix<T,VelocitySize> Inertia_Matrix(){
        assert(initialized);
        Matrix<T,VelocitySize,1> inertia;
        inertia<<moi.translation,moi.rotation;
        inertia*=-1;
        return inertia.asDiagonal();
    }

    AlignedBox<T,3> Bounding_Box(){
        TV minimum=TV::Constant(std::numeric_limits<T>::max());
        TV maximum=TV::Constant(std::numeric_limits<T>::lowest());
        for(auto& substructure : substructures){
            TV rotated_extent=(frame.orientation*substructure.capsule_extent).cwiseAbs();
            minimum=minimum.cwiseMin(frame.position-rotated_extent-TV::Constant(substructure.radius));
            maximum=maximum.cwiseMax(frame.position+rotated_extent+TV::Constant(substructure.radius));}
        return AlignedBox<T,3>(minimum,maximum);
    }

    SUBSTRUCTURE<TV>& Substructure(const int index){return substructures[index];}
    const SUBSTRUCTURE<TV>& Substructure(const int index) const {return substructures[index];}
    std::vector<SUBSTRUCTURE<TV>>& Substructures(){return substructures;}
    const std::vector<SUBSTRUCTURE<TV>>& Substructures() const {return substructures;}


    template<class Archive>
    void serialize(Archive& archive) {archive(CEREAL_NVP(name),
            CEREAL_NVP(frame),
            CEREAL_NVP(error),
            CEREAL_NVP(moi),
            CEREAL_NVP(twist),
            CEREAL_NVP(substructures),
            CEREAL_NVP(initialized));}

    void Initialize_Inertia(const T eta);
    DEFINE_TYPE_NAME("RIGID_STRUCTURE")
};
}
#endif
