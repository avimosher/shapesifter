#ifndef __FORCE__
#define __FORCE__

#include <Utilities/TYPE_UTILITIES.h>
#include <memory>
#include <vector>
#include <Eigen/SparseCore>
#include <osg/Group>


namespace Mechanics{
template<class TV> class DATA;
template<class TV> class FORCE_TYPE;


template<class TV>
class FORCE:public std::vector<std::shared_ptr<FORCE_TYPE<TV>>>
{
    typedef typename TV::Scalar T;

public:
    FORCE(){}
    ~FORCE(){}

    template<class SUBTYPE> std::shared_ptr<SUBTYPE> Find_Or_Create() {
        Finder<std::shared_ptr<FORCE_TYPE<TV>>> finder={SUBTYPE::Static_Name()};
        auto found=std::find_if(this->begin(),this->end(),finder);
        if(found==this->end()){
            auto force=std::make_shared<SUBTYPE>();
            this->push_back(force);
            return force;
        }
        return std::static_pointer_cast<SUBTYPE>(*found);
    }

    T Compute_Dt(DATA<TV>& data,FORCE<TV>& force,const T target_time);
    int Force_DOF() const;
    void Pack_Forces(Matrix<T,Dynamic,1>& forces);
    void Unpack_Forces(const Matrix<T,Dynamic,1>& forces);
    void Viewer(const DATA<TV>& data,osg::Group*& root);
};
}
#endif
