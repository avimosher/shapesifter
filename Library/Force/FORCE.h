#ifndef __FORCE__
#define __FORCE__

#include <Force/FORCE_TYPE.h>
#include <Utilities/CEREAL_HELPERS.h>
#include <Utilities/TYPE_UTILITIES.h>
#include <memory>
#include <vector>
#include <Eigen/SparseCore>
#include <osg/Group>


namespace Mechanics{
template<class TV> class DATA;
template<class TV> class FORCE_TYPE;
template<class T> class STORED_FORCE;

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
            return force;}
        return std::static_pointer_cast<SUBTYPE>(*found);
    }

    template<class Archive>
    void serialize(Archive& archive)
    {for(int i=0;i<(*this).size();i++){(*this)[i]->Archive(archive);}}

    T Compute_Dt(DATA<TV>& data,FORCE<TV>& force,const T target_time);
    void Pack_Forces(STORED_FORCE<T>& stored_force) const;
    void Unpack_Forces(const STORED_FORCE<T>& stored_force);
    void Increment_Forces(const STORED_FORCE<T>& stored_force,T ratio);
    void Viewer(const DATA<TV>& data,osg::Group*& root);
    bool Equations_Changed() const;
};

}

namespace cereal
{
template<class Archive,class TV>
struct specialize<Archive,Mechanics::FORCE<TV>,cereal::specialization::member_serialize>{};
}
#endif
