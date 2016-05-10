#ifndef __DATA__
#define __DATA__

#include <Data/DATA_TYPE.h>
#include <Data/RANGE.h>
#include <Utilities/EIGEN_HELPERS.h>
#include <Utilities/TYPE_UTILITIES.h>
#include <memory>
#include <Eigen/Geometry>
#include <unordered_map>

namespace Mechanics{
template<class T> class RANDOM;

template<class TV>
class DATA:public std::vector<std::shared_ptr<DATA_TYPE<TV>>>
{
    typedef typename TV::Scalar T;

public:
    std::unordered_map<std::string,T> globals;
    RANDOM<T>& random;
    RANGE<TV> domain;
    bool periodic;

    DATA();
    ~DATA();

    TV Wrap(const TV& unwrapped) const{
        if(periodic){return domain.Wrap(unwrapped);} return unwrapped;
    }

    TV Minimum_Offset(const TV& X1,const TV& X2) const{
        if(periodic){return cwiseMinMag((X2-X1).eval(),cwiseMinMag((X2-(X1-domain.Edge_Lengths())).eval(),((X2-domain.Edge_Lengths())-X1).eval()));}
        return X2-X1;
    }

    template<class SUBTYPE> std::shared_ptr<SUBTYPE> Find_Or_Create() {
        Finder<std::shared_ptr<DATA_TYPE<TV>>> finder={SUBTYPE::Static_Name()};
        auto found=std::find_if(this->begin(),this->end(),finder);
        if(found==this->end()){
            auto data=std::make_shared<SUBTYPE>();
            this->push_back(data);
            return data;}
        return std::static_pointer_cast<SUBTYPE>(*found);
    }

    template<class SUBTYPE> std::shared_ptr<SUBTYPE> Find() const {
        Finder<std::shared_ptr<DATA_TYPE<TV>>> finder={SUBTYPE::Static_Name()};
        auto found=std::find_if(this->begin(),this->end(),finder);
        return std::static_pointer_cast<SUBTYPE>(*found);
    }

    template<class SUBTYPE> int Index() const {
        Finder<std::shared_ptr<DATA_TYPE<TV>>> finder={SUBTYPE::Static_Name()};
        auto found=std::find_if(this->begin(),this->end(),finder);
        return std::distance(this->begin(),found);
    }

    int Velocity_DOF() const;
    int Position_DOF() const;
    void Pack_Velocities(Matrix<T,Dynamic,1>& velocities) const;
    void Unpack_Velocities(const Matrix<T,Dynamic,1>& velocities);
    void Pack_Positions(Matrix<T,Dynamic,1>& positions) const;
    void Unpack_Positions(const Matrix<T,Dynamic,1>& positions);
    void Store_Errors(const Matrix<T,Dynamic,1>& errors);
    void Step();
    void Eliminate_Rotation();
    DEFINE_TYPE_NAME("DATA")
};
}
#endif
