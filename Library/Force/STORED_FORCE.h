#ifndef __STORED_FORCES__
#define __STORED_FORCES__

#include <Utilities/TYPE_UTILITIES.h>

namespace Mechanics{

template<class T>
struct FORCE_REFERENCE
{
    Matrix<T,Dynamic,1> value;

    FORCE_REFERENCE(){}
    virtual ~FORCE_REFERENCE(){};
    virtual int Size() const{return value.size();}
    virtual const Matrix<T,Dynamic,1>& Vector() const{return value;}
    DEFINE_TYPE_NAME("FORCE_REFERENCE");
};

template<class T>
struct STORED_FORCE:public std::vector<std::shared_ptr<FORCE_REFERENCE<T>>>
{
    template<class SUBTYPE> std::shared_ptr<SUBTYPE> Find_Or_Create() {
        Finder<std::shared_ptr<FORCE_REFERENCE<T>>> finder={SUBTYPE::Static_Name()};
        auto found=std::find_if(this->begin(),this->end(),finder);
        if(found==this->end()){
            auto reference=std::make_shared<SUBTYPE>();
            this->push_back(reference);
            return reference;
        }
        return std::static_pointer_cast<SUBTYPE>(*found);
    }

    template<class SUBTYPE> std::shared_ptr<SUBTYPE> Find() const {
        Finder<std::shared_ptr<FORCE_REFERENCE<T>>> finder={SUBTYPE::Static_Name()};
        return std::static_pointer_cast<SUBTYPE>(*std::find_if(this->begin(),this->end(),finder));
    }

    int Size() const{return this->size();}

    void Resize(int size){
        if(Size()!=size){this->resize(size);}
    }

    void setZero(){
        for(auto& stored_force:(*this)){stored_force->value.setZero();}
    }

    Matrix<T,Dynamic,1> Vector() const
    {
        int total_size=0;
        for(int i=0;i<Size();i++){total_size+=(*this)[i]->Size();}
        Matrix<T,Dynamic,1> vector(total_size);
        int current_position=0;
        for(int i=0;i<Size();i++){
            int size=(*this)[i]->Size();
            vector.block(current_position,0,size,1)=(*this)[i]->Vector();
            current_position+=size;}
        return vector;
    }

    void Set(const Matrix<T,Dynamic,1>& input_values){
        int current_position=0;
        for(int i=0;i<Size();i++){
            int force_size=(*this)[i]->value.size();
            (*this)[i]->value=input_values.block(current_position,0,force_size,1);
            current_position+=force_size;
        }
    }
};
}
#endif
