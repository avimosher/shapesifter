#ifndef __RANGE__
#define __RANGE__

namespace Mechanics{


template<class TV>
class RANGE:public Matrix<TV,2,1>
{
public:
    TV Wrap(const TV& X) const
    {
        TV wrapped;TV edges=Edge_Lengths();
        for(int i=0;i<TV::RowsAtCompileTime;i++){
            if(X[i]>=(*this)[0][i]){
                wrapped[i]=X[i]-((int)((X[i]-(*this)[0][i])/edges[i]))*edges[i];}
            else{
                wrapped[i]=X[i]+((int)(((*this)[1][i]-X[i])/edges[i]))*edges[i];}}
        return wrapped;
    }

    TV Edge_Lengths() const
    {return (*this)[1]-(*this)[0];}
};
}

#endif
