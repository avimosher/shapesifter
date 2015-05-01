#ifndef __RANGE__
#define __RANGE__

namespace Mechanics{


template<class TV>
class RANGE
{
public:
    TV minimum_corner;
    TV maximum_corner;

    TV Wrap(const TV& X) const
    {
        TV wrapped;TV edges=Edge_Lengths();
        for(int i=0;i<TV::RowsAtCompileTime;i++){
            if(X[i]>=minimum_corner[i]){
                wrapped[i]=X[i]-((int)((X[i]-minimum_corner[i])/edges[i]))*edges[i];}
            else{
                wrapped[i]=X[i]+((int)((maximum_corner[i]-X[i])/edges[i]))*edges[i];}}
        return wrapped;
    }

    TV Edge_Lengths() const
    {return maximum_corner-minimum_corner;}
};
}

#endif
