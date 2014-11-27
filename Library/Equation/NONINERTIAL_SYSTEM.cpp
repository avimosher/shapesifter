///////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson/Mosher.
///////////////////////////////////////////////////////////////////////
// Class NONINERTIAL_SYSTEM
///////////////////////////////////////////////////////////////////////
#include <Evolution/NONINERTIAL_SYSTEM.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> NONINERTIAL_SYSTEM<TV>::
NONINERTIAL_SYSTEM()
{
}
///////////////////////////////////////////////////////////////////////
template<class TV> typename TV::Scalar NONINERTIAL_SYSTEM<TV>::
Solve(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time)
{
    const int newton_iterations=40;
    const int solve_iterations=200;
    

    for(int i=0;i<newton_iterations;i++){
        Compute(data,force,dt,time);

        BiCGSTAB<SparseMatrix<T> > solver;
        solver.compute(A);

        x=solver.solve(b);

    }
}
GENERIC_TYPE_DEFINITION(NONINERTIAL_SYSTEM)
