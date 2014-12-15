//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class EVOLUTION_STEP
///////////////////////////////////////////////////////////////////////
#include <Evolution/EVOLUTION_STEP.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> EVOLUTION_STEP<TV>::
EVOLUTION_STEP()
{}
///////////////////////////////////////////////////////////////////////
template<class TV> EVOLUTION_STEP<TV>::
~EVOLUTION_STEP()
{}
///////////////////////////////////////////////////////////////////////
template<class TV> bool EVOLUTION_STEP<TV>::
Satisfied(DATA<TV>& data,FORCE<TV>& force,const T dt,const T time)
{
    return false;
}
///////////////////////////////////////////////////////////////////////
GENERIC_TYPE_DEFINITION(EVOLUTION_STEP)