//#####################################################################
// Copyright 2014, Avi Robinson-Mosher.
//#####################################################################
// Class DRIVER
//#####################################################################
#ifndef __DRIVER__
#define __DRIVER__

namespace Mechanics{
template<class TV> class DATA;
template<class TV> class EVOLUTION;
template<class TV> class FORCE;

template<class TV>
class DRIVER
{
    typedef typename TV::Scalar T;

public:
    DATA<TV>& data;
    EVOLUTION<TV>& evolution;
    FORCE<TV>& force;
    int output_number;
    int current_frame;
    const T time;

    DRIVER(DATA<TV>& data,EVOLUTION<TV>& evolution,FORCE<TV>& force);
    ~DRIVER();

    void Initialize() PHYSBAM_OVERRIDE;
    void Advance_One_Time_Step(const T target_time,bool &done);
    void Advance_To_Target_Time(const T target_time) PHYSBAM_OVERRIDE;
    void Execute_Main_Program() PHYSBAM_OVERRIDE;
};
}
#endif
