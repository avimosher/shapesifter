//#####################################################################
// Copyright 2014, Avi Robinson-Mosher.
//#####################################################################
// Class DRIVER
//#####################################################################
#ifndef __DRIVER__
#define __DRIVER__

namespace Mechanics{
template<class TV> class DATA;

template<class TV>
class DRIVER
{
    typedef typename TV::SCALAR T;

public:
    DATA<TV>& data;
    int output_number;

    DRIVER(DATA<TV>& data);
    ~DRIVER();

    void Initialize() PHYSBAM_OVERRIDE;
    void Advance_One_Time_Step(const T target_time,bool &done);
    void Advance_To_Target_Time(const T target_time) PHYSBAM_OVERRIDE;
    void Execute_Main_Program() PHYSBAM_OVERRIDE;
};
}
#endif
