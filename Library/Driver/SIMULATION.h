//#####################################################################
// Copyright 2014, Avi Robinson-Mosher.
//#####################################################################
// Class SIMULATION
//#####################################################################
#ifndef __SIMULATION__
#define __SIMULATION__

#include <Utilities/TYPE_UTILITIES.h>
#include <memory>

namespace Mechanics{
template<class TV> class DATA;
template<class TV> class EVOLUTION;
template<class TV> class FORCE;

template<class TV>
class SIMULATION
{
    typedef typename TV::Scalar T;
public:
    DATA<TV>& data;
    EVOLUTION<TV>& evolution;
    FORCE<TV>& force;
    int current_frame;
    int restart_frame;
    T time;
    T last_time;
    bool restart;
    std::string output_directory;

    SIMULATION();
    ~SIMULATION();

    void Write(const int frame);
    void Read(const int frame);

};
}
#endif