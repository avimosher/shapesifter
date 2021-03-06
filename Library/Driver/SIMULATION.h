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
    int output_number;
    T dt;
    T time;
    T last_time;
    bool restart;
    bool substeps;
    bool write;
    std::string output_directory;
    std::string title;

    SIMULATION();
    ~SIMULATION();

    void Set_Restart(int restart_frame)
    {this->restart_frame=restart_frame;restart=true;}

    void Write(const std::string& frame_title);
    bool Read(const int frame);
};
}
#endif
