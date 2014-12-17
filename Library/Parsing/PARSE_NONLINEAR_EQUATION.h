#ifndef __PARSE_NONLINEAR_EQUATION__
#define __PARSE_NONLINEAR_EQUATION__

#include <Equation/NONLINEAR_EQUATION.h>
#include <json/json.h>

namespace Mechanics{
template<class TV> class SIMULATION;

template<class TV>
class PARSE_NONLINEAR_EQUATION
{
public:
    PARSE_NONLINEAR_EQUATION(){};
    ~PARSE_NONLINEAR_EQUATION(){};

    static void Parse(Json::Value& node,SIMULATION<TV>& simulation);
    static std::string Static_Name(){
        return NONLINEAR_EQUATION<TV>::Static_Name();
    }
};
}
#endif
