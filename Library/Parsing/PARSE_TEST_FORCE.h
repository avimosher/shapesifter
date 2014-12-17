#ifndef __PARSE_TEST_FORCE__
#define __PARSE_TEST_FORCE__

#include <Force/TEST_FORCE.h>
#include <json/json.h>

namespace Mechanics{
template<class TV> class SIMULATION;

template<class TV>
class PARSE_TEST_FORCE
{
public:
    PARSE_TEST_FORCE(){};
    ~PARSE_TEST_FORCE(){};

    static void Parse(Json::Value& node,SIMULATION<TV>& simulation);
    static std::string Static_Name(){
        return TEST_FORCE<TV>::Static_Name();
    }
};
}
#endif
