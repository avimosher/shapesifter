#ifndef __PARSE_TEST_DATA__
#define __PARSE_TEST_DATA__

#include <Data/TEST_DATA.h>
#include <json/json.h>

namespace Mechanics{
template<class TV> class SIMULATION;

template<class TV>
class PARSE_TEST_DATA
{
public:
    PARSE_TEST_DATA(){};
    ~PARSE_TEST_DATA(){};

    static void Parse(Json::Value& node,SIMULATION<TV>& simulation);
    static std::string Static_Name(){
        return TEST_DATA<TV>::Static_Name();
    }
};
}
#endif
