#ifndef __PARSE_TEST_DATA__
#define __PARSE_TEST_DATA__

#include <json/json.h>

namespace Mechanics{
template<class TV> class DATA;

template<class TV>
class PARSE_TEST_DATA
{
public:
    PARSE_TEST_DATA(){};
    ~PARSE_TEST_DATA(){};

    static void Parse(Json::Value& node,DATA<TV>& data);
    static std::string Static_Name(){
        return "TEST_DATA";
    }
};
}
#endif
