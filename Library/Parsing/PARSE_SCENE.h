#ifndef __PARSE_SCENE__
#define __PARSE_SCENE__

#include <Utilities/TYPE_UTILITIES.h>
#include <memory>

namespace Mechanics{
template<class TV> class DATA;

template<class TV>
class PARSE_SCENE
{
    typedef typename TV::Scalar T;

public:
    PARSE_SCENE(){};
    ~PARSE_SCENE(){};

    static bool Parse_Scene(std::istream& input,DATA<TV>& data);
};
}
#endif
