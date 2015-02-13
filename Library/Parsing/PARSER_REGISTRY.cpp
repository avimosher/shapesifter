#include <Analysis/AGGREGATOR.h>
#include <Analysis/PREDICATE.h>
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/TYPE_UTILITIES.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV,class PARSED_TYPE> std::map<std::string,std::shared_ptr<PARSED_TYPE> (*)(Json::Value&,SIMULATION<TV>&)>& PARSER_REGISTRY<TV,PARSED_TYPE>::
Parsers()
{
    static std::map<std::string,Parse_Function> parser_registry;
    return parser_registry;
}
///////////////////////////////////////////////////////////////////////
#define SCALAR_TYPE_DEFINITION(TYPE,T,RET)           \
    template class TYPE<Matrix<T,3,1>,RET>;

#define PARSER_TYPE_DEFINITION(RET) \
    SCALAR_TYPE_DEFINITION(PARSER_REGISTRY,float,void) \
    SCALAR_TYPE_DEFINITION(PARSER_REGISTRY,double,void)

PARSER_TYPE_DEFINITION(void)

#define COMPLEX_SCALAR_TYPE_DEFINITION(TYPE,T,RET)           \
    template class TYPE<Matrix<T,3,1>,RET<Matrix<T,3,1>>>;

#define COMPLEX_PARSER_TYPE_DEFINITION(RET) \
    COMPLEX_SCALAR_TYPE_DEFINITION(PARSER_REGISTRY,float,RET)  \
    COMPLEX_SCALAR_TYPE_DEFINITION(PARSER_REGISTRY,double,RET)

COMPLEX_PARSER_TYPE_DEFINITION(AGGREGATOR)
COMPLEX_PARSER_TYPE_DEFINITION(PREDICATE)
