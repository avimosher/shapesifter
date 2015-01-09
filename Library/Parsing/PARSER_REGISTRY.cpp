//////////////////////////////////////////////////////////////////////
// Copyright 2014, Avi Robinson-Mosher.
///////////////////////////////////////////////////////////////////////
// Class PARSER_REGISTRY
///////////////////////////////////////////////////////////////////////
#include <Parsing/PARSER_REGISTRY.h>
#include <Utilities/TYPE_UTILITIES.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV,class PARSED_TYPE> std::map<std::string,PARSED_TYPE* (*)(Json::Value&,SIMULATION<TV>&)>& PARSER_REGISTRY<TV,PARSED_TYPE>::
Parsers()
{
    static std::map<std::string,Parse_Function> parser_registry;
    return parser_registry;
}
///////////////////////////////////////////////////////////////////////
#define SCALAR_TYPE_DEFINITION(TYPE,T,RET)           \
    template class TYPE<Matrix<T,1,1>,RET>;          \
    template class TYPE<Matrix<T,2,1>,RET>;          \
    template class TYPE<Matrix<T,3,1>,RET>;

SCALAR_TYPE_DEFINITION(PARSER_REGISTRY,float,void)
SCALAR_TYPE_DEFINITION(PARSER_REGISTRY,double,void)
