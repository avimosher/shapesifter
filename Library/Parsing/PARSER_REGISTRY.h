//#####################################################################
// Copyright 2014, Avi Robinson-Mosher.
//#####################################################################
// Class PARSER_REGISTRY
//#####################################################################
#ifndef __PARSER_REGISTRY__
#define __PARSER_REGISTRY__

#include <json/json.h>
#include <map>
#include <iostream>

namespace Mechanics{
template<class TV> class SIMULATION;

template<class TV>
class PARSER_REGISTRY
{
    typedef void (*Parse_Function)(Json::Value&,SIMULATION<TV>&);
    static std::map<std::string,Parse_Function>& Parsers();
public:
    PARSER_REGISTRY(){};
    ~PARSER_REGISTRY(){};

    static void Parse(Json::Value& node,SIMULATION<TV>& data)
    {std::cout<<node["type"].asString()<<std::endl;
        (*Parsers()[node["type"].asString()])(node,data);}
    template<class T_PARSER>static void Register();
};

template<class TV> template<class T_PARSER> void PARSER_REGISTRY<TV>::
Register()
{
    Parsers()[T_PARSER::Static_Name()]=T_PARSER::Parse;
}

#define REGISTER_PARSER_SCALAR(TYPE,T,d) PARSER_REGISTRY<Matrix<T,1,d>>::Register<TYPE<Matrix<T,1,d>>>();
#define REGISTER_PARSER_GENERIC(TYPE,T) REGISTER_PARSER_SCALAR(TYPE,T,1);REGISTER_PARSER_SCALAR(TYPE,T,2);REGISTER_PARSER_SCALAR(TYPE,T,3);
#define REGISTER_PARSER(TYPE) \
    namespace Mechanics{ \
    bool Register_##TYPE##_Parser() \
    { \
        REGISTER_PARSER_GENERIC(TYPE,float);REGISTER_PARSER_GENERIC(TYPE,double); \
    return true; \
    } \
    static bool registered=Register_##TYPE##_Parser(); \
    };

#define DEFINE_AND_REGISTER_PARSER(TYPE)                                \
    namespace Mechanics{                                                \
    template<class TV> class PARSE_##TYPE                               \
    {                                                                   \
      public:                                                           \
        PARSE_##TYPE(){};                                               \
            static void Parse(Json::Value& node,SIMULATION<TV>& simulation); \
            static std::string Static_Name(){                           \
                return TYPE<TV>::Static_Name();                         \
            }                                                           \
    };                                                                  \
    }                                                                   \
    GENERIC_TYPE_DEFINITION(PARSE_##TYPE)                               \
    REGISTER_PARSER(PARSE_##TYPE)                                       \
    template<class TV> void PARSE_##TYPE<TV>::Parse(Json::Value& node,SIMULATION<TV>& simulation)


}
#endif
