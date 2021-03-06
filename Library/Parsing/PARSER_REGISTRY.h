#ifndef __PARSER_REGISTRY__
#define __PARSER_REGISTRY__

#include <Data/ROTATION.h>
#include <Evolution/EVOLUTION.h>
#include <Utilities/LOG.h>
#include <iostream>
#include <map>
#include <memory>
#include <json/json.h>

namespace Mechanics{
template<class TV> class SIMULATION;

template<class TV,class PARSED_TYPE>
class PARSER_REGISTRY
{
    typedef std::shared_ptr<PARSED_TYPE> (*Parse_Function)(Json::Value&,SIMULATION<TV>&);
    static std::map<std::string,Parse_Function>& Parsers();
public:
    PARSER_REGISTRY(){};
    ~PARSER_REGISTRY(){};

    static std::shared_ptr<PARSED_TYPE> Parse(Json::Value& node,SIMULATION<TV>& data)
    {
        std::string node_type=node["type"].asString();
        LOG::cout<<node_type<<std::endl;
        if(Parsers().find(node_type)==Parsers().end()){std::cout<<"Parsing failed for node type: "<<node_type<<std::endl;exit(-1);}
        return (*Parsers()[node_type])(node,data);
    }
    template<class T_PARSER>static void Register();
};

template<class TV,class PARSED_TYPE> template<class T_PARSER> void PARSER_REGISTRY<TV,PARSED_TYPE>::
Register()
{
    Parsers()[T_PARSER::Static_Name()]=T_PARSER::Parse;
}

#define REGISTER_PARSER_SCALAR(TYPE,T,d,PARSED_TYPE) PARSER_REGISTRY<Matrix<T,d,1>,PARSED_TYPE>::Register<TYPE<Matrix<T,d,1>>>();
#define REGISTER_PARSER_GENERIC(TYPE,T,PARSED_TYPE) \
    REGISTER_PARSER_SCALAR(TYPE,T,3,PARSED_TYPE);
#define REGISTER_PARSER(TYPE,PARSED_TYPE)       \
    namespace Mechanics{                        \
    bool Register_##TYPE##_Parser() \
    { \
        REGISTER_PARSER_GENERIC(TYPE,double,PARSED_TYPE); \
    return true; \
    } \
    static bool registered=Register_##TYPE##_Parser(); \
    };

#define DEFINE_AND_REGISTER_PARSER(TYPE,PARSED_TYPE)                    \
    namespace Mechanics{                                                \
    template<class TV> class PARSE_##TYPE                               \
    {                                                                   \
        typedef typename TV::Scalar T;                                  \
      public:                                                           \
        PARSE_##TYPE(){};                                               \
            static std::shared_ptr<PARSED_TYPE> Parse(Json::Value& node,SIMULATION<TV>& simulation); \
            static std::string Static_Name(){                           \
                return TYPE<TV>::Static_Name();                         \
            }                                                           \
    };                                                                  \
    }                                                                   \
    GENERIC_TYPE_DEFINITION(PARSE_##TYPE)                               \
    REGISTER_PARSER(PARSE_##TYPE,PARSED_TYPE)                                      \
    template<class TV> std::shared_ptr<PARSED_TYPE> PARSE_##TYPE<TV>::Parse(Json::Value& node,SIMULATION<TV>& simulation)

#define REGISTER_TEMPLATE_PARSER_SCALAR(TYPE,T,d,PARSED_TYPE) PARSER_REGISTRY<Matrix<T,d,1>,PARSED_TYPE<Matrix<T,d,1>>>::Register<TYPE<Matrix<T,d,1>>>();
#define REGISTER_TEMPLATE_PARSER_GENERIC(TYPE,T,PARSED_TYPE) \
    REGISTER_TEMPLATE_PARSER_SCALAR(TYPE,T,3,PARSED_TYPE);
#define REGISTER_TEMPLATE_PARSER(TYPE,PARSED_TYPE)                      \
    namespace Mechanics{                                                \
    bool Register_##TYPE##_Parser()                                     \
    {                                                                   \
        REGISTER_TEMPLATE_PARSER_GENERIC(TYPE,double,PARSED_TYPE); \
        return true;                                                    \
    }                                                                   \
    static bool TYPE##_registered=Register_##TYPE##_Parser();                  \
    };

#define DEFINE_AND_REGISTER_TEMPLATE_PARSER(TYPE,PARSED_TYPE)           \
    namespace Mechanics{                                                \
    template<class TV> class PARSE_##TYPE                               \
    {                                                                   \
      public:                                                           \
        PARSE_##TYPE(){};                                               \
            static std::shared_ptr<PARSED_TYPE<TV>> Parse(Json::Value& node,SIMULATION<TV>& simulation); \
            static std::string Static_Name(){                           \
                return TYPE<TV>::Static_Name();                         \
            }                                                           \
    };                                                                  \
    }                                                                   \
    GENERIC_TYPE_DEFINITION(PARSE_##TYPE)                               \
    REGISTER_TEMPLATE_PARSER(PARSE_##TYPE,PARSED_TYPE)                  \
    template<class TV> std::shared_ptr<PARSED_TYPE<TV>> PARSE_##TYPE<TV>::Parse(Json::Value& node,SIMULATION<TV>& simulation)

template<class T> void Parse_Scalar(Json::Value& node,T& value,const T& default_value=T())
{
    if(!node.isNull()){value=node.asDouble();}
    else{value=default_value;}
}

inline void Parse_String(Json::Value& node,std::string& value,const std::string& default_value="")
{
    if(!node.isNull()){value=node.asString();}
    else{value=default_value;}
}

template<class TV> void Parse_Vector(Json::Value& node,TV& vector,const TV& default_vector=TV())
{
    if(!node.isNull()){for(int i=0;i<vector.size();i++){vector[i]=node[i].asDouble();}}
    else{vector=default_vector;}
}

template<class TV> void Parse_Rotation(Json::Value& node,ROTATION<TV>& orientation,const ROTATION<TV>& default_orientation=ROTATION<TV>::ORIENTATION::Identity())
{
    if(!node.isNull()){
        typename TV::Scalar angle;Parse_Scalar(node["angle"],angle);
        TV axis;Parse_Vector(node["axis"],axis);
        orientation=ROTATION<TV>::From_Rotation_Vector((angle*M_PI/180)*axis.normalized());
    }
    else{orientation=default_orientation;}
}

template<class TV> void Set_Vector(const TV& vector,Json::Value& node)
{
    node=Json::Value(Json::arrayValue);
    node.resize(vector.size());
    for(int i=0;i<vector.size();i++){node[i]=vector[i];}
}

}
#endif
