#include <fstream>
#include <string>
#include <cerrno>
#include <iostream>
#include <Eigen/Dense>
#include <json/json.h>

using Eigen::MatrixXd;

std::string get_file_contents(const char* filename)
{
    std::ifstream in(filename,std::ios::in|std::ios::binary);
    if(in) {
        std::string contents;
        in.seekg(0,std::ios::end);
        contents.reserve(in.tellg());
        in.seekg(0,std::ios::beg);
        contents.assign((std::istreambuf_iterator<char>(in)),std::istreambuf_iterator<char>());
        in.close();
        return(contents);
    }
    throw(errno);
}

int main()
{
    MatrixXd m(2,2);
    m(0,0)=3;
    m(1,0)=2.5;
    m(0,1)=-1;
    m(1,1)=m(1,0)+m(0,1);
    std::cout<<m<<std::endl;

    Json::Value root;
    Json::Reader reader;
    bool parsingSuccess=reader.parse(get_file_contents("test.json"),root);
    std::cout<<parsingSuccess<<std::endl;
    const Json::Value objects=root["objects"];
    for (int index=0;index<objects.size();++index){
        std::cout<<objects[index].asString()<<std::endl;
    }
}
