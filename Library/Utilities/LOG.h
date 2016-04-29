#ifndef __LOG__
#define __LOG__

#include <iostream>

namespace Mechanics{

class NULLBUFFER:public std::streambuf
{
public:
    int overflow(int c){return c;}
};

class LOG
{
public:
    static void Output(bool print){
        if(print){LOG::cout.rdbuf(std::cout.rdbuf());}
        else{LOG::cout.rdbuf(new NULLBUFFER());}
    }

    static std::ostream cout;
};

}
#endif
