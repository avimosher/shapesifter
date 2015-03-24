#ifndef __LOG__
#define __LOG__

#include <iostream>

namespace Mechanics{

class NULLBUFFER:public std::streambuf
{
public:
    int overflow(int c){return c;}
};

class LOG:public std::ostream
{
public:
    static void Output(bool print){
        if(print){cout.rdbuf(std::cout.rdbuf());}
        else{cout.rdbuf(new NULLBUFFER());}
    }

    static LOG cout;
};

}
#endif
