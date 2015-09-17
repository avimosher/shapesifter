#ifndef __CIRCULAR_STACK__
#define __CIRCULAR_STACK__

namespace Mechanics{

template<class T>
class circular_stack
{
    std::vector<T> data;
    int top_index;
public:
    circular_stack(int size):data(size),top_index(-1)
    {}

    void push(const T& item)
    {data[++top_index%data.size()]=item;}

    T pop()
    {return data[top_index--%data.size()];}
};
}
#endif
