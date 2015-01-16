#ifndef __HASHING__
#define __HASHING__

namespace std
{
template<class A,class B>
class hash<std::pair<A,B>>
{
  public:
    size_t operator()(const pair<A,B>& p) const {
        static hash<A> first_hash;
        static hash<B> second_hash;
        size_t h1=first_hash(p.first);
        return h1^(second_hash(p.second)+0x9e3779b9+(h1<<6)+(h1>>2));
    }
};
}

#endif
