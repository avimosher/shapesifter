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

template<class T>
inline void hash_combine(std::size_t& seed,T const& v)
{
    seed^=std::hash<T>()(v)+0x9e3779b9+(seed<<6)+(seed>>2);
}

template<class A,class B,class C>
class hash<tuple<A,B,C>>
{
  public:
    size_t operator()(const tuple<A,B,C>& p) const {
        std::size_t hash_result=std::hash<A>()(std::get<0>(p));
        hash_combine(hash_result,std::get<1>(p));
        hash_combine(hash_result,std::get<2>(p));
        return hash_result;
    }
};    
}

#endif
