#ifndef __HASHING__
#define __HASHING__

#include <functional>

namespace std
{
template<class T>
inline void hash_combine(std::size_t& seed,T const& v)
{
    seed^=std::hash<T>()(v)+0x9e3779b9+(seed<<6)+(seed>>2);
}

template<class A,class B>
class hash<std::pair<A,B>>
{
  public:
    size_t operator()(const pair<A,B>& p) const {
        size_t h=std::hash<A>()(p.first);
        hash_combine(h,p.second);
        return h;
    }
};

template<class Tuple,size_t Index=tuple_size<Tuple>::value-1>
struct HashImplementation
{
    static void apply(size_t& seed,Tuple const& tt){
        HashImplementation<Tuple,Index-1>::apply(seed,tt);
        hash_combine(seed,std::get<Index>(tt));
    }
};

template<class Tuple>
struct HashImplementation<Tuple,0>
{
    static void apply(size_t& seed,Tuple const& tt){
        hash_combine(seed,std::get<0>(tt));
    }
};

// stackoverflow 7110301
template<typename ... TT>
struct hash<tuple<TT...>>
{
    size_t operator()(tuple<TT...> const& tt) const {
        size_t seed=0;
        HashImplementation<tuple<TT...>>::apply(seed,tt);
        return seed;
    }
};

template<class A,size_t n>
class hash<array<A,n>>
{
  public:
    size_t operator()(const array<A,n>& p) const {
        std::size_t hash_result=std::hash<A>()(p[0]);
        for(int i=1;i<n;i++){hash_combine(hash_result,p[i]);}
        return hash_result;
    }
};

template<class A>
class hash<array<A,2>>
{
  public:
    size_t operator()(const array<A,2>& p) const {
        std::size_t hash_result=std::hash<A>()(p[0]);
        hash_combine(hash_result,p[1]);
        return hash_result;
    }
};

}

#endif
