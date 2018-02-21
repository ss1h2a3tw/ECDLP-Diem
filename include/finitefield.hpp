#pragma once
#include <bitset>
#include <algorithm>
#include <type_traits>
#include <ostream>
template <size_t N>
class GF2n{
public:
    std::bitset<2*N> x;
    const std::bitset<2*N> irr;
    GF2n<N>(const std::bitset<N> y,const std::bitset<N+1>& yirr):x(convert(y)),irr(convert(yirr)){}
    GF2n<N>(const std::bitset<N> y,const std::bitset<2*N>& yirr):x(convert(y)),irr(yirr){}
    GF2n<N>(const std::bitset<2*N> y,const std::bitset<N+1>& yirr):x(y),irr(convert(yirr)){}
    GF2n<N>(const std::bitset<2*N> y,const std::bitset<2*N>& yirr):x(y),irr(yirr){}
    GF2n<N>(const GF2n<N>&)=default;
    GF2n<N>(GF2n<N>&&)=default;
    GF2n<N>& operator=(const GF2n<N>&)=default;
    GF2n<N>& operator=(GF2n<N>&&)=default;
    GF2n<N> operator+(const GF2n<N> y)const{
        std::bitset<2*N> tmp=x^y.x;
        return GF2n<N>(tmp,irr);
    }
    GF2n<N> operator*(const GF2n<N> y)const{
        std::bitset<2*N> tmp;
        for(size_t i = N-1 ;; i--){
            tmp<<=1;
            if(x[i])tmp^=y.x;
            if(i==0)break;
        }
        for(size_t i = 2*N-2 ; i >= N ; i --){
            if(tmp[i]){
                tmp^=irr<<(i-N);
            }
        }
        return GF2n<N>(tmp,irr);
    }
    bool operator==(const GF2n<N> y)const{
        return x==y.x && irr == y.irr;
    }
    bool operator!=(const GF2n<N> y)const{
        return !((*this)==y);
    }
private:
    template<size_t A>
    std::bitset<2*N> convert(std::bitset<A> k){
        const size_t S=std::min(A,2*N);
        std::bitset<2*N> out;
        for(size_t i = 0 ; i < S ; i ++){
            out[i]=k[i];
        }
        return out;
    }
};
template<size_t N>
std::ostream& operator<<(std::ostream& os,const GF2n<N> x){
    for(size_t i = N-1 ;; i --){
        os << x.x[i];
        if(i==0)break;
    }
    return os;
}
