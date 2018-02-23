#pragma once
#include <bitset>
#include <algorithm>
#include <type_traits>
#include <ostream>
#include <cassert>
#include <string>
#include <utility>
template <size_t N,const char* IRR>
class GF2n{
public:
    static constexpr std::bitset<2*N> zero{0};
    static constexpr std::bitset<2*N> one{1};
    //this is [0..n-1] is 1
    static constexpr std::bitset<2*N> all=[](){
        std::bitset<2*N> tmp;
        for(size_t i = 0 ; i < N ; i ++){
            tmp[i]=1;
        }
        return tmp;
    };
    static const std::bitset<2*N> irr;
    std::bitset<2*N> x;
    GF2n<N,IRR>(const std::string& y):x(y){}
    GF2n<N,IRR>(const std::bitset<N>& y):x(convert(y)){}
    GF2n<N,IRR>(const std::bitset<2*N>& y):x(y){}
    GF2n<N,IRR>(const GF2n<N,IRR>&)=default;
    GF2n<N,IRR>(GF2n<N,IRR>&&)=default;
    GF2n<N,IRR>& operator=(const GF2n<N,IRR>&)=default;
    GF2n<N,IRR>& operator=(GF2n<N,IRR>&&)=default;
    GF2n<N,IRR> operator+(const GF2n<N,IRR>& y)const{
        std::bitset<2*N> tmp=x^y.x;
        return GF2n<N,IRR>(tmp);
    }
    GF2n<N,IRR> operator*(const GF2n<N,IRR>& y)const{
        assert(irr==y.irr);
        return GF2n<N,IRR>(bsmul(x,y.x));
    }
    GF2n<N,IRR> operator/(const GF2n<N,IRR>& y)const{
        return GF2n<N,IRR>(bsdiv(x,y.x));
    }
    bool operator==(const GF2n<N,IRR>& y)const{
        return x==y.x;
    }
    bool operator!=(const GF2n<N,IRR>& y)const{
        return !((*this)==y);
    }
private:
    std::bitset<2*N> bsmul(const std::bitset<2*N>& l,const std::bitset<2*N>& r)const{
        std::bitset<2*N> tmp;
        for(size_t i = N-1 ;; i--){
            tmp<<=1;
            if(l[i])tmp^=r;
            if(i==0)break;
        }
        for(size_t i = 2*N-2 ; i >= N ; i --){
            if(tmp[i]){
                tmp^=irr<<(i-N);
            }
        }
        return tmp;
    }
    std::bitset<2*N> bsdiv(const std::bitset<2*N>& l,const std::bitset<2*N>& r)const{
        auto k=inv(r);
        return bsmul(l,k);
    }
    template<size_t A>
    std::bitset<2*N> convert(const std::bitset<A>& k)const{
        const size_t S=std::min(A,2*N);
        std::bitset<2*N> out;
        for(size_t i = 0 ; i < S ; i ++){
            out[i]=k[i];
        }
        return out;
    }
    size_t bslen(const std::bitset<2*N>& bs)const{
        size_t len = N+1;
        for(size_t i = N ;; i --){
            if(bs[i])break;
            len--;
            if(i==0)break;
        }
        return len;
    }
    std::pair<std::bitset<2*N>,std::bitset<2*N>> bsdivqr(std::bitset<2*N>&& r,const std::bitset<2*N>& d)const{
        //calulating p = d*q+r
        //because r calculating the result of substraction of p
        //so we just read p as r
        size_t plen = bslen(r);
        size_t dlen = bslen(d);
        assert(d!=0);
        if(plen<dlen){
            return {zero,r};
        }
        std::bitset<2*N> q;
        for(size_t i = plen-dlen ;; i --){
            if(r[i+dlen-1]){
                q[i]=1;
                r^=d<<i;
            }
            if(i==0)break;
        }
        return {q,r};
    }
    std::pair<std::bitset<2*N>,std::bitset<2*N>> bsdivqr(const std::bitset<2*N>& tr,const std::bitset<2*N>& d)const{
        std::bitset<2*N> r(tr);
        return bsdivqr(std::move(r),d);
    }
    std::bitset<N*2> inv(const std::bitset<2*N>& k)const{
        if(k==one)return k;
        auto [q,r] = bsdivqr(irr,k);
        return bsmul(inv(r),q);
    }
};
template<size_t N,const char* IRR>
const std::bitset<2*N> GF2n<N,IRR>::irr = [](){
    const std::bitset<2*N> tmp(IRR);
    assert(tmp[N]==1);
    for(size_t i = N+1 ; i < 2*N ; i ++){
        assert(tmp[i]==0);
    }
    return tmp;
}();
template<size_t N,const char* IRR>
std::ostream& operator<<(std::ostream& os,const GF2n<N,IRR> x){
    for(size_t i = N-1 ;; i --){
        os << x.x[i];
        if(i==0)break;
    }
    return os;
}
// Using elliptic curve with characteristic=2 , and j(E) != 0
// So the equation will be y^2 + xy = x^3 + a2x^2 + a6
template<class F>
class EC{
    F x,y;
    bool inf;
    const F a2,a6;
    EC<F>(const EC<F>&)=default;
    EC<F>(EC<F>&&)=default;
    EC<F>& operator=(const EC<F>&)=default;
    EC<F>& operator=(EC<F>&&)=default;
    EC<F>(const F& ix,const F& iy,bool iinf,const F& ia2,const F& ia6):x(ix),y(iy),inf(iinf),a2(ia2),a6(ia6){
        assert(a6.x!=F::zero);
        if(!inf){
            assert(y*y+x*y==x*x*x+a2*x*x+a6);
        }
    }
    bool operator==(const EC<F>& r)const{
        assert(a2==r.a2&&a6==r.a6);
        if(inf!=r.inf)return false;
        if(inf) return true;
        return x==r.x&&y==r.y;
    }
    EC<F> operator-()const{
        return EC<F>(x,x+y,inf,a2,a6);
    }
    EC<F> operator+(const EC<F>& r)const{
        assert(a2==r.a2&&a6==r.a6);
        if(r.inf){
            return *this;
        }
        if(inf){
            return r;
        }
        if(*this==-r){
            return EC<F>(x,y,true,a2,a6);
        }
        if(*this==r){
        }
        else{
        }
    }

};
