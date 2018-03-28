#pragma once
#include <bitset>
#include <algorithm>
#include <type_traits>
#include <ostream>
#include <cassert>
#include <string>
#include <utility>
#include <iostream>

struct GF2{
    bool x;
    static constexpr bool one{true};
    GF2 operator-()const;
    GF2 operator-(const GF2& r)const;
    GF2& operator-=(const GF2&r);
    GF2 operator+(const GF2& r)const;
    GF2& operator+=(const GF2&r);
    GF2 operator*(const GF2& r)const;
    GF2& operator*=(const GF2& r);
    GF2 sqrt()const;
    GF2 pow(unsigned long long p)const;
    bool operator==(const GF2& r)const;
    bool operator!=(const GF2& r)const;
    bool iszero()const;
};
std::ostream& operator<<(std::ostream& os,const GF2 x);

template <size_t N,const char* IRR>
class GF2n{
public:
    using F = GF2n<N,IRR>;
    static constexpr size_t n = N;
    static constexpr std::bitset<2*N> zero{0};
    static constexpr std::bitset<2*N> one{1};
    //this is [0..n-1] is 1
    static constexpr std::bitset<2*N> all=[](){
        std::bitset<2*N> tmp;
        for(size_t i = 0 ; i < N ; i ++){
            tmp[i]=1;
        }
        return tmp;
    }();
    static const std::bitset<2*N> irr;
    std::bitset<2*N> x;
    GF2n<N,IRR>():x(){}
    GF2n<N,IRR>(const std::string& y):x(y){}
    GF2n<N,IRR>(const std::bitset<N>& y):x(convert(y)){}
    GF2n<N,IRR>(const std::bitset<2*N>& y):x(y){}
    GF2n<N,IRR>(std::bitset<2*N>&& y):x(std::move(y)){}
    GF2n<N,IRR>(unsigned long long y):x(y){};
    GF2n<N,IRR>(const F&)=default;
    GF2n<N,IRR>(F&&)=default;
    F& operator=(const F&)=default;
    F& operator=(F&&)=default;
    F operator-(const F& y)const{
        return *this+y;
    }
    F operator-(F&& y)const{
        return *this+std::move(y);
    }
    F operator-()const{
        return *this;
    }
    F& operator-=(const F& y){
        *this+=y;
        return *this;
    }
    F operator+(const F& y)const{
        std::bitset<2*N> tmp=x^y.x;
        return F(tmp);
    }
    F& operator+=(const F& y){
        x^=y.x;
        return *this;
    }
    F operator*(const F& y)const{
        return F(bsmul(x,y.x));
    }
    F& operator*=(const F& y){
        *this=*this*y;
        return *this;
    }
    F operator/(const F& y)const{
        return F(bsdiv(x,y.x));
    }
    F& operator/=(const F& y){
        *this = *this/y;
        return *this;
    }
    F pow(unsigned long long p)const{
        auto t=*this;
        auto now=F{one};
        while(p){
            if(p&1){
                now*=t;
            }
            p>>=1;
            t=t*t;
        }
        return now;
    }
    F sqrt()const{
        auto t = *this;
        for(size_t i = 0 ; i < N-1 ; i ++){
            t=t*t;
        }
        return t;
    }
    F trace()const{
        F t=x;
        F res{};
        for(size_t i = 0 ; i < n ; i ++){
            res+=t;
            t=t*t;
        }
        return res;
    }
    F halftrace()const{
        F t=x;
        F res{};
        for(size_t i = 0 ; i <= (n-1)/2 ; i ++){
            res+=t;
            t=t*t;
            t=t*t;
        }
        return res;
    }
    bool operator==(const F& y)const{
        return x==y.x;
    }
    bool operator!=(const F& y)const{
        return !((*this)==y);
    }
    bool iszero()const{
        return x==zero;
    }
private:
    constexpr std::bitset<2*N> bsmul(const std::bitset<2*N>& l,const std::bitset<2*N>& r)const{
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
    constexpr std::bitset<2*N> bsdiv(const std::bitset<2*N>& l,const std::bitset<2*N>& r)const{
        const auto k=inv(r);
        return bsmul(l,k);
    }
    template<size_t A>
    constexpr std::bitset<2*N> convert(const std::bitset<A>& k)const{
        const size_t S=std::min(A,2*N);
        std::bitset<2*N> out;
        for(size_t i = 0 ; i < S ; i ++){
            out[i]=k[i];
        }
        return out;
    }
    constexpr size_t bslen(const std::bitset<2*N>& bs)const{
        size_t len = N+1;
        for(size_t i = N ;; i --){
            if(bs[i])break;
            len--;
            if(i==0)break;
        }
        return len;
    }
    constexpr std::pair<std::bitset<2*N>,std::bitset<2*N>> bsdivqr(std::bitset<2*N>&& r,const std::bitset<2*N>& d)const{
        //calulating p = d*q+r
        //because r calculating the result of substraction of p
        //so we just read p as r
        const size_t plen = bslen(r);
        const size_t dlen = bslen(d);
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
    constexpr std::pair<std::bitset<2*N>,std::bitset<2*N>> bsdivqr(const std::bitset<2*N>& tr,const std::bitset<2*N>& d)const{
        std::bitset<2*N> r(tr);
        return bsdivqr(std::move(r),d);
    }
    constexpr std::bitset<N*2> inv(const std::bitset<2*N>& k)const{
        assert(k!=zero);
        if(k==one)return k;
        const auto [q,r] = bsdivqr(irr,k);
        return bsmul(inv(r),q);
    }
};
template<size_t N,const char* IRR>
const std::bitset<2*N> __attribute__ ((init_priority (200))) GF2n<N,IRR>::irr = [](){
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
template<class Field,const Field& A2,const Field& A6>
class EC{
public:
    using F=Field;
    using E=EC<F,A2,A6>;
    static const F &a2,&a6;
    F x,y;
    bool inf;
    EC<F,A2,A6>():x(),y(),inf(true){}
    EC<F,A2,A6>(const E&)=default;
    EC<F,A2,A6>(E&&)=default;
    EC<F,A2,A6>& operator=(const E&)=default;
    EC<F,A2,A6>& operator=(E&&)=default;
    template <typename J,typename K>
    EC<F,A2,A6>(J && ix,K && iy,bool iinf):x(std::forward<J>(ix)),y(std::forward<K>(iy)),inf(iinf){
        if(!inf){
            assert(y*y+x*y==x*x*x+a2*x*x+a6);
        }
    }
    static F gety(const F& x){
        static const F goodbasis=[](){
            F tmp{};
            for(size_t i = 0 ; i < F::n ; i ++){
                tmp.x[i]=1;
                if(tmp.trace().x==F::one){
                    break;
                }
                tmp.x[i]=0;
            }
            return tmp;
        }();
        //if can't find, throw 0
        if(x.x==F::zero){
            return a6.sqrt();
        }
        F b = x;
        F c = -(x*x*x + a2*x*x + a6);
        F tb = c/(b*b);
        if(tb.trace().x!=F::zero){
            throw 0;
        }
        if(F::n&1){
            //odd
            tb = tb.halftrace();
            tb *= b;
            return tb;
        }
        else{
            //even
            F res{};
            for(size_t i=0 ; i <= F::n-2 ; i ++){
                F d=goodbasis;
                for(size_t j = 0 ; j <= i ; j ++){
                    d=d*d;
                }
                F mul={};
                for(size_t j = i + 1 ; j <= F::n-1 ; j ++){
                    mul+=d;
                    d=d*d;
                }
                res += mul*tb;
                tb=tb*tb;
            }
            res *= b;
            return res;
        }
    }
    bool operator==(const E& r)const{
        if(inf!=r.inf)return false;
        if(inf) return true;
        return x==r.x&&y==r.y;
    }
    E operator-()const{
        if(inf)return *this;
        return E(x,x+y,inf);
    }
    E operator-(const E& r)const{
        return *this+(-r);
    }
    E& operator-=(const E& r){
        *this = *this-r;
        return *this;
    }
    E operator+(const E& r)const{
        assert(a2==r.a2&&a6==r.a6);
        if(r.inf){
            return *this;
        }
        if(inf){
            return r;
        }
        if(*this==-r){
            return E(x,y,true);
        }
        if(*this==r){
            const auto tmp = x*x;
            auto x3 = tmp+a6/tmp;
            auto y3 = tmp + (x + y/x)*x3 + x3;
            return E(std::move(x3),std::move(y3),false);
        }
        //*this!=r,this!=-r,none of them is inf
        const auto lam = (y+r.y)/(x+r.x);
        auto x3 = lam*lam + lam + x + r.x + a2;
        auto y3 = lam*(x+x3) + x3 + y;
        return E(std::move(x3),std::move(y3),false);
    }
    E& operator+=(const E& r){
        *this = *this+r;
        return *this;
    }
    //Need to add bignum support later
    E operator*(unsigned long long k){
        auto ans = E(F::zero,F::zero,true);
        auto x = *this;
        while(k){
            if(k&1){
                ans+=x;
            }
            k>>=1;
            x=x+x;
        }
        return ans;
    }
};
template<class F,const F& A2,const F& A6>
const F& EC<F,A2,A6>::a2 = A2;
template<class F,const F& A2,const F& A6>
const F& EC<F,A2,A6>::a6 = [](){
    assert(A6!=F::zero);
    return A6;
}();
template<class F,const F& A2,const F& A6>
std::ostream& operator<<(std::ostream& os,const EC<F,A2,A6> x){
    os << '(' << x.x << ',' << x.y << ')' ;
    return os;
}
