#pragma once
#include "finitefield.hpp"
#include <set>
#include <vector>
#include <array>
#include <algorithm>
#include <initializer_list>
#include <map>
#include <cstdio>

template <size_t M,class F>
class Poly{
public:
    using Term=std::array<int,M>;
    using P = Poly<M,F>;
    std::map<Term,F> f;
    Poly():f(){}
    Poly(std::initializer_list<std::pair<const Term,F>> l):f(l){}
    Poly(const P&) = default;
    Poly(P&&) = default;
    P& operator=(const P&) = default;
    P& operator=(P&&) = default;
    Poly(std::map<Term,F>&& x):f(x){};
    Poly(const std::map<Term,F>& x):f(x){};
    bool operator==(const P& r)const{
        //Assert that no value in map is zero
        return f==r.f;
    }
    bool operator!=(const P& r)const{
        return !(*this == r);
    }
    P operator-()const{
        auto t = *this;
        for(auto& p:t.f){
            p.second=-p.second;
        }
        return t;
    }
    P operator-(const P& r)const{
        return *this+(-r);
    }
    P& operator-=(const P& r){
        *this = *this-r;
        return *this;
    }
    P operator+(P&& r)const{
        for(const auto& [k,v]:f){
            r.f[k]+=v;
        }
        clear_zero(r.f);
        return r;
    }
    P operator+(const P& r)const{
        auto t = r;
        return *this+std::move(t);
    }
    P& operator+=(const P& r){
        for(const auto& [k,v]:r.f){
            f[k]+=v;
        }
        clear_zero(f);
        return *this;
    }
    P operator*(const P& r)const{
        std::map<Term,F> m;
        for(const auto& [x,xv]:f){
            for(const auto& [y,yv]:r.f){
                auto tmp = x;
                for(size_t i = 0 ; i < M ; i ++){
                    tmp[i]+=y[i];
                }
                m[tmp]+=xv*yv;
            }
        }
        clear_zero(m);
        return P{std::move(m)};
    }
    P& operator*=(const P& r){
        *this=*this*r;
        return *this;
    }
    F eval(const std::array<F,M>& val)const{
        F res{};
        for(const auto& [t,s]:f){
            res+=evalTerm(t,val)*s;
        }
        return res;
    }
private:
    F evalTerm(const Term& x,const std::array<F,M>& val)const{
        F tmp{F::one};
        for(size_t i = 0 ; i < M ; i ++){
            tmp*=val[i].pow(x[i]);
        }
        return tmp;
    }
    void clear_zero(std::map<Term,F> &m)const{
        for(auto it = m.begin() ; it != m.end() ;){
            if(it->second.iszero()){
                it = m.erase(it);
            }
            else it ++;
        }
    }
};

template <size_t M,class F>
std::ostream& operator<<(std::ostream& os,const Poly<M,F>& x){
    for(auto [t,s]:x.f){
        os << s << " * ";
        for(int i = 0 ; i < M ; i ++){
            os << "x" << i << "^" << t[i];
            if(i!=M-1)os << " * ";
            else os << " ";
        }
    }
    return os;
}
