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
    Poly(std::initializer_list<std::pair<const Term,F>> l):f(l){}
    Poly(const P&) = default;
    Poly(P&&) = default;
    Poly(std::map<Term,F>&& x):f(x){};
    Poly(const std::map<Term,F>& x):f(x){};
    P& operator=(P&&) = default;
    bool operator==(const P& r)const{
        //Assert that no value in map is zero
        return f==r.f;
    }
    P operator+(P&& r)const {
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
private:
    void clear_zero(std::map<Term,F> &m)const{
        for(auto it = m.begin() ; it != m.end() ;){
            if(it->second.iszero()){
                it = m.erase(it);
            }
            else it ++;
        }
    }
};
