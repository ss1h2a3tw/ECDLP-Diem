#pragma once
#include <array>
#include <vector>
#include <cstdio>
#include "poly.hpp"

template<size_t K,size_t maxDeg,class P>
std::array<Poly<K*P::m,GF2>,P::F::n> weilDescent (const P& p){
    using std::array;
    using std::vector;
    using F = typename P::F;
    const constexpr size_t n = P::F::n;
    const constexpr size_t m = P::m;
    using W = array<Poly<K*m,GF2>,n>;
    const static array<vector<int>,2*n-1>multab = [](){
        array<vector<int>,2*n-1> res{};
        F base[n]={};
        for(size_t i = 0 ; i < n ; i ++){
            base[i].x[i]=1;
            res[i].push_back(i);
        }
        for(size_t i = 1 ; i < n ; i ++){
            auto tmp=base[n-1]*base[i];
            for(size_t j = 0 ; j < n ; j ++){
                if(tmp.x[j]){
                    res[n-1+i].push_back(j);
                }
            }
        }
        return res;
    }();
    auto mul = [](const W& l,const W& r){
        W res{};
        for(size_t i = 0 ; i < n ; i ++){
            for(size_t j = 0 ; j < n ; j ++){
                auto tmp = l[i]*r[j];
                Poly<K*m,GF2> af{};
                for(const auto& [t,s] : tmp.f){
                    auto tt =t;
                    for(size_t idx = 0 ; idx < K*m ; idx ++){
                        if(tt[idx])tt[idx]=1;
                    }
                    af.addTerm(tt,s);
                }
                af.clearZero();
                for(int idx:multab[i+j]){
                    res[idx]+=af;
                }
            }
        }
        return res;
    };
    auto makeVar = [](size_t x){
        W res{};
        for(size_t i = 0 ; i < K ; i ++){
            array<int,K*m> tmp{};
            tmp[x*K+i]=1;
            res[i].addTerm(tmp,GF2{true});
        }
        return res;
    };
    auto convScaler = [](const F& s){
        W res{};
        array<int,K*m> tmp{};
        for(size_t i = 0 ; i < n ; i ++){
            if(s.x[i]){
                res[i].addTerm(tmp,GF2{true});
            }
        }
        return res;
    };
    static const array<array<W,maxDeg+1>,m> varTab = [&](){
        const W one{Poly<K*m,GF2>{{array<int,K*m>{},GF2{true}}}};
        array<array<W,maxDeg+1>,m> res;
        for(size_t i = 0 ; i < m ; i ++){
            res[i][0]=one;
        }
        for(size_t i = 0 ; i < m ; i ++){
            res[i][0]=one;
            res[i][1]=makeVar(i);
            for(size_t j = 2 ; j <= maxDeg ; j ++){
                res[i][j]=mul(res[i][j-1],res[i][1]);
            }
        }
        return res;
    }();
    W ret{};
    int cnt=0;
    for(const auto [t,s]:p.f){
        auto tmp = convScaler(s);
        for(size_t i = 0 ; i < m ; i ++){
            tmp=mul(tmp,varTab[i][t[i]]);
        }
        for(size_t i = 0 ; i < n ; i ++){
            ret[i]+=tmp[i];
        }
        cnt++;
    }
    return ret;
}
