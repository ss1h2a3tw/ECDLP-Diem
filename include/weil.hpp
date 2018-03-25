#pragma once
#include <array>
#include <vector>
#include <cstdio>
#include "poly.hpp"

template<size_t K,class P>
std::array<Poly<K*P::m,GF2>,P::F::n> weilDescent (const P& p){
    using std::array;
    using std::vector;
    using F = typename P::F;
    const constexpr size_t n = P::F::n;
    const constexpr size_t m = P::m;
    const static array<array<vector<int>,n>,n>multab = [](){
        array<array<vector<int>,n>,n> res{};
        F base[n]={};
        for(size_t i = 0 ; i < n ; i ++){
            base[i].x[i]=1;
        }
        for(size_t i = 0 ; i < n ; i ++){
            for(size_t j = 0 ; j < n ; j ++){
                auto tmp=base[i]*base[j];
                for(size_t k = 0 ; k < n ; k ++){
                    if(tmp.x[k]){
                        res[i][j].push_back(k);
                    }
                }
            }
        }
        return res;
    }();
    auto mul = [](const array<Poly<K*m,GF2>,n>& l,const array<Poly<K*m,GF2>,n>& r){
        //printf("mul\n");
        array<Poly<K*m,GF2>,n> res{};
        for(size_t i = 0 ; i < n ; i ++){
            for(size_t j = 0 ; j < n ; j ++){
                auto tmp = l[i]*r[j];
                Poly<K*m,GF2> af{};
                for(const auto& [t,s] : tmp.f){
                    auto tt =t;
                    for(size_t idx = 0 ; idx < K*m ; idx ++){
                        if(tt[i])tt[i]=1;
                    }
                    af.addTerm(tt,s);
                }
                af.clearZero();
                for(int idx:multab[i][j]){
                    res[idx]+=af;
                }
            }
        }
        //printf("Finish mul\n");
        return res;
    };
    auto pow = [mul](const array<Poly<K*m,GF2>,n>& l,int x){
        array<Poly<K*m,GF2>,n> res{};
        // make res = 1
        res[0].addTerm(array<int,K*m>{},GF2{true});
        auto now = l;
        while(x){
            if(x&1){
                res=mul(res,now);
            }
            now = mul(now,now);
            x>>=1;
        }
        return res;
    };
    auto makeVar = [](size_t x){
        //printf("makeVar\n");
        array<Poly<K*m,GF2>,n> res{};
        for(size_t i = 0 ; i < K ; i ++){
            array<int,K*m> tmp{};
            //printf("!!m%d x%d K%d i%d = %d\n",(int)m,(int)x,(int)K,(int)i,(int)(x*K+i));
            tmp[x*K+i]=1;
            res[i].addTerm(tmp,GF2{true});
        }
        //printf("Finish makeVar\n");
        return res;
    };
    auto convScaler = [](const F& s){
        array<Poly<K*m,GF2>,n> res{};
        array<int,K*m> tmp{};
        for(size_t i = 0 ; i < n ; i ++){
            if(s.x[i]){
                res[i].addTerm(tmp,GF2{true});
            }
        }
        return res;
    };
    array<Poly<K*m,GF2>,n> ret{};
    int cnt=0;
    for(const auto [t,s]:p.f){
        auto tmp = convScaler(s);
        for(size_t i = 0 ; i < m ; i ++){
            tmp=mul(tmp,pow(makeVar(i),t[i]));
        }
        for(size_t i = 0 ; i < n ; i ++){
            ret[i]+=tmp[i];
        }
        cnt++;
    }
    return ret;
}
