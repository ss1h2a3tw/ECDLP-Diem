#pragma once
#include <queue>
#include <utility>
#include <vector>
#include <algorithm>
#include <cassert>
template<class P>
std::vector<P> grobnerGF2(std::vector<P> pv){
    using std::queue;
    using std::pair;
    using std::array;
    using std::max;
    constexpr size_t m = P::m;
    queue<pair<size_t,size_t>> q;
    for(size_t i = 0 ; i < pv.size() ; i ++){
        for(size_t j = i+1 ; j < pv.size() ; j ++){
            q.push({i,j});
        }
    }
    while(!q.empty()){
        auto [x,y] = q.front();
        q.pop();
        const auto& l = pv[x];
        const auto& r = pv[y];
        array<int,m> lcm,ml,mr;
        for(size_t i = 0 ; i < m ; i ++){
            lcm[i]=max(l.f.rbegin()->first[i],r.f.rbegin()->first[i]);
            ml[i]=lcm[i]-l.f.rbegin()->first[i];
            mr[i]=lcm[i]-r.f.rbegin()->first[i];
        }
        P aml,amr;
        for(const auto& [t,s]:l.f){
            auto tt = t;
            for(size_t i = 0 ; i < m ; i ++){
                tt[i]+=ml[i];
                if(tt[i])tt[i]=1;
            }
            aml.addTerm(tt,s);
        }
        for(const auto& [t,s]:r.f){
            auto tt = t;
            for(size_t i = 0 ; i < m ; i ++){
                tt[i]+=mr[i];
                if(tt[i])tt[i]=1;
            }
            amr.addTerm(tt,s);
        }
        aml.clearZero();
        amr.clearZero();
        aml-=amr;
        for(const auto& tp:pv){
            if(aml.f.size()==0)break;
            assert(tp.f.size()>0);
            bool nl=false;
            array<int,m> mp;
            for(size_t i = 0 ; i < m ; i ++){
                if(aml.f.rbegin()->first[i]<tp.f.rbegin()->first[i])nl=true;
                mp[i]=aml.f.rbegin()->first[i]-tp.f.rbegin()->first[i];
            }
            if(nl)continue;
            P ttp;
            for(const auto& [t,s]:tp.f){
                auto tt=t;
                for(size_t i = 0 ; i < m ; i ++){
                    tt[i]+=mp[i];
                    if(tt[i])tt[i]=1;
                }
                ttp.addTerm(tt,s);
            }
            ttp.clearZero();
            aml-=ttp;
        }
        if(aml.f.size()>0){
            pv.push_back(aml);
            for(size_t i = 0 ; i < pv.size()-1 ; i ++){
                q.push({i,pv.size()-1});
            }
        }
    }
    return pv;
}
