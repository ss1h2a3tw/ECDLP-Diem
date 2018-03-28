#include <iostream>
#include <random>
#include <unordered_set>
#include <unordered_map>
#include <list>
#include "finitefield.hpp"
#include "poly.hpp"
#include "semaev.hpp"
#include "weil.hpp"
#include "grobner.hpp"

using namespace std;
template <size_t N,size_t M,class F>
vector<Poly<M,F>> conv(const array<Poly<M,F>,N>& a){
    vector<Poly<M,F>> v;
    for(size_t i = 0 ; i < N ; i ++){
        if(a[i].f.size())v.push_back(a[i]);
    }
    return v;
}

const int N=10;
const int n=5,m=2;
//x^10 + x^6 + x^5 + x^3 + x^2 + x + 1
const char irr[]="10001101111";
using F = GF2n<N,irr>;
//x^9 + x^8 + x^5 + x^3 + x
const F a2 = F{"1100101010"};
//x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1
const F a6 = F{"0011111111"};
using E = EC<F,a2,a6>;

unsigned long long get_rand(unsigned long long m){
    static random_device r;
    static mt19937_64 e(r());
    uniform_int_distribution<unsigned long long> u(1,m);
    return u(e);
}
template <size_t M>
list<pair<Poly<M,GF2>,bitset<M>>> convlist(const vector<Poly<M,GF2>>& v){
    using PP = pair<Poly<M,GF2>,bitset<M>>;
    vector<PP> tv;
    for(auto p:v){
        bitset<M> have{};
        for(auto [t,_]:p.f){
            (void)_;
            for(size_t i = 0 ; i < n*m ; i ++){
                if(t[i])have[i]=true;
            }
        }
        tv.push_back({p,have});
    }
    auto cmp = [](const PP& x,const PP& y){
        return x.second.count() < y.second.count();
    };
    sort(tv.begin(),tv.end(),cmp);
    list<PP> l(tv.begin(),tv.end());
    return l;
}
E gen_point(){
    while(1){
        auto x = F{get_rand((1<<N)-1)};
        try{
            auto y = E::gety(x);
            E p(move(x),move(y),false);
            if((p+p).inf!=true){
                return p;
            }
        }
        catch (int e){
        }
    }
}
unordered_map<bitset<2*N>,bitset<2*N>> allpoint;
unordered_map<bitset<2*N>,size_t> allpoint_id;

void init_allpoint(){
    for(unsigned long long i = 0 ; i < (1ll<<(n)) ; i ++){
        try{
            F x{i};
            F y=E::gety(x);
            allpoint[x.x]=y.x;
            allpoint_id[x.x]=allpoint.size()-1;
        }
        catch(int e){
        }
    }
    printf("Init fin size: %d\n",(int)allpoint.size());
}

template <size_t M>
array<GF2,M> solve(list<pair<Poly<M,GF2>,bitset<M>>> l){
    array<GF2,M> result{};
    bitset<M> solved{};
    size_t solvedcnt=0;
    auto addans = [&](size_t idx,GF2 val){
        result[idx]=val;
        if(solved[idx]==false)solvedcnt++;
        solved[idx]=true;
        solvedcnt++;
        for(auto it=l.begin();it!=l.end();){
            auto& x=*it;
            if(x.second[idx]==false){
                it++;
                continue;
            }
            x.first=x.first.partialEvalNoShift(idx,val);
            x.second[idx]=false;
            if(x.second.count()==0){
                it=l.erase(it);
            }
            else{
                it++;
            }
        }
    };
    auto it = l.begin();
    array<GF2,M> tmp{};
    while(it!=l.end()&&solvedcnt!=M){
        if(it->second.count()==1){
            size_t idx=0;
            for(size_t i = 0 ; i < M ; i ++){
                if(it->second[i]){
                    idx=i;
                    break;
                }
            }
            {
                if(it->first.eval(tmp)==GF2{false}){
                    addans(idx,GF2{false});
                    it=l.begin();
                    continue;
                }
            }
            {
                tmp[idx]=GF2{true};
                if(it->first.eval(tmp)==GF2{false}){
                    tmp[idx]=GF2{false};
                    assert(it->first.eval(GF2{true})==GF2{true});
                    addans(idx,GF2{true});
                    it=l.begin();
                    continue;
                }
                tmp[idx]=GF2{false};
            }
        }
        it++;
    }
    assert(solvedcnt==M);
    return result;
}

template <size_t M>
array<long long,(1<<n)> solvePoint(array<F,M> pxa,const E& point){
    E p[M]{};
    array<long long,(1<<n)> res{};
    for(size_t i = 0 ; i < M ; i ++){
        p[i].x=pxa[i];
    }
    for(unsigned long long X = 0 ; X < (1<<M) ; X ++){
        for(size_t i = 0 ; i < M ; i ++){
            if((1ull<<i)&X){
                p[i].y=F{(allpoint[p[i].x.x])^(p[i].x.x)};
                res[allpoint_id[p[i].x.x]]=-1;
            }
            else{
                p[i].y=F{(allpoint[p[i].x.x])};
                res[allpoint_id[p[i].x.x]]=1;
            }
        }
        E tmp=-point;
        for(size_t i = 0 ; i < M ; i ++){
            tmp+=p[i];
        }
        if(tmp.inf==true){
            return res;
        }
    }
    assert(false);
    return res;
}

int main (){
    init_allpoint();
    E p = gen_point();
    auto ans=get_rand(1e9);
    E q = p*ans;
    cout << "P: " << p << endl << "Px" << ans << "=Q: " << q << endl;
    const auto& f=semaev_GF2n<m,E>;
    unordered_map<bitset<2*N>,unordered_set<bitset<2*N>>> picked;
    int cnt=0;
    while(1){
        cnt++;
        const auto a=get_rand(1e9);
        const auto b=get_rand(1e9);
        auto r = p*a+q*b;
        if(picked[r.x.x].count(r.y.x))continue;
        picked[r.x.x].insert(r.y.x);
        auto s = f.partialEval(0,r.x);
        cout << "s:" << s << endl;
        auto w = weilDescent<n,(1<<(m-2))>(s);
        auto g = grobnerGF2(conv(move(w)));
        auto l = convlist(move(g));
        assert(l.front().second.count()!=0);
        if(l.front().second.count()>1){
            continue;
        }
        auto gf2res=solve(move(l));
        array<F,m-1> pxa;
        for(size_t i = 0 ; i < m-1 ; i ++){
            for(size_t j = 0 ; j < n ; j ++){
                pxa[i].x[j]=gf2res[i*n+j].x;
            }
        }
        array<long long,(1<<n)> scaler = solvePoint(pxa,r);
        (void) scaler;
        cout << "fin,try: " << cnt << endl;
        cnt=0;
        getchar();
    }
}

