#include <iostream>
#include <vector>
#include <cassert>
#include "finitefield.hpp"
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

void basicPlus(){
    //x^3 + x + 1
    static const char irr [] = "1011";
    using F = GF2n<3,irr>;
    using A = array<int,2>;
    Poly<2,F> p{{A{0,1},F::one},{A{1,0},F::one}};
    auto v=conv(weilDescent<3,2>(p));
    auto gb=grobnerGF2(v);
    (void)gb;
}
void basicMul(){
    //x^3 + x + 1
    static const char irr [] = "1011";
    using F = GF2n<3,irr>;
    using A = array<int,2>;
    Poly<2,F> p{{A{1,1},F::one}};
    auto v = conv(weilDescent<3,2>(p));
    auto gb = grobnerGF2(v);
    (void)gb;
}
//x^10 + x^6 + x^5 + x^3 + x^2 + x + 1
static const char irr[] = "10001101111";
//x^9 + x^5 + x^4 + x^3
static const GF2n<10,irr> a2("1000111000");
//x^8 + x^3
static const GF2n<10,irr> a6("0100001000");
void semaev(){
    using F = GF2n<10,irr>;
    using E = EC<F,a2,a6>;
    const auto& f4 = semaev_GF2n<3,E>;
    auto p = weilDescent<2,2>(f4);
    auto v=conv(p);
    auto gb = grobnerGF2(v);
    (void)gb;
}
int main (){
    basicPlus();
    basicMul();
    printf("Basic Passed\n");
    semaev();
    printf("Semaev generated\n");
}
