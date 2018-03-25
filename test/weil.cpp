#include <iostream>
#include <cassert>
#include "finitefield.hpp"
#include "semaev.hpp"
#include "weil.hpp"

using namespace std;
//x^3 + x + 1
//x^9 + x^5 + x^4 + x^3
//const GF2n<10,irr> a2("1000111000");
//x^8 + x^3
//const GF2n<10,irr> a6("0100001000");
void basicPlus(){
    static const char irr [] = "1011";
    using F = GF2n<3,irr>;
    using bF = GF2;
    bF bo{true};
    using A = array<int,2>;
    using bA = array<int,6>;
    Poly<2,F> p{{A{0,1},F::one},{A{1,0},F::one}};
    array<Poly<6,bF>,3> res = {
        Poly<6,bF>{{bA{1,0,0,0,0,0},bo},{bA{0,0,0,1,0,0},bo}},
        Poly<6,bF>{{bA{0,1,0,0,0,0},bo},{bA{0,0,0,0,1,0},bo}},
        Poly<6,bF>{{bA{0,0,1,0,0,0},bo},{bA{0,0,0,0,0,1},bo}}};
    assert(weilDescent<3>(p)==res);
}
void basicMul(){
    static const char irr [] = "1011";
    using F = GF2n<3,irr>;
    using bF = GF2;
    bF bo{true};
    using A = array<int,2>;
    using bA = array<int,6>;
    Poly<2,F> p{{A{1,1},F::one}};
    array<Poly<6,bF>,3> res = {
        Poly<6,bF>{{bA{1,0,0,1,0,0},bo},{bA{0,1,0,0,0,1},bo},{bA{0,0,1,0,1,0},bo}},
        Poly<6,bF>{{bA{0,1,0,1,0,0},bo},{bA{1,0,0,0,1,0},bo},{bA{0,0,1,0,0,1},bo},{bA{0,1,0,0,0,1},bo},{bA{0,0,1,0,1,0},bo}},
        Poly<6,bF>{{bA{0,0,1,1,0,0},bo},{bA{0,1,0,0,1,0},bo},{bA{1,0,0,0,0,1},bo},{bA{0,0,1,0,0,1},bo}}};
    assert(weilDescent<3>(p)==res);
}
int main (){
    //using F = GF2n<10,irr>;
    //using E = EC<F,a2,a6>;
    //const auto& f5 = semaev_GF2n<5,E>;
    //printf("!!%d",(int)f5.f.size());
    //getchar();
    //getchar();
    //getchar();
    //auto p = weilDescent<5>(f5);
    //(void)p;
    basicPlus();
    basicMul();
}
