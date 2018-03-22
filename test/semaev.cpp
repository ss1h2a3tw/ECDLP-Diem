#include <cassert>
#include <iostream>
#include "finitefield.hpp"
#include "semaev.hpp"
#include "determinant.hpp"

using namespace std;
//x^10 + x^6 + x^5 + x^3 + x^2 + x + 1
const char irr[] = "10001101111";
//x^9 + x^5 + x^4 + x^3
const GF2n<10,irr> __attribute__((init_priority(200))) a2("1000111000");
//x^8 + x^3
const GF2n<10,irr> __attribute__((init_priority(200))) a6("0100001000");
int main (){
    using F = GF2n<10,irr>;
    using E = EC<F,a2,a6>;
    E inf(F::zero,F::zero,true);
    //x^8 + x^7 + x^6 + 1 : x^9 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x
    E x(F("0111000001"),F("1011111110"),false);
    E y = x*2;
    E z = x*3;
    E a = x*4;
    E b = x*5;
    E c = x*6;
    assert(y.inf==false);
    assert(z.inf==false);
    assert(a.inf==false);
    assert(b.inf==false);
    assert(c.inf==false);
    assert((x+y+a+b).inf==false);
    const auto& f2 = semaev_GF2n<2,E>;
    const auto& f3 = semaev_GF2n<3,E>;
    const auto& f4 = semaev_GF2n<4,E>;
    const auto& f5 = semaev_GF2n<5,E>;
    //const auto& f6 = semaev_GF2n<6,E>;
    assert(f2.eval({x.x,(-x).x})==F{});
    assert(f3.eval({x.x,y.x,(x+y).x})==F{});
    assert(f4.eval({x.x,y.x,a.x,(x+y+a).x})==F{});
    assert(f5.eval({x.x,y.x,a.x,b.x,(-(x+y+a+b)).x})==F{});
    //assert(f6.eval({x.x,y.x,a.x,b.x,c.x,(-(x+y+a+b+c)).x})==F{});
}
