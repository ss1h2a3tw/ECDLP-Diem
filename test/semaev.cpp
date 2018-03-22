#include <cassert>
#include <iostream>
#include "semaev.hpp"
#include "determinant.hpp"

using namespace std;
//x^10 + x^6 + x^5 + x^3 + x^2 + x + 1
const char irr[] = "10001101111";
//x^9 + x^5 + x^4 + x^3
const GF2n<10,irr> a2("1000111000");
//x^8 + x^3
const GF2n<10,irr> a6("0100001000");
int main (){
    using F = GF2n<10,irr>;
    using E = EC<F,a2,a6>;
    E inf(F::zero,F::zero,true);
    //x^8 + x^7 + x^6 + 1 : x^9 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x
    E x(F("0111000001"),F("1011111110"),false);
    E y = x*2;
    E z = -(x+y);
    E a = x*3;
    E b = -(x+y+a);
    const auto& f2 = semaev_GF2n<2,E>;
    const auto& f3 = semaev_GF2n<3,E>;
    const auto& f4 = semaev_GF2n<4,E>;
    cout << "a6^2" << a6.pow(2) << endl;
    assert(f2.eval({x.x,(-x).x})==F{});
    assert(f3.eval({x.x,y.x,z.x})==F{});
    assert(a.inf==false);
    assert(b.inf==false);
    assert((x+y+a+b).inf);
    assert(f4.eval({x.x,y.x,a.x,b.x})==F{});
}
