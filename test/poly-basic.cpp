#include <cassert>
#include <cstdio>
#include "poly.hpp"
#include "semaev.hpp"

using namespace std;
const char irr[] = "10001101111";
using F = GF2n<10,irr>;
const F a2("1000111000");
const F a6("0100001000");
using E = EC<F,a2,a6>;
int main (){
    using F = GF2n<10,irr>;
    //x^8 + x^7 + x^6 + x^5 + x
    F x("0111100010");
    //x^9 + x^6 + x^5 + x^4 + x + 1
    F y("1001110011");
    //x^9 + x^8 + x^7 + x^4 + 1
    F p("1110010001");
    //x^9 + x^8 + x^4 + x + 1
    F m("1100010011");
    using A = array<int,3>;
    using P = Poly<3,F>;
    P px{{A{0,0,0},x}};
    P py{{A{0,0,0},y},{A{1,0,0},y}};
    P pp{{A{0,0,0},p},{A{1,0,0},y}};
    P pm{{A{0,0,0},m},{A{1,0,0},m}};
    P ppp = px;
    ppp+=py;
    P pmp = px;
    pmp*=py;
    assert(px+py==pp);
    assert(ppp==pp);
    assert(px*py==pm);
    assert(pmp==pm);
}
