#include <cstdio>
#include <cassert>
#include <bitset>
#include <iostream>
#include "finitefield.hpp"
using namespace std;

//x^10 + x^6 + x^5 + x^3 + x^2 + x + 1
const char irr[] = "10001101111";
void testgf(){
    //x^8 + x^7 + x^6 + x^5 + x
    GF2n<10,irr> x(bitset<10>("0111100010"));
    //x^9 + x^6 + x^5 + x^4 + x + 1
    GF2n<10,irr> y(bitset<10>("1001110011"));
    //lets try string
    //x^9 + x^8 + x^7 + x^4 + 1
    GF2n<10,irr> p("1110010001");
    //x^9 + x^8 + x^4 + x + 1
    GF2n<10,irr> m("1100010011");
    //x^9 + x^8 + x^6 + x^5 + x^4 + x^2 + x + 1
    GF2n<10,irr> d("1101110111");
    GF2n<10,irr> o("0000000001");
    assert(x+y==p);
    assert(x-y==p);
    assert(x==-x);
    assert(x*y==m);
    assert(x/y==d);
    assert(x/o==x);
    assert(y/o==y);
    auto tmp = x;
    tmp+=y;
    assert(tmp==p);
    tmp = x;
    tmp-=y;
    assert(tmp==p);
    tmp=x;
    tmp*=y;
    assert(tmp==m);
    tmp=x;
    tmp/=y;
    assert(tmp==d);
    assert(x.pow(5)==x*x*x*x*x);
    cout << "Finite Field Basic Passed" << endl;
}
//x^9 + x^5 + x^4 + x^3
const GF2n<10,irr> a2("1000111000");
//x^8 + x^3
const GF2n<10,irr> a6("0100001000");
void testec(){
    using F = GF2n<10,irr>;
    using E = EC<F,a2,a6>;
    E inf(F::zero,F::zero,true);
    //x^8 + x^7 + x^6 + 1 : x^9 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x
    E x(F("0111000001"),F("1011111110"),false);
    // -x  = x^8 + x^7 + x^6 + 1 : x^9 + x^8 + x^5 + x^4 + x^3 + x^2 + x + 1
    E nx(F("0111000001"),F("1100111111"),false);
    // 2*x = x^8 + x^7 + x^6 + x^2 + x + 1 : x^8 + x^7 + x^2
    E dx(F("0111000111"),F("0110000100"),false);
    assert(inf==inf);
    assert(inf==-inf);
    assert(inf==inf+inf);
    assert(x+inf==x);
    assert(x==x);
    assert(x==-nx);
    assert(-x==nx);
    assert(x+nx==inf);
    assert(x+x==dx);
    assert(x*2==dx);
    assert(x*0==inf);
    assert(x*1==x);
    assert(x*4==dx*2);
    //order = 1002
    assert(x*1003==x);
    assert(x*(1002ull*999999999999ull+1)==x);
    assert(dx-x==x);
    //x^9 + x^8 + x^7 + x^5 + x^4 + x^3 + 1 : x^8 + x^5 + x^3 + x + 1
    E y(F("1110111001"),F("0100101011"),false);
    //x+y = x^8 + x^3 + x^2 : x^9 + x^7 + x^5
    E a(F("0100001100"),F("1010100000"),false);
    //x-y = (x^9 + x^8 + x^7 + x^5 + x^2 + 1 : x^7 + x^6 + x^3 + x^2 + 1
    E s(F("1110100101"),F("0011001101"),false);
    assert(x+y==a);
    assert(x-y==s);
    cout << "Elliptic Curve Basic Passed" << endl;
}
int main (){
    testgf();
    testec();
}
