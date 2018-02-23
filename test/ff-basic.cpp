#include <cstdio>
#include <cassert>
#include <bitset>
#include <iostream>
#include "finitefield.hpp"
using namespace std;

const char irr[] = "10001101111";
int main (){
    //x^10 + x^6 + x^5 + x^3 + x^2 + x + 1
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
    assert(x*y==m);
    assert(x/y==d);
    assert(x/o==x);
    assert(y/o==y);
    cout << "Finite Field Basic Passed" << endl;
}
