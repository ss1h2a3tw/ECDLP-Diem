#include <cstdio>
#include <cassert>
#include <bitset>
#include <iostream>
#include "finitefield.hpp"
using namespace std;

int main (){
    //x^10 + x^6 + x^5 + x^3 + x^2 + x + 1
    const bitset<11> irr("10001101111");
    //x^8 + x^7 + x^6 + x^5 + x
    GF2n<10> x(bitset<10>("0111100010"),irr);
    //x^9 + x^6 + x^5 + x^4 + x + 1
    GF2n<10> y(bitset<10>("1001110011"),irr);
    //lets try string
    //x^9 + x^8 + x^7 + x^4 + 1
    GF2n<10> p("1110010001","10001101111");
    //x^9 + x^8 + x^4 + x + 1
    GF2n<10> m("1100010011","10001101111");
    //x^9 + x^8 + x^6 + x^5 + x^4 + x^2 + x + 1
    GF2n<10> d("1101110111","10001101111");
    assert(x+y==p);
    assert(x*y==m);
    assert(x/y==d);
    cout << "Finite Field Basic Passed" << endl;
}
