#include <cstdio>
#include <cassert>
#include <bitset>
#include <iostream>
#include "finitefield.hpp"
using namespace std;

int main (){
    const bitset<11> irr("10001101111");
    GF2n<10> x(bitset<10>("0111100010"),irr);
    GF2n<10> y(bitset<10>("1001110011"),irr);
    GF2n<10> p(bitset<10>("1110010001"),irr);
    GF2n<10> m(bitset<10>("1100010011"),irr);
    assert(x+y==p);
    assert(x*y==m);
    cout << "Finite Field Basic Passed" << endl;
}
