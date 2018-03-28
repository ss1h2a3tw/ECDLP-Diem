#include<iostream>
#include<cassert>
#include"gauss.hpp"
using namespace std;

int main (){
    gaussElimination<3> g;
    gaussElimination<3>::solutionType ret;
    ret = g.addPoly(gaussElimination<3>::polyType({8,7},array<long long,3>{1,2,3}));
    assert((ret == pair<long long,long long>{0,0}));
    ret = g.addPoly(gaussElimination<3>::polyType({18,15},array<long long,3>{3,3,3}));
    assert((ret == pair<long long,long long>{0,0}));
    ret = g.addPoly(gaussElimination<3>::polyType({4,3},array<long long,3>{2,0,2}));
    assert((ret == pair<long long,long long>{0,0}));
    ret = g.addPoly(gaussElimination<3>::polyType({2,1},array<long long,3>{2,2,3}));
    assert((ret == pair<long long,long long>{12,11}));
}
