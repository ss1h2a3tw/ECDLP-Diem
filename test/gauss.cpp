#include<iostream>
#include<cassert>
#include"gauss.hpp"
using namespace std;

int main (){
    {
        gaussElimination<3> g;
        gaussElimination<3>::solutionType ret;
        ret = g.addPoly(gaussElimination<3>::polyType({8,7},array<long long,3>{1,2,3}));
        assert((ret == pair<long long,long long>{0,0}));
        ret = g.addPoly(gaussElimination<3>::polyType({18,15},array<long long,3>{3,3,3}));
        // -6 -6 = 0 -3 -6
        assert((ret == pair<long long,long long>{0,0}));
        ret = g.addPoly(gaussElimination<3>::polyType({4,3},array<long long,3>{2,0,2}));
        // -12 -11 = 0 -4 -4
        // -12 -9 = 0 0 12
        assert((ret == pair<long long,long long>{0,0}));
        ret = g.addPoly(gaussElimination<3>::polyType({2,1},array<long long,3>{2,2,3}));
        // -14 -13 = 0 -2 -3
        // -30 -27 = 0 0 3
        // -108 -99 = 0 0 0
        assert(( ret == pair<long long,long long>{-108,-99} ||
                    ret == pair<long long,long long>{108,99} ));
    }
    {
        gaussElimination<3> g;
        gaussElimination<3>::solutionType ret;
        ret = g.addPoly(gaussElimination<3>::polyType({8,7},array<long long,3>{0,0,0}));
        assert(( ret == pair<long long,long long>{-8,-7} ||
                    ret == pair<long long,long long>{8,7} ));
    }
    {
        gaussElimination<3> g;
        gaussElimination<3>::solutionType ret;
        ret = g.addPoly(gaussElimination<3>::polyType({3,2},array<long long,3>{1,0,0}));
        ret = g.addPoly(gaussElimination<3>::polyType({8,5},array<long long,3>{0,1,0}));
        ret = g.addPoly(gaussElimination<3>::polyType({8,2},array<long long,3>{0,0,1}));
        ret = g.addPoly(gaussElimination<3>::polyType({20,10},array<long long,3>{2,1,1}));
        assert(( ret == pair<long long,long long>{-2,-1} ||
                    ret == pair<long long,long long>{2,1} ));
    }
    {
        gaussElimination<3> g;
        gaussElimination<3>::solutionType ret;
        ret = g.addPoly(gaussElimination<3>::polyType({3,2},array<long long,3>{1,0,0}));
        ret = g.addPoly(gaussElimination<3>::polyType({8,5},array<long long,3>{0,1,0}));
        ret = g.addPoly(gaussElimination<3>::polyType({8,2},array<long long,3>{0,0,1}));
        ret = g.addPoly(gaussElimination<3>::polyType({20,10},array<long long,3>{2,0,1}));
        assert(( ret == pair<long long,long long>{-6,-4} ||
                    ret == pair<long long,long long>{6,4} ));
    }
}
