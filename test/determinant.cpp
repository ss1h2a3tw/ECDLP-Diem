#include <cassert>
#include <array>
#include "determinant.hpp"

using namespace std;
int main (){
    using sA = array<int,2>;
    array<array<int,2>,2> sm{
        sA{1,2},
        sA{3,4}
    };
    assert(determinant(sm)==-2);
    using mA = array<int,3>;
    array<array<int,3>,3> mm{
        mA{1,2,3},
        mA{4,5,6},
        mA{7,8,9}
    };
    assert(determinant(mm)==0);

    using A = array<int,5>;
    array<array<int,5>,5> m{
        A{1,2,4,3,0},
        A{2,1,-1,1,3},
        A{4,-1,-2,5,1},
        A{7,3,6,2,1},
        A{1,0,-1,1,1}
    };
    assert(determinant(m)==-34);
}
