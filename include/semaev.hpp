#pragma once
#include "poly.hpp"
#include <cstdio>
#include <iostream>
template<size_t N,class EC>
const Poly<N,typename EC::F> __attribute__((init_priority(65535))) semaev_GF2n = [](){
    using std::array;
    using F = typename EC::F;
    using Term = typename Poly<N,F>::Term;
    const size_t J = N/2-1;
    const size_t DL = 1<<(N-J-2);
    const size_t DR = 1<<(J+2-2);
    const size_t MS = DL+DR;
    array<array<Poly<N,F>,MS>,MS> m{};
    array<Poly<N,F>,DL+1> l{};
    array<Poly<N,F>,DR+1> r{};
    const auto& lf = semaev_GF2n<N-J,EC>;
    const auto& rf = semaev_GF2n<J+2,EC>;
    for(const auto& [t,s]:lf.f){
        Term now{0};
        for(size_t i = 1 ; i < N-J ; i ++){
            now[i-1]=t[i];
        }
        l[t[0]]+=Poly<N,F>{{now,s}};
    }
    for(const auto& [t,s]:rf.f){
        Term now{0};
        for(size_t i = 1 ; i < J+2 ; i ++){
            now[N-J-1+i-1]=t[i];
        }
        r[t[0]]+=Poly<N,F>{{now,s}};
    }
    for(size_t i = 0 ; i < DR ; i ++){
        for(size_t j = 0 ; j < DL+1 ; j ++){
            m[i][i+j]=l[DL-j];
        }
    }
    for(size_t i = 0 ; i < DL ; i ++){
        for(size_t j = 0 ; j < DR+1 ; j ++){
            m[i+DR][i+j]=r[DR-j];
        }
    }
    return determinant(m);
}();
template<class EC>
const Poly<2,typename EC::F> __attribute__((init_priority(65535))) semaev_GF2n<2,EC> = [](){
    using std::array;
    using A=array<int,2>;
    using F = typename EC::F;
    const auto& one = F::one;
    A x{1},y{0,1};
    return Poly<2,F>{{x,one},{y,one}};
}();
template<class EC>
const Poly<3,typename EC::F> __attribute__((init_priority(65535))) semaev_GF2n<3,EC> = [](){
    using std::array;
    using A=array<int,3>;
    using F = typename EC::F;
    const auto& one = F::one;
    const auto& a6 = EC::a6;
    A t0{2,2,0},t1{0,2,2},t2{2,0,2},t3{1,1,1},t4{0,0,0};
    return Poly<3,F>{{t0,one},{t1,one},{t2,one},{t3,one},{t4,a6}};
}();
