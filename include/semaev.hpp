#pragma once
#include "poly.hpp"
#include <cstdio>
#include <iostream>
template<size_t N,class EC>
const Poly<N,typename EC::F> semaev_GF2n = [](){
    using std::array;
    using F = typename EC::F;
    using Term = typename Poly<N,F>::Term;
    const size_t J = N/2-1;
    array<array<Poly<N,F>,N>,N> m{};
    array<Poly<N,F>,N-J> l{};
    array<Poly<N,F>,J+2> r{};
    const auto& lf = semaev_GF2n<N-J,EC>;
    const auto& rf = semaev_GF2n<J+2,EC>;
    for(const auto& [t,s]:lf.f){
        Term now{};
        for(size_t i = 1 ; i < N-J ; i ++){
            now[i-1]=t[i];
        }
        l[t[0]]+=Poly<N,F>{{now,s}};
    }
    printf("L (%d):\n",N-J);
    for(int i = 0 ; i < N-J ; i ++){
        printf("%d pow\n",i);
        for(auto [t,s]:l[i].f){
            for(int i = 0 ; i < N ; i ++){
                printf("%d-",t[i]);
            }
            printf("\n");
        }
    }
    for(const auto& [t,s]:rf.f){
        Term now{};
        for(size_t i = 1 ; i < J+2 ; i ++){
            now[N-J-1+i-1]=t[i];
        }
        r[t[0]]+=Poly<N,F>{{now,s}};
    }
    printf("R (%d):\n",J+2);
    for(int i = 0 ; i < J+2 ; i ++){
        printf("%d pow\n",i);
        for(auto [t,s]:r[i].f){
            for(int i = 0 ; i < N ; i ++){
                printf("%d-",t[i]);
            }
            printf("\n");
        }
    }
    for(size_t i = 0 ; i < J+1 ; i ++){
        for(size_t j = 0 ; j < N-J ; j ++){
            m[i][i+j]=l[N-J-1-j];
        }
    }
    for(size_t i = 0 ; i < N-J-1 ; i ++){
        for(size_t j = 0 ; j < J+2 ; j ++){
            m[i+J+1][i+j]=r[J+1-j];
        }
    }
    for(int i = 0 ; i < N ; i ++){
        for(int j = 0 ; j < N ; j ++){
            printf("(%d,%d):\n",i,j);
            for(auto [t,s]:m[i][j].f){
                std::cout << s << "* (";
                for(int k = 0 ; k < N ; k ++){
                    std::cout << t[k] << "-";
                }
                std::cout << ")\n";
            }
        }
    }
    auto ttt = determinant(m);
    printf("res\n");
    for(auto [t,s]:ttt.f){
        std::cout << s << "* (";
        for(int k = 0 ; k < N ; k ++){
            std::cout << t[k] << "-";
        }
        std::cout << ")\n";
    }
    return ttt;
}();
template<class EC>
const Poly<2,typename EC::F> semaev_GF2n<2,EC> = [](){
    using std::array;
    using A=array<int,2>;
    using F = typename EC::F;
    const auto& one = F::one;
    A x{1},y{0,1};
    return Poly<2,F>{{x,one},{y,one}};
}();
template<class EC>
const Poly<3,typename EC::F> semaev_GF2n<3,EC> = [](){
    using std::array;
    using A=array<int,3>;
    using F = typename EC::F;
    const auto& one = F::one;
    const auto& a6 = EC::a6;
    A t0{2,2,0},t1{0,2,2},t2{2,0,2},t3{1,1,1},t4{0,0,0};
    return Poly<3,F>{{t0,one},{t1,one},{t2,one},{t3,one},{t4,a6}};
}();
