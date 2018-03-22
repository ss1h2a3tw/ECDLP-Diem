#pragma once
#include <array>

template<typename F>
F determinant(const std::array<std::array<F,1>,1>& m){
    return m[0][0];
}
template<typename F,size_t M>
F determinant(const std::array<std::array<F,M>,M>& m){
    using std::array;
    array<array<F,M-1>,M-1> tmp{};
    F ret{};
    for(size_t i = 0 ; i < M ; i ++){
        for(size_t r = 0 ; r < M-1 ; r ++){
            size_t nowc = 0;
            for(size_t c = 0 ; c < M-1 ; c ++){
                if(nowc==i)nowc++;
                tmp[r][c]=m[r+1][nowc];
                nowc++;
            }
        }
        if(i&1){
            ret-=m[0][i]*determinant(tmp);
        }
        else{
            ret+=m[0][i]*determinant(tmp);
        }
    }
    return ret;
}
