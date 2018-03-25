#include "finitefield.hpp"

GF2 GF2::operator-()const{
    return *this;
}
GF2 GF2::operator-(const GF2& r)const{
    return GF2{x^r.x};
}
GF2& GF2::operator-=(const GF2& r){
    x^=r.x;
    return *this;
}
GF2 GF2::operator+(const GF2& r)const{
    return GF2{x^r.x};
}
GF2& GF2::operator+=(const GF2&r ){
    x^=r.x;
    return *this;
}
GF2 GF2::operator*(const GF2& r)const{
    return GF2{x&r.x};
}
GF2& GF2::operator*=(const GF2&r ){
    this->x&=r.x;
    return *this;
}
bool GF2::operator==(const GF2& r)const{
    return x==r.x;
}
bool GF2::operator!=(const GF2& r)const{
    return x!=r.x;
}
bool GF2::iszero()const{
    return !x;
}
std::ostream& operator<<(std::ostream& os,const GF2 x){
    os << (int)x.x;
    return os;
}
