#include<array>
#include<vector>
#include<utility>
#include<algorithm>

template<size_t M>
class gaussElimination{
    public:
        using LL = long long;
        using solutionType = std::pair<LL, LL>;
        class polyType{
            public:
                polyType(const std::pair<LL, LL>& _lhs, std::array<LL,M>&& _rhs):lhs(_lhs), rhs(move(_rhs)){
                    firstItemPos = _firstItemPos();
                }
                polyType(const std::pair<LL, LL>& _lhs, const std::array<LL,M>& _rhs):lhs(_lhs), rhs(_rhs){
                    firstItemPos = _firstItemPos();
                }
                polyType() = default;
                polyType(polyType&&) = default;
                polyType(const polyType&) = default;
                polyType& operator=(polyType&&) = default;
                polyType& operator=(const polyType&) = default;
                int order(const polyType& p)const{
                    return firstItemPos - p.firstItemPos;
                }
                void doEli(const polyType& p){
                    if (firstItemPos != p.firstItemPos || firstItemPos == M) return;
                    LL l = lcm( rhs[firstItemPos], p.rhs[p.firstItemPos] );
                    LL sc = (-l)/rhs[firstItemPos];
                    LL pc = l/p.rhs[p.firstItemPos];
                    for (size_t i=0; i<M; i++){
                        rhs[i] *= sc;
                        rhs[i] += p.rhs[i] * pc;
                    }
                    lhs.first *= sc;
                    lhs.first += p.lhs.first * pc;
                    lhs.second *= sc;
                    lhs.second += p.lhs.second * pc;
                    firstItemPos = _firstItemPos();
                }
                size_t getFirstItemPos()const{ return firstItemPos; }
                solutionType getSolution()const{ return lhs; }
            private:
                LL lcm(LL a, LL b)const{
                    LL tmp = gcd(a, b);
                    a /= tmp;
                    return a * b;
                }
                LL gcd(LL a, LL b)const{
                    return b ? gcd(b, a%b) : a;
                }
                size_t _firstItemPos(){ return std::find_if(rhs.begin(), rhs.end(), [](LL x){return x!=0;}) - rhs.begin(); }
                std::pair<LL,LL> lhs;
                std::array<LL,M> rhs;
                size_t firstItemPos;
        };
        solutionType addPoly(const std::pair<LL,LL>& p, std::array<LL,M>&& a){
            return addPoly(polyType(p, move(a)));
        }
        solutionType addPoly(const std::pair<LL,LL>& p, const std::array<LL,M>& a){
            return addPoly(polyType(p, a));
        }
        solutionType addPoly(polyType p){
            for (size_t i=0; i<polys.size(); i++){
                int tmp = polys[i].order(p);
                if (tmp > 0){
                    polys.insert(polys.begin()+i, std::move(p));
                    return {0,0};
                }
                else if (!tmp){
                    p.doEli(polys[i]);
                }
            }
            polys.push_back(std::move(p));
            if (polys.back().getFirstItemPos() == M) return polys.back().getSolution();
            return {0,0};
        }
    private:
        std::vector<polyType> polys;
};
