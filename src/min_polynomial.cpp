#include "min_polynomial.hpp"

#include "doctest.h"

Polynomial divide_exponents(Polynomial A, int p){
    std::vector<int> output(std::max(A.degree,1), 0);
    int deg = std::div(A.degree, p).quot;

    for (int i = 0; i < deg+1; i++){
        output[i] = A.p[i * p];
    }
    return Polynomial(output);
}

std::vector<std::tuple<int, Polynomial>> squarefree_factorization(Polynomial A, int p) {
    int e = 1;

    Polynomial T0 = (A*inverse_mod_p(A.lead_coeff, p)).pure_mod(p);

    std::vector<std::tuple<int, Polynomial>> output;

    step2:
    if(T0.degree == 0) {
        return output;
    }

    Polynomial T = Polynomial::gcd_mod_p(T0, T0.derivative(), p).pure_mod(p);
    Polynomial V = std::get<0>(Polynomial::division_mod_p(T0,T, p));
    int k = 0;

    step3:
    if(V.degree == 0) {
        T0 = divide_exponents(T, p);
        e = p*e;
        goto step2;
    }
    
    k += 1;
    if(k % p == 0) {
        T = std::get<0>(Polynomial::division_mod_p(T, V, p));
        k += 1;
    }

    Polynomial W = Polynomial::gcd_mod_p(T, V, p);
    Polynomial Aek = std::get<0>(Polynomial::division_mod_p(V, W, p));
    V = W;
    T = std::get<0>(Polynomial::division_mod_p(T, V, p));

    if(Aek.degree != 0) {
        output.push_back(std::make_tuple(e*k, Aek));
    }
    goto step3;
}

TEST_CASE("divide_exponents") {
    Polynomial p1 = Polynomial({0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15});
    Polynomial divided_p = Polynomial({0,5,10,15});
    CHECK(divide_exponents(p1, 5) == divided_p);
}
TEST_CASE("squarefree_factorization") {
    Polynomial p1 = Polynomial({3,3,4,0,2,2,1,0,0});

    Polynomial p2 = Polynomial({3,2,4});
    auto p2_tup = std::make_tuple(1,p2);
    Polynomial p3 = Polynomial({1,1,3});
    auto p3_tup = std::make_tuple(2,p3);
    auto predicted_output = std::vector({p2_tup, p3_tup});

    std::vector<std::tuple<int, Polynomial>> output = squarefree_factorization(p1, 5);

    CHECK(output == predicted_output);
}