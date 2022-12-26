#include "min_polynomial.hpp"

#include "doctest.h"
#include "polynomialmod.hpp"

PolynomialMod divide_exponents(PolynomialMod A){
    std::vector<int> output(std::max(A.degree,1), 0);
    int deg = std::div(A.degree, A.p).quot;

    for (int i = 0; i < deg+1; i++){
        output[i] = A.poly.p[i * A.p];
    }
    return PolynomialMod(output, A.p);
}

std::vector<std::tuple<int, PolynomialMod>> squarefree_factorization(PolynomialMod A) {
    int p = A.p;
    int e = 1;

    PolynomialMod T0 = A/A.lead_coeff;

    std::vector<std::tuple<int, PolynomialMod>> output;

    step2:
    if(T0.degree == 0) {
        return output;
    }

    PolynomialMod T = PolynomialMod::gcd(T0, T0.derivative());
    PolynomialMod V = T0/T;
    int k = 0;

    step3:
    if(V.degree == 0) {
        T0 = divide_exponents(T);
        e = p*e;
        goto step2;
    }
    
    k += 1;
    if(k % p == 0) {
        T = T/V;
        k += 1;
    }

    PolynomialMod W = PolynomialMod::gcd(T, V);
    PolynomialMod Aek = V/W;
    V = W;
    T = T/V;

    if(Aek.degree != 0) {
        output.push_back(std::make_tuple(e*k, Aek));
    }
    goto step3;
}


/*
std::vector<std::tuple<int, PolynomialMod>> distinct_degree_factorization(PolynomialMod A) {
    int p = A.p;

    std::vector<std::tuple<int, PolynomialMod>> output;

    PolynomialMod V = A;
    PolynomialMod W = PolynomialMod({0,1}, p);
    int d = 0;

    int e = V.degree;
    if(d+1 > e/2) {
        if(e > 0) {
            output.push_back(std::make_tuple(e, V));
        }
    }
    else {
        d += 1;
    }
}*/

TEST_CASE("divide_exponents") {
    int p = 5;
    PolynomialMod p1 = PolynomialMod({0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}, p);
    PolynomialMod divided_p = PolynomialMod({0,5,10,15}, p);
    CHECK(divide_exponents(p1) == divided_p);
}

TEST_CASE("squarefree_factorization") {
    int p = 5;
    PolynomialMod p1 = PolynomialMod({3,3,4,0,2,2,1,0,0}, p);

    PolynomialMod p2 = PolynomialMod({3,2,4}, p);
    auto p2_tup = std::make_tuple(1,p2);
    PolynomialMod p3 = PolynomialMod({1,1,3}, p);
    auto p3_tup = std::make_tuple(2,p3);
    auto predicted_output = std::vector({p2_tup, p3_tup});

    std::vector<std::tuple<int, PolynomialMod>> output = squarefree_factorization(p1);

    CHECK(output == predicted_output);
}