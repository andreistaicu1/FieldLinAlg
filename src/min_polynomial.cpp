#include "min_polynomial.hpp"

#include "doctest.h"

Polynomial divide_exponents(Polynomial A, int p){
    std::vector<int> output(A.degree, 0);
    int deg = std::div(A.degree, p).quot;

    for (int i = 0; i < deg; i++){
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

    Polynomial T = Polynomial::gcd(T0, T0.derivative()).pure_mod(p);
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

    Polynomial W = Polynomial::gcd(T, V);
    Polynomial Aek = std::get<0>(Polynomial::division_mod_p(V, W, p));
    V = W;
    T = std::get<0>(Polynomial::division_mod_p(T, V, p));

    if(Aek.degree != 0) {
        output.push_back(std::make_tuple(e*k, Aek));
    }
    goto step3;
}