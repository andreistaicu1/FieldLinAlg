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
    T0.print();

    std::vector<std::tuple<int, Polynomial>> output;

    step2:
    std::cout << "step2" << std::endl;
    T0.print();
    if(T0.degree == 0) {
        return output;
    }

    Polynomial T = Polynomial::gcd(T0, T0.derivative()).pure_mod(p);
    T.print();
    Polynomial V = std::get<0>(Polynomial::division_mod_p(T0,T, p));
    V.print();
    int k = 0;

    step3:
    std::cout << "step3" << std::endl;
    V.print();
    if(V.degree == 0) {
        std::cout << "V.degree == 0" << std::endl;
        T.print();
        T0 = divide_exponents(T, p);
        T0.print();
        e = p*e;
        goto step2;
    }
    
    k += 1;
    std::cout << "k: " << k << std::endl;
    if(k % p == 0) {
        std::cout << "k div by p" << std::endl;
        T = std::get<0>(Polynomial::division_mod_p(T, V, p));
        k += 1;
    }

    Polynomial W = Polynomial::gcd(T, V);
    W.print();
    Polynomial Aek = std::get<0>(Polynomial::division_mod_p(V, W, p));
    Aek.print();
    V = W;
    T = std::get<0>(Polynomial::division_mod_p(T, V, p));
    T.print();

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
    //Polynomial p1 = Polynomial({3,3,4,0,2,2,1,0,0});
    Polynomial p1 = Polynomial({4,2,3,1,1});
    Polynomial deriv = p1.derivative();
    p1.print();
    deriv.print();
    Polynomial gcdp = Polynomial::gcd(p1, deriv);
    gcdp.print();
    /*std::vector<std::tuple<int, Polynomial>> output = squarefree_factorization(p1, 5);
    std::cout << "factored" << std::endl;
    for(auto elem : output) {
        std::cout << std::get<0>(elem) << std::endl;;
        std::get<1>(elem).print();
    }*/
}