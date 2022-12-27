#include "min_polynomial.hpp"

#include "doctest.h"
#include "polynomialmod.hpp"

#include <random>

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


std::vector<std::tuple<int, PolynomialMod>> distinct_degree_factorization(PolynomialMod A) {
    int p = A.p;

    std::vector<std::tuple<int, PolynomialMod>> output;

    PolynomialMod V = A;
    PolynomialMod X = PolynomialMod({0,1}, p);
    PolynomialMod W = X;
    int d = 0;

    while(true) {
        int e = V.degree;
        std::cout << "d: " << d << std::endl;
        std::cout << "e: " << e << std::endl;
        if(d+1 > e/2) {
            if(e > 0) {
                output.push_back(std::make_tuple(e, V/V.lead_coeff));
            }

            for (auto [d, B] : output){
                std::cout << "degree: " << d << std::endl;
                B.print();
            }

            return output;
        }
        else {
            d += 1;
            W = W.exp(p).mod(V);
        }

        std::cout << "W: " << std::endl;
        W.print();

        PolynomialMod Ad = PolynomialMod::gcd(W - X, V);

        std::cout << "Ad: " << std::endl;
        Ad.print();

        output.push_back(std::make_tuple(d, Ad/Ad.lead_coeff));

        if(!(Ad == PolynomialMod({1,0},p))) {
            V = V/Ad;
            W = W.mod(V);
        }

        std::cout << "V: " << std::endl;
        V.print();

        if(W == PolynomialMod({0}, p)) {
            //throw;
        }
    }
}

unsigned SEED = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine GENERATOR(SEED);

PolynomialMod random_polynomial(int degree, int p) {
    std::vector<int> output(degree+1, 0);
    output[degree] = 1;

    std::uniform_int_distribution<int> distribution(0,p);
    for(auto i = 0; i < degree; ++i) {
        output[i] = distribution(GENERATOR);
    }

    return PolynomialMod(output, p);
}

std::vector<PolynomialMod> cantor_zassenhaus_split_2(PolynomialMod A, int d) {
    int p = A.p;

    int k = A.degree / d;

    if(k == 1) {
        return {A};
    }
    PolynomialMod T = PolynomialMod({0,1}, p);

    PolynomialMod B = PolynomialMod(p);
    while(true) {
        PolynomialMod C = T;
        for(auto i = 0; i < d-1; ++i) {
            C = (T+C.exp(2)).mod(A);
        }

        B = PolynomialMod::gcd(A, C);
        if(!(B.degree == 0 || B.degree == A.degree)) {
            break;
        }
        T = T*PolynomialMod({0,0,1}, p);
    }

    auto output1 = cantor_zassenhaus_split_2(A/B, d);
    auto output2 = cantor_zassenhaus_split_2(B, d);
    output1.insert(output1.end(), output2.begin(), output2.end());  // Optimize by moving in place

    return output1;
}

std::vector<PolynomialMod> cantor_zassenhaus_split(PolynomialMod input_poly, int d) {
    int p = input_poly.p;

    if(p == 2) {
        return cantor_zassenhaus_split_2(input_poly, d);
    }
    
    PolynomialMod A = input_poly;

    int k = A.degree / d;
    if (k == 1){
        return {A/A.lead_coeff};
    }

    PolynomialMod B = PolynomialMod(p);
    while(true){
        int exponent = std::div(std::pow(p, d) - 1, 2).quot;
        PolynomialMod T = random_polynomial((2 * d) - 1, p).exp(exponent);
        B = PolynomialMod::gcd(A, T - 1);
        if (!(B.degree == 0 || B.degree == A.degree)){
            break;
        }
    }

    auto output1 = cantor_zassenhaus_split(A / B, d);
    auto output2 = cantor_zassenhaus_split(B, d);
    output1.insert(output1.end(), output2.begin(), output2.end());  // Optimize by moving in place

    return output1;
}

PolynomialMod create_min_polynomial(int p, int degree) {
    int leading_index = std::pow(p, degree);
    std::vector<int> initial_vec(leading_index + 1, 0);

    initial_vec[leading_index] = 1;
    initial_vec[1] = 1;

    PolynomialMod initial = PolynomialMod(initial_vec, p);

    std::vector<std::tuple<int, PolynomialMod>> squarefree = squarefree_factorization(initial);
    std::cout << "Length: " << squarefree.size() << std::endl;

    for(auto [d, A] : squarefree) {
        if(A.degree < degree) {
            continue;
        }
        std::vector<std::tuple<int, PolynomialMod>> dd_factor = distinct_degree_factorization(A); // INFINITE LOOP

        for(auto [deg, B] : dd_factor) {
            if (deg == degree) {
                return cantor_zassenhaus_split(B, deg)[0]; // Wish us luck
            }
        }
    }
}

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

TEST_CASE("distinct degree factorization") {
    /*
    int p = 5;
    PolynomialMod p1 = PolynomialMod({1,1,2,1,1}, p);
    
    auto output = distinct_degree_factorization(p1);
    for(auto [i, Ai] : output) {
        std::cout << i << std::endl;
        Ai.print();
    }*/
    int p = 5;
    PolynomialMod A = PolynomialMod({0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, p);
    auto output = distinct_degree_factorization(A);
    std::cout << "we made it into the test case!" << std::endl;
    std::cout << "and look it has length: " << output.size() << std::endl;
    for(auto [i, Ai] : output) {
        std::cout << i << std::endl;
        Ai.print();
    }
}

TEST_CASE("cantor_zassenhaus_split") {
    SUBCASE("p != 2") {
        int p = 5;
        PolynomialMod p1 = PolynomialMod({1,0,1}, p);

        auto output = cantor_zassenhaus_split(p1, 1);
        std::cout << "length: " << output.size() << std::endl;
        for(PolynomialMod elem : output) {
            //elem.print();
        }
    }
    SUBCASE("p == 2") {
        int p = 2;
        PolynomialMod p1 = PolynomialMod({0,1,1}, p);

        auto output = cantor_zassenhaus_split(p1, 1);
        std::cout << "length: " << output.size() << std::endl;
        for(PolynomialMod elem : output) {
            //elem.print();
        }
    }
}

TEST_CASE("create_min_polynomial") {
    auto min_polynomial = create_min_polynomial(5, 2);
    //min_polynomial.print();
}