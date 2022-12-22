#include <vector>
#include <memory>
#include <cmath>
#include <numeric>

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <cstdlib>

#include <iostream>
#include <tuple>
#include <unordered_map>

#include "min_polynomial.cpp"
#include "polynomial.hpp"

#include "doctest.h"

class ff_element;

const Polynomial ZERO_P = Polynomial(Polynomial::ZERO);
const Polynomial ONE_P = Polynomial(Polynomial::ONE);

class finite_field {
public:
    int p_char, degree;
    Polynomial min_polynomial;

    finite_field() {
        this->p_char = 0;
        this->degree = 0;
        this->min_polynomial = Polynomial();
    }
    
    finite_field(int p_char, int degree) : p_char(p_char), degree(degree) {
        this->min_polynomial = Polynomial::zero_polynomial(1);
    };

    bool operator== (finite_field const &other) const {
        return this->p_char == other.p_char && this->degree == other.degree;
    }

    int inverse_mod_p(int val) const {
        if (val == 1){
            return val;
        }
        int inverse = euclidean_algorithm(p_char, val, std::make_pair(1, 0), std::make_pair(0, 1));
        return (((inverse % p_char) + p_char) % p_char);
    }

    int euclidean_algorithm(int a, int b, std::tuple<int, int> compute_a, std::tuple<int, int> compute_b) const {
        std::div_t div_n = std::div(a, b);
        int q = div_n.quot;
        int r = div_n.rem;

        int compute_r_0 = std::get<0>(compute_a) - (std::get<0>(compute_b) * q);
        int compute_r_1 = std::get<1>(compute_a) - (std::get<1>(compute_b) * q);

        if (r == 1) {
            return compute_r_1;
        }

        auto compute_r = std::make_pair(compute_r_0, compute_r_1);
        return euclidean_algorithm(b, r, compute_b, compute_r);
    }
    
};

class ff_element {
public:
    Polynomial polynomial;
    finite_field ff;
    
    ff_element(Polynomial polynomial, finite_field ff) {
        assert(polynomial.degree <= ff.degree);
        polynomial.mod(ff.p_char);
        this->polynomial = polynomial.trim_end();
        this->ff = ff;
    };

    ff_element(finite_field ff) : polynomial(Polynomial::zero_polynomial(ff.degree)), ff(ff) {};

    bool operator==(ff_element const &other) const {
        return polynomial==other.polynomial && ff==other.ff;
    }

    ff_element operator= (ff_element const &other) {
        assert(this->ff == other.ff);

        this->polynomial = other.polynomial;
        return (*this);
    }

    ff_element operator+ (ff_element const &other) const {
        assert(this->ff == other.ff);

        Polynomial output = this->polynomial + other.polynomial;
        output.mod(this->ff.p_char);

        return ff_element(output, ff);
    };

    ff_element operator- (ff_element const &other) const {
        assert(this->ff == other.ff);

        Polynomial output = this->polynomial - other.polynomial;
        output.mod(this->ff.p_char);

        return ff_element(output, this->ff);
    };

    ff_element operator* (ff_element const &other) const {
        assert(this->ff == other.ff);

        std::tuple<Polynomial, Polynomial, int> triple = Polynomial::pseudo_div(polynomial*other.polynomial, ff.min_polynomial);

        Polynomial remainder = std::get<1>(triple);
        remainder.mod(this->ff.p_char);

        int c = std::get<2>(triple);
        int c_inverse = ff.inverse_mod_p(c % ff.p_char);

        Polynomial output = c_inverse*remainder;
        output.mod(this->ff.p_char);

        return ff_element(output, this->ff);
    };

    ff_element operator/ (ff_element const &other) const {
        assert(this->ff == other.ff);
        assert(!(ZERO_P == other.polynomial));

        ff_element inverse = other.make_inverse();

        return (*this) * inverse;
    };

    ff_element make_inverse() const {
        if (this->polynomial == ONE_P){
            return *this;
        }
        return ff_element(euclidean_algorithm(ff.min_polynomial, polynomial), ff);
    }
    
    Polynomial euclidean_algorithm(Polynomial elemA, Polynomial elemB) const {
        Polynomial A = elemA;
        Polynomial B = elemB;

        if(B.degree > A.degree) {
            return euclidean_algorithm(elemB, elemA);
        }

        std::tuple<Polynomial, Polynomial> compute_A = std::make_pair(Polynomial::ONE, Polynomial::ZERO);
        std::tuple<Polynomial, Polynomial> compute_B = std::make_pair(Polynomial::ZERO, Polynomial::ONE);

        assert(!(B == Polynomial::ZERO));

        while(true) {
            auto division = Polynomial::pseudo_div(A, B);

            int delta = std::get<2>(division);
            int delta_inverse = ff.inverse_mod_p(((delta % ff.p_char)+ ff.p_char) % ff.p_char);

            Polynomial Q = delta_inverse*std::get<0>(division);
            Polynomial R = delta_inverse*std::get<1>(division);
            Q.mod(ff.p_char);
            R.mod(ff.p_char);

            Polynomial compute_r_0 = std::get<0>(compute_A) - std::get<0>(compute_B)*Q;
            Polynomial compute_r_1 = std::get<1>(compute_A) - std::get<1>(compute_B)*Q;
            compute_r_0.mod(ff.p_char);
            compute_r_1.mod(ff.p_char);

            if(R.degree == 0) {
                return ff.inverse_mod_p(R.lead_coeff)*compute_r_1;
            }

            compute_B = std::make_pair(compute_r_0, compute_r_1);

            A = B;
            B = R;
        }
    }
};

TEST_CASE("finite field") {
    SUBCASE("p_char = 5, degree = 1") {
        finite_field ff = finite_field(5, 1);

        SUBCASE("inverse_mod_p") {
            CHECK(ff.inverse_mod_p(1) == 1);
            CHECK(ff.inverse_mod_p(2) == 3);
            CHECK(ff.inverse_mod_p(4) == 4);
            CHECK(ff.inverse_mod_p(39) == 4);
            CHECK(ff.inverse_mod_p(113) == 2);
        }

        SUBCASE("euclidean_algorithm") {

        }
    }
}

TEST_CASE("ff_element") {
    finite_field ff = finite_field(5, 3);
    ff.min_polynomial = Polynomial({1,2,0,1});

    Polynomial p1 = Polynomial({1,2,3});
    ff_element elem1 = ff_element(p1, ff);

    Polynomial p2 = Polynomial({2,4,6});
    ff_element elem2 = ff_element(p2, ff);

    Polynomial p3 = Polynomial({2,1,2});
    ff_element elem3 = ff_element(p3, ff);

    ff_element zero_e = ff_element(ZERO_P, ff);

    SUBCASE("adding ff_element") {
        CHECK((elem1 + elem1) == elem2);
    }

    SUBCASE("subtracting ff_element") {
        CHECK((elem1 - elem1) == zero_e);
    }

    SUBCASE("multiplication") {
        Polynomial p4 = Polynomial({0,0,3});
        ff_element elem4 = ff_element(p4, ff);

        CHECK(elem1*elem3 == elem4);
    }

    SUBCASE("division") {
        Polynomial p4 = Polynomial({2,2,4});
        ff_element elem4 = ff_element(p4, ff);

        CHECK(elem1/elem3 == elem4);
    }
}
