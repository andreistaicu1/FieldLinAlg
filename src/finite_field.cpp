#include "finite_field.hpp"
#include "doctest.h"

finite_field::finite_field() : p_char(0), degree(0), min_polynomial(PolynomialMod(0)) {};

finite_field::finite_field(int p_char, int degree) : p_char(p_char), degree(degree), min_polynomial(PolynomialMod(p_char)) {};

bool finite_field::operator== (finite_field const &other) const {
    return this->p_char == other.p_char && this->degree == other.degree;
}
    
ff_element::ff_element(PolynomialMod polynomial, finite_field ff) : polynomial(polynomial), ff(ff) {};

ff_element::ff_element(finite_field ff) : polynomial(PolynomialMod(ff.p_char)), ff(ff) {};

bool ff_element::operator==(ff_element const &other) const {
    return polynomial==other.polynomial && ff==other.ff;
}

ff_element ff_element::operator= (ff_element const &other) {
    assert(this->ff == other.ff);

    this->polynomial = other.polynomial;
    return (*this);
}

ff_element ff_element::operator+ (ff_element const &other) const {
    assert(this->ff == other.ff);

    return ff_element(polynomial + other.polynomial, ff);
};

ff_element ff_element::operator- (ff_element const &other) const {
    assert(this->ff == other.ff);

    return ff_element(polynomial - other.polynomial, this->ff);
};

ff_element ff_element::operator* (ff_element const &other) const {
    assert(this->ff == other.ff);

    PolynomialMod R = (polynomial*other.polynomial).mod(ff.min_polynomial);

    return ff_element(R, ff);
};

ff_element ff_element::operator/ (ff_element const &other) const {
    assert(this->ff == other.ff);
    assert(!(PolynomialMod({0}, ff.p_char) == other.polynomial));

    ff_element inverse = other.make_inverse();

    return (*this) * inverse;
};

ff_element ff_element::make_inverse() const {
    if (this->polynomial == PolynomialMod({1}, ff.p_char)){
        return *this;
    }
    return ff_element(euclidean_algorithm(ff.min_polynomial, polynomial), ff);
}

PolynomialMod ff_element::euclidean_algorithm(PolynomialMod elemA, PolynomialMod elemB) const {
    int p = ff.p_char;

    PolynomialMod A = elemA;
    PolynomialMod B = elemB;

    if(B.degree > A.degree) {
        return euclidean_algorithm(elemB, elemA);
    }

    std::tuple<PolynomialMod, PolynomialMod> compute_A = std::make_pair(PolynomialMod({1}, p), PolynomialMod(p));
    std::tuple<PolynomialMod, PolynomialMod> compute_B = std::make_pair(PolynomialMod(p), PolynomialMod({1}, p));

    assert(!(B == PolynomialMod(p)));

    while(true) {
        auto [Q, R] = PolynomialMod::div(A, B);

        PolynomialMod compute_r_0 = std::get<0>(compute_A) - std::get<0>(compute_B)*Q;
        PolynomialMod compute_r_1 = std::get<1>(compute_A) - std::get<1>(compute_B)*Q;

        if(R.degree == 0) {
            return compute_r_1/R.lead_coeff;
        }

        compute_B = std::make_pair(compute_r_0, compute_r_1);

        A = B;
        B = R;
    }
}

TEST_CASE("ff_element") {
    int p = 5;
    finite_field ff = finite_field(p, 3);
    ff.min_polynomial = PolynomialMod({1,2,0,1}, p);

    PolynomialMod p1 = PolynomialMod({1,2,3}, p);
    ff_element elem1 = ff_element(p1, ff);

    PolynomialMod p2 = PolynomialMod({2,4,6}, p);
    ff_element elem2 = ff_element(p2, ff);

    PolynomialMod p3 = PolynomialMod({2,1,2}, p);
    ff_element elem3 = ff_element(p3, ff);

    ff_element zero_e = ff_element(PolynomialMod({0}, p), ff);

    SUBCASE("adding ff_element") {
        CHECK((elem1 + elem1) == elem2);
    }

    SUBCASE("subtracting ff_element") {
        CHECK((elem1 - elem1) == zero_e);
    }

    SUBCASE("multiplication") {
        PolynomialMod p4 = PolynomialMod({0,0,3}, p);
        ff_element elem4 = ff_element(p4, ff);

        CHECK(elem1*elem3 == elem4);
    }

    SUBCASE("division") {
        PolynomialMod p4 = PolynomialMod({2,2,4}, p);
        ff_element elem4 = ff_element(p4, ff);

        CHECK(elem1/elem3 == elem4);
    }
}
