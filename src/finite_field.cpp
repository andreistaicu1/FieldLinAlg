#include "finite_field.hpp"
#include "doctest.h"

finite_field::finite_field() {
    this->p_char = 0;
    this->degree = 0;
    this->min_polynomial = Polynomial();
}

finite_field::finite_field(int p_char, int degree) : p_char(p_char), degree(degree) {
    this->min_polynomial = Polynomial::zero_polynomial(1);
};

bool finite_field::operator== (finite_field const &other) const {
    return this->p_char == other.p_char && this->degree == other.degree;
}
    
ff_element::ff_element(Polynomial polynomial, finite_field ff) {
    assert(polynomial.degree <= ff.degree);
    polynomial.mod(ff.p_char);
    this->polynomial = polynomial.trim_end();
    this->ff = ff;
};

ff_element::ff_element(finite_field ff) : polynomial(Polynomial::zero_polynomial(ff.degree)), ff(ff) {};

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

    Polynomial output = this->polynomial + other.polynomial;
    output.mod(this->ff.p_char);

    return ff_element(output, ff);
};

ff_element ff_element::operator- (ff_element const &other) const {
    assert(this->ff == other.ff);

    Polynomial output = this->polynomial - other.polynomial;
    output.mod(this->ff.p_char);

    return ff_element(output, this->ff);
};

ff_element ff_element::operator* (ff_element const &other) const {
    assert(this->ff == other.ff);

    std::tuple<Polynomial, Polynomial> output = Polynomial::division_mod_p(polynomial*other.polynomial, ff.min_polynomial, ff.p_char);
    Polynomial remainder = std::get<1>(output);

    return ff_element(remainder, ff);
};

ff_element ff_element::operator/ (ff_element const &other) const {
    assert(this->ff == other.ff);
    assert(!(Polynomial({0}) == other.polynomial));

    ff_element inverse = other.make_inverse();

    return (*this) * inverse;
};

ff_element ff_element::make_inverse() const {
    if (this->polynomial == Polynomial({1})){
        return *this;
    }
    return ff_element(euclidean_algorithm(ff.min_polynomial, polynomial), ff);
}

Polynomial ff_element::euclidean_algorithm(Polynomial elemA, Polynomial elemB) const {
    Polynomial A = elemA;
    Polynomial B = elemB;

    if(B.degree > A.degree) {
        return euclidean_algorithm(elemB, elemA);
    }

    std::tuple<Polynomial, Polynomial> compute_A = std::make_pair(Polynomial({1}), Polynomial());
    std::tuple<Polynomial, Polynomial> compute_B = std::make_pair(Polynomial(), Polynomial({1}));

    assert(!(B == Polynomial()));

    while(true) {
        auto division = Polynomial::pseudo_div(A, B);

        int delta = std::get<2>(division);
        int delta_inverse = inverse_mod_p(((delta % ff.p_char)+ ff.p_char) % ff.p_char, ff.p_char);

        Polynomial Q = delta_inverse*std::get<0>(division);
        Polynomial R = delta_inverse*std::get<1>(division);
        Q.mod(ff.p_char);
        R.mod(ff.p_char);

        Polynomial compute_r_0 = std::get<0>(compute_A) - std::get<0>(compute_B)*Q;
        Polynomial compute_r_1 = std::get<1>(compute_A) - std::get<1>(compute_B)*Q;
        compute_r_0.mod(ff.p_char);
        compute_r_1.mod(ff.p_char);

        if(R.degree == 0) {
            return inverse_mod_p(R.lead_coeff, ff.p_char)*compute_r_1;
        }

        compute_B = std::make_pair(compute_r_0, compute_r_1);

        A = B;
        B = R;
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

    ff_element zero_e = ff_element(Polynomial({0}), ff);

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
