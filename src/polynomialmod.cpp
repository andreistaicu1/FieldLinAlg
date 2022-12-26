#include "polynomialmod.hpp"
#include "integers.hpp"
#include "doctest.h"

PolynomialMod::PolynomialMod(Polynomial poly, int p) : p(p) {
    this->poly = poly.pure_mod(p);
    lead_coeff = this->poly.lead_coeff;
    degree = this->poly.degree;
}

PolynomialMod::PolynomialMod(std::vector<int> vec, int p) : p(p) {
    this->poly = Polynomial(vec).pure_mod(p);
    lead_coeff = this->poly.lead_coeff;
    degree = this->poly.degree;
}

PolynomialMod::PolynomialMod(int p) : p(p) {
    this->poly = Polynomial().pure_mod(p);
    lead_coeff = this->poly.lead_coeff;
    degree = this->poly.degree;
}

PolynomialMod PolynomialMod::operator=(const PolynomialMod& other) {
    poly = other.poly;
    p = other.p;
    lead_coeff = other.lead_coeff;
    degree = other.degree;

    return (*this);
}

bool PolynomialMod::operator==(const PolynomialMod& other) const {
    return (poly == other.poly && p == other.p);
}

PolynomialMod PolynomialMod::operator+(const PolynomialMod& other) const {
    return PolynomialMod((poly+other.poly).pure_mod(p), p);
}
PolynomialMod PolynomialMod::operator-(const PolynomialMod& other) const {
    return PolynomialMod((poly-other.poly).pure_mod(p), p);
}

PolynomialMod PolynomialMod::operator*(const PolynomialMod& other) const {
    return PolynomialMod(poly*other.poly, p);
}

PolynomialMod PolynomialMod::operator/(const PolynomialMod& B) const {
    return std::get<0>(div(*this, B));
}

PolynomialMod PolynomialMod::derivative() const {
    return PolynomialMod(poly.derivative().pure_mod(p), p);
}

std::tuple<PolynomialMod, PolynomialMod> PolynomialMod::div(const PolynomialMod& A, const PolynomialMod& B) {
    int p = A.p;

    PolynomialMod R = A;
    PolynomialMod Q = PolynomialMod(p);

    int lead_B = B.lead_coeff;
    int lead_B_inv = inverse_mod_p(lead_B, p);

    while(true) {
        if(R.degree < B.degree || R == PolynomialMod(p)) {
            return std::make_pair(Q,R);
        }

        int lead_R = R.lead_coeff;
        int lead_R_div_lead_B = mod_p(lead_R*lead_B_inv, p);

        std::vector<int> S_vec(R.degree-B.degree+1,0);
        S_vec[R.degree-B.degree] = lead_R_div_lead_B;
        PolynomialMod S = PolynomialMod(S_vec, p);

        Q = Q + S;
        R = R - (S*B);
    }
}

PolynomialMod PolynomialMod::gcd(const PolynomialMod& Aorig, const PolynomialMod& Borig) {
    int p = Aorig.p;

    PolynomialMod A = Aorig;
    PolynomialMod B = Borig;

    while(true) {
        if(B == PolynomialMod(p)) {
            return A;
        }

        auto [Q, R] = div(A, B);
        A = B;
        B = R;
    }
}

PolynomialMod PolynomialMod::mod(const PolynomialMod& modulus) {
    auto [Q, R] = div((*this), modulus);
    return R;
}

void PolynomialMod::print() const {
    poly.print();
}

PolynomialMod PolynomialMod::operator*(int scalar) const {
    return PolynomialMod(poly*scalar, p);
}

PolynomialMod PolynomialMod::operator/(int scalar) const {
    int inverse = inverse_mod_p(scalar, p);
    return PolynomialMod(poly*inverse, p);
}

PolynomialMod operator*(int scalar, const PolynomialMod& poly) {
    return poly*scalar;
}

PolynomialMod PolynomialMod::exp(int n) const {
    assert(n > 0);

    PolynomialMod y = PolynomialMod({1}, p);
    if(n == 0) {
        return y;
    }

    int N = n;
    PolynomialMod z = (*this);

    while(true) {
        if(N % 2 == 1) {
            y = z * y;
        }

        N = N/2;

        if(N == 0) {
            return y;
        }
        z = z*z;
    }
}

TEST_CASE("exp") {
    int p = 7;
    PolynomialMod p1 = PolynomialMod({1,2,3}, 7);
    PolynomialMod output = PolynomialMod({1,6,0,2,0,5,6}, 7);

    CHECK(p1.exp(3) == output);
}