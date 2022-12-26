#ifndef POLYNOMIALMOD_H
#define POLYNOMIALMOD_H

#include "polynomial.hpp"

#include <tuple>

class PolynomialMod {
public:
    Polynomial poly;
    int p;
    int lead_coeff;
    int degree;

    PolynomialMod(Polynomial poly, int p);
    PolynomialMod(std::vector<int> vec, int p);
    PolynomialMod(int p);

    PolynomialMod operator=(const PolynomialMod&);
    bool operator==(const PolynomialMod&) const;
    PolynomialMod operator+(const PolynomialMod&) const;
    PolynomialMod operator-(const PolynomialMod&) const;
    PolynomialMod operator*(const PolynomialMod&) const;
    PolynomialMod operator/(const PolynomialMod&) const;
    PolynomialMod derivative() const;

    static std::tuple<PolynomialMod, PolynomialMod> div(const PolynomialMod&, const PolynomialMod&);
    static PolynomialMod gcd(const PolynomialMod&, const PolynomialMod&);

    PolynomialMod mod(const PolynomialMod& modulus);

    PolynomialMod operator*(int scalar) const;
    PolynomialMod operator/(int scalar) const;

    PolynomialMod exp(int n) const;

    void print() const;
};

extern PolynomialMod operator*(int scalar, const PolynomialMod& poly);
#endif