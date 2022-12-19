#ifndef polynomialH
#define polynomialH
#include <tuple>
#include <vector>
#include "doctest.h"

const std::vector<int> ZERO = {0};
const std::vector<int> ONE = {1};

class Polynomial {
public:
    std::vector<int> p;
    int degree;
    int lead_coeff;

    Polynomial(std::vector<int> p);
    Polynomial();

    bool operator==(const Polynomial& other) const;
    void operator=(const Polynomial& other);
    Polynomial operator+(const Polynomial& other);
    Polynomial operator-(const Polynomial& other);
    Polynomial operator*(const Polynomial& other);
    Polynomial trim_end() const;
    Polynomial pad_to_length(int length);
    Polynomial operator*(int scalar) const;
    Polynomial operator/(int scalar) const;
    void mod(int p);
    int cont();
    int size();
    void print();
    static Polynomial zero_polynomial(int length);
    static std::tuple<Polynomial, Polynomial, int> pseudo_div(Polynomial A, Polynomial B);
    static Polynomial gcd(Polynomial A, Polynomial B);
};

Polynomial operator*(int scalar, const Polynomial& poly);
#endif