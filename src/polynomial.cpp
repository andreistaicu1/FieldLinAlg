#include <vector>
#include <memory>
#include <cmath>
#include <numeric>

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <iostream>

#include "min_polynomial.cpp"

#include <tuple>

#include "polynomial.hpp"

#include "doctest.h"

Polynomial::Polynomial(std::vector<int> p) : p(p) {
    int last = 0;
    for(auto i = 0; i < this->p.size(); ++i) {
        if(this->p[i] != 0) {
            last = i;
        }
    }
    this->degree = last;
    this->lead_coeff = this->p[this->degree];
};

Polynomial::Polynomial(){
    std::vector<int> output(1, 0);
    this->p = output;
    this->degree = 0;
    this->lead_coeff = 0;
}

bool Polynomial::operator==(const Polynomial& other) const {

    if (this->degree != other.degree){
        return false;
    }

    for (auto i = 0; i < degree+1; i++) {
        if (this->p[i] != other.p[i]) return false;
    }

    return true;
}

void Polynomial::operator=(const Polynomial& other) {
    p = other.p;
    degree = other.degree;
    lead_coeff = other.lead_coeff;
}

Polynomial Polynomial::operator+(const Polynomial& other) const {
    auto max = std::max(this->p.size(), other.p.size());

    std::vector<int> output(max,0);
    for(auto i = 0; i < max; ++i) {
        int first = 0;
        if(this->p.size() > i) {
            first = this->p[i];
        }
        
        int second = 0;
        if(other.p.size() > i) {
            second = other.p[i];
        }

        output[i] = first+second;
    }

    return Polynomial(output);
};

Polynomial Polynomial::operator-(const Polynomial& other) const {
    auto max = std::max(this->p.size(), other.p.size());

    std::vector<int> output(max,0);
    for(auto i = 0; i < max; ++i) {
        int first = 0;
        if(this->p.size() > i) {
            first = this->p[i];
        }
        
        int second = 0;
        if(other.p.size() > i) {
            second = other.p[i];
        }

        output[i] = first-second;
    }

    return Polynomial(output);
};

Polynomial Polynomial::operator*(const Polynomial& other) const {
    
    int new_degree = this->degree + other.degree + 1;
    std::vector<int> output(new_degree, 0);

    for(auto n = 0; n < new_degree; n++){
        for(auto i = 0; i <= n; i++){
            if (i <= this->degree && n - i <= other.degree){
                output[n] += this->p[i] * other.p[n-i];
            }
        }
    }
    
    return Polynomial(output);
}

Polynomial Polynomial::operator*(int scalar) const {
    std::vector<int> output(this->p.size(), 0);

    for(auto i = 0; i < this->p.size(); ++i) {
        output[i] = this->p[i]*scalar;
    }

    return Polynomial(output);
}

Polynomial Polynomial::operator/(int scalar) const {
    std::vector<int> output(this->p.size(), 0);

    for(auto i = 0; i < p.size(); ++i) {
        output[i] = p[i]/scalar;
    }

    return output;
}

Polynomial Polynomial::trim_end() const {
    std::vector<int> output;
    for(auto i = 0; i < degree+1; ++i) {
        output.push_back(this->p[i]);
    }

    return Polynomial(output);
}

Polynomial Polynomial::pad_to_length(int length) {
    std::vector<int> output;
    for(auto i = 0; i < length; ++i) {
        if(i < this->p.size()) {
            output.push_back(this->p[i]);
        }
        else {
            output.push_back(0);
        }
    }
    return Polynomial(output);
}

void Polynomial::mod(int p) {
    for (auto i = 0; i < this->p.size(); i++){
        this->p[i] = ((this->p[i] % p) + p) % p;
    }
}

int Polynomial::size() {
    return this->p.size();
}

int Polynomial::cont() {
    int result = p[0];
    for(auto i = 1; i < degree; ++i) {
        result = std::gcd(result, p[i]);
    }

    return result;
}
void Polynomial::print() const {
    for(auto elem : this->p) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;
}

Polynomial Polynomial::zero_polynomial(int length) {
    std::vector<int> output(length, 0);
    return Polynomial(output);
}

/*
Polynomial Polynomial::euclidean_algorithm(Polynomial p1, Polynomial p2, std::tuple<Polynomial, Polynomial> compute_p1, 
    std::tuple<Polynomial, Polynomial> compute_p2){

        auto output = polynomial_div(p1, p2);
        Polynomial quotient = std::get<0>(output);
        Polynomial remainder = std::get<1>(output);

        // This function should only be input relatively prime things
        assert(!(remainder == zero_polynomial(1)));

        std::vector<int> one_vec{1};
        Polynomial one = Polynomial(one_vec);

        // Checks if the remainder is one
        if (remainder == one){
            return std::get<1>(compute_p1) - (std::get<1>(compute_p2) * quotient);
        }

        Polynomial compute_rem_a = std::get<0>(compute_p1) - (std::get<0>(compute_p2) * quotient);
        Polynomial compute_rem_b = std::get<1>(compute_p1) - (std::get<1>(compute_p2) * quotient);
        return euclidean_algorithm(p2, remainder, compute_p2, std::make_pair(compute_rem_a, compute_rem_b));
}
*/

Polynomial operator*(int scalar, const Polynomial& poly) {
    return poly*scalar;
}

// Returns a tuple (Q, R) such that
// d^(m-n+1)*A = B*Q + R
// where m is the degree of A, n is the degree of B, d is the leading coefficient of B,
// and the degree of R < degree of B
std::tuple<Polynomial, Polynomial, int> Polynomial::pseudo_div(Polynomial A, Polynomial B) {
    int m = A.degree;
    int n = B.degree;
    int d = B.lead_coeff;

    if(m < n) {
        return std::make_tuple(zero_polynomial(1), A, 1);
    }

    Polynomial R = A;
    Polynomial Q = Polynomial({0});
    int e = m - n + 1;

    while(true) {
        std::cout << "loop" << std::endl;
        Q.print();
        R.print();
        B.print();
        if(R.degree < B.degree) {
            int q = std::pow(d, e);
            Q = q*Q;
            R = q*R;

            return std::make_tuple(Q, R, std::pow(d, m-n+1));
        }

        std::vector<int> S_vec(R.degree-B.degree+1, 0);
        S_vec[R.degree-B.degree] = R.lead_coeff;
        Polynomial S = Polynomial(S_vec);

        Q = d*Q + S;
        R = d*R - S*B;
        e = e-1;
    }
}

TEST_CASE("pseudo_div") {
    Polynomial p1 = Polynomial({1,3});
    Polynomial p2 = Polynomial({4,0,0,0});
    auto triple = Polynomial::pseudo_div(p1, p2);
    std::get<0>(triple).print();
    std::get<1>(triple).print();
    std::cout << std::get<2>(triple) << std::endl;
}

// Returns the GCD of A and B
Polynomial Polynomial::gcd(Polynomial A, Polynomial B) {
    if(B.degree > A.degree) {
        return gcd(B, A);
    }

    if(B == zero_polynomial(1)) {
        return A;
    }

    int a = A.cont(); // Do we even need GCD since we are mod p?
    int b = B.cont();
    int d = std::gcd(a, b);

    A = A/a;
    B = B/b;

    while(true) {
        Polynomial R = std::get<1>(pseudo_div(A, B));

        if(R == zero_polynomial(1)) {
            return d*B;
        }
        if(R.degree == 0) {
            return d*Polynomial({1});
        }

        A = B;
        B = R/R.cont();
    }
}
