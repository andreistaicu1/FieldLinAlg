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

class finite_field {
public:
    int p_char, degree;
    std::vector<int> min_polynomial;
    
    finite_field(int p_char, int degree) : p_char(p_char), degree(degree) {
        std::vector<int> temp(10, 0);
        this->min_polynomial = temp;
    };

    bool operator== (finite_field const &other) {
        return this->p_char == other.p_char && this->degree == other.degree;
    }

    int inverse_mod_p(int val) {
        assert(val % this->p_char == val);
        assert(val > 0);

        int inverse = euclidean_algorithm(p_char, val, std::make_pair(1, 0), std::make_pair(0, 1));
        return (((inverse % this->p_char) + this->p_char) % this->p_char);
    }

    /*
    * Given p and any x in {0, ..., p-1} it computes "a" such that ax = 1 mod p
    * 
    * It first finds q and r such that a = qb + r, if r = 1 it returns (-q) if not it recurses
    * on the input (b, r, (0, 1), (0, 1), (1, -q)), each tuple shows how to make the first or second
    * argument as a linear combination of a and b at each step.
    * 
    * Clearly at the final step we have a linear combination of a b to get 1. We only care about the
    * coefficient of the "b", so we return that.
    */
    int euclidean_algorithm(int a, int b, std::tuple<int, int> compute_a, std::tuple<int, int> compute_b){
        int q = std::floor(a / b);
        int r = a - (b * q);

        int compute_r_0 = std::get<0>(compute_a) - (std::get<0>(compute_b) * q);
        int compute_r_1 = std::get<1>(compute_b) - (std::get<1>(compute_b) * q);

        if (r == 1) {
            return compute_r_1;
        }

        auto compute_r = std::make_pair(compute_r_0, compute_r_1);
        return euclidean_algorithm(b, r, compute_b, compute_r);
    }
    
};

class ff_element {
public:
    std::vector<int> polynomial;
    finite_field ff;

    ff_element(std::vector<int> polynomial, finite_field ff) : polynomial(polynomial), ff(ff) {};
    ff_element(finite_field ff) : ff(ff) {
        std::vector<int> zeroes(ff.degree, 0);
        this->polynomial = zeroes;
    };

    ff_element operator+ (ff_element const &other) {
        assert(this->ff == other.ff);

        ff_element output = ff_element(this->ff);
        for (auto i = 0; i < this->ff.degree; i++){
            output.polynomial[i] = (this->polynomial[i]+other.polynomial[i]) % this->ff.p_char;
        }
        
        return output;
    };

    ff_element operator- (ff_element const &other) {
        assert(this->ff == other.ff);
        
        ff_element output = ff_element(this->ff);
        for (auto i = 0; i < this->ff.degree; i++){
            output.polynomial[i] = ((this->polynomial[i] - other.polynomial[i]) % this->ff.p_char + this->ff.p_char) % this->ff.p_char;
        }

        return output;
    };

    ff_element operator* (ff_element const &other) {
        assert(this->ff == other.ff);
        
        ff_element output = ff_element(this->ff);
    };

    ff_element operator/ (ff_element const &other) {
        assert(this->ff == other.ff);
        
        ff_element output = ff_element(this->ff);
        
    };
};