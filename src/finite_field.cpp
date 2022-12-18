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


struct PolynomialHasher {
    int operator()(const Polynomial &polynomial) const {
        const std::vector<int> trimmed = polynomial.trim_end().p;
        int hash = trimmed.size();
        for(auto &i : trimmed) {
            hash ^= i + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

class finite_field {
public:
    int p_char, degree;
    Polynomial min_polynomial;
    std::unordered_map<Polynomial, void *, PolynomialHasher> elements;
    std::unordered_map<void *, void *> inverses;

    ~finite_field() {
        for(auto it : elements) {
            free(it.second);
        }
    }

    finite_field() {
        this->p_char = 0;
        this->degree = 0;
        this->min_polynomial = Polynomial();
    }
    
    finite_field(int p_char, int degree) : p_char(p_char), degree(degree) {
        this->min_polynomial = Polynomial::zero_polynomial(1);
    };

    bool operator== (finite_field const &other) {
        return this->p_char == other.p_char && this->degree == other.degree;
    }

    void set_element(Polynomial poly, void *element){
        this->elements[poly] = element;
    }

    void * get_element(Polynomial poly){
        if (elements.find(poly) != elements.end()){
            return elements[poly];
        }
        return nullptr;
    }

    void set_inverse(void *element, void *inverse){
        this->inverses[element] = inverse;
    }

    void * get_inverse(void *element){
        if (inverses.find(element) != inverses.end()){
            return inverses[element];
        }
        return nullptr;
    }

    int inverse_mod_p(int val) {
        if (val == 1){
            return val;
        }
        int inverse = euclidean_algorithm(p_char, val, std::make_pair(1, 0), std::make_pair(0, 1));
        return (((inverse % p_char) + p_char) % p_char);
    }

    int euclidean_algorithm(int a, int b, std::tuple<int, int> compute_a, std::tuple<int, int> compute_b){
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
private:
    ff_element(Polynomial polynomial, finite_field ff) {
        assert(polynomial.p.size() <= ff.degree);
        this->polynomial = polynomial;
        this->ff = ff;
    };

    ff_element(finite_field ff) : polynomial(Polynomial::zero_polynomial(ff.degree)), ff(ff) {};

public:
    Polynomial polynomial;
    finite_field ff;

    static ff_element get_ff_element(Polynomial polynomial, finite_field ff) {
        ff_element * element = (ff_element *) ff.get_element(polynomial);
        if(element != nullptr) {
            return *element;
        }

        ff_element* new_elem = new ff_element(polynomial, ff);
        ff.set_element(polynomial, new_elem);
        return *new_elem;
    }

    ff_element operator+ (ff_element const &other) {
        assert(this->ff == other.ff);

        Polynomial output = this->polynomial + other.polynomial;
        output.mod(this->ff.p_char);

        return get_ff_element(output, ff);
    };

    ff_element operator- (ff_element const &other) {
        assert(this->ff == other.ff);

        Polynomial output = this->polynomial - other.polynomial;
        output.mod(this->ff.p_char);

        return get_ff_element(output, this->ff);
    };

    ff_element operator* (ff_element const &other) {
        assert(this->ff == other.ff);

        Polynomial p = Polynomial::zero_polynomial(1);

        std::tuple<Polynomial, Polynomial, int> triple = Polynomial::pseudo_div(polynomial*other.polynomial, ff.min_polynomial);

        Polynomial remainder = std::get<1>(triple);
        int c = std::get<2>(triple);
        remainder.mod(this->ff.p_char);

        Polynomial output = remainder * (this->ff.inverse_mod_p(c % this->ff.p_char));
        output.mod(this->ff.p_char);

        return get_ff_element(output, this->ff);
    };

    ff_element operator/ (ff_element const &other) {
        assert(this->ff == other.ff);
        
        ff_element output = ff_element(this->ff);
        
    };
};