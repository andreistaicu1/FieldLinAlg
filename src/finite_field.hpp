#ifndef FINITEFIELD_H
#define FINITEFIELD_H

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

#include "polynomial.hpp"
#include "min_polynomial.hpp"
#include "integers.hpp"

class finite_field {
public:
    int p_char, degree;
    Polynomial min_polynomial;

    finite_field(int p_char, int degree);
    finite_field();

    bool operator==(finite_field const &other) const;
};

class ff_element {
public:
    Polynomial polynomial;
    finite_field ff;

    ff_element(Polynomial polynomial, finite_field ff);
    ff_element(finite_field ff);

    bool operator== (ff_element const &other) const;
    ff_element operator= (ff_element const &other);
    ff_element operator+ (ff_element const &other) const;
    ff_element operator- (ff_element const &other) const;
    ff_element operator* (ff_element const &other) const;
    ff_element operator/ (ff_element const &other) const;
    ff_element make_inverse() const;
    Polynomial euclidean_algorithm(Polynomial elemA, Polynomial elemB) const;
};

#endif