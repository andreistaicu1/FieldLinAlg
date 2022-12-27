#ifndef MINPOLYNOMIAL_H
#define MINPOLYNOMIAL_H

#include <vector>
#include <tuple>
#include <iostream>
#include "polynomialmod.hpp"
#include "integers.hpp"

extern std::vector<std::tuple<int, PolynomialMod>> squarefree_factorization(PolynomialMod A);
extern std::vector<PolynomialMod> cantor_zassenhaus_split(PolynomialMod input_poly, int d);

#endif