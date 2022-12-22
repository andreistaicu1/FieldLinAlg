#ifndef MINPOLYNOMIAL_H
#define MINPOLYNOMIAL_H

#include <vector>
#include <tuple>
#include "polynomial.hpp"
#include "integers.hpp"

extern std::vector<std::tuple<int, Polynomial>> squarefree_factorization(Polynomial A, int p);

#endif