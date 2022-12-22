#ifndef INTEGERS_H
#define INTEGERS_H
#include <tuple>
#include "doctest.h"

extern int euclidean_algorithm(int a, int b, std::tuple<int, int> compute_a, std::tuple<int, int> compute_b, int p);
extern int inverse_mod_p(int val, int p);
extern int mod_p(int val, int p);

#endif