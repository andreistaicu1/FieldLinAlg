#include "integers.hpp"
#include "doctest.h"

int mod_p(int val, int p) {
    return ((val % p) + p) % p;
}

int euclidean_algorithm(int a, int b, std::tuple<int, int> compute_a, std::tuple<int, int> compute_b, int p) {
    std::div_t div_n = std::div(a, b);
    int q = div_n.quot;
    int r = div_n.rem;

    int compute_r_0 = std::get<0>(compute_a) - (std::get<0>(compute_b) * q);
    int compute_r_1 = std::get<1>(compute_a) - (std::get<1>(compute_b) * q);

    if (r == 1) {
        return compute_r_1;
    }

    auto compute_r = std::make_pair(compute_r_0, compute_r_1);
    return euclidean_algorithm(b, r, compute_b, compute_r, p);
}

int inverse_mod_p(int val, int p) {
    assert(mod_p(val, p) != 0);

    if (mod_p(val, p) == 1){
        return val;
    }
    int inverse = euclidean_algorithm(p, mod_p(val, p), std::make_pair(1, 0), std::make_pair(0, 1), p);
    return mod_p(inverse, p);
}

TEST_CASE("integers") {
    SUBCASE("inverse_mod_p") {
        CHECK(inverse_mod_p(1, 5) == 1);
        CHECK(inverse_mod_p(2, 5) == 3);
        CHECK(inverse_mod_p(4, 5) == 4);
        CHECK(inverse_mod_p(39, 5) == 4);
        CHECK(inverse_mod_p(113, 5) == 2);
    }

    SUBCASE("mod_p") {
        CHECK(mod_p(-3, 5) == 2);
    }

    SUBCASE("euclidean_algorithm") {
    }
}