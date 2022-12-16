#ifndef MIN_POLYNOMIAL_H
#define MIN_POLYNOMIAL_H

#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include "assert.h"

/*
* Given a degree, it finds the minimum polynomial
*/

size_t *find_min_polynomial(size_t degree, size_t p_char);

#endif //MIN_POLYNOMIAL_H