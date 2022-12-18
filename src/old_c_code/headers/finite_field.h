//
// Created by Andrei Staicu on 12/15/22.
//

#ifndef FINITE_FIELD_H
#define FINITE_FIELD_H

#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include "assert.h"
#include <math.h>

typedef struct finite_field finite_field_t;

typedef struct ff_element ff_element_t;

/*
* Creates finite field
*/
finite_field_t *init_finite_field();

/*
* Frees finite field
*/
void free_finite_field(finite_field_t *finite_field);

/*
* Creates finite field element
*/
void init_ff_element(finite_field_t *finite_field);

/*
* Checks if two given finite fields are equal
*/ 
bool ff_equal(finite_field_t *ff_1, finite_field_t *ff_2);

/*
* Adds two elements of the finite field
*/
ff_element_t add(ff_element_t a, ff_element_t b);

/*
* Subtracts two elements of the finite field
*/
void subtract(ff_element_t a, ff_element_t b);

/*
* Multiplies two elements of the finite field
*/
void multiply(ff_element_t a, ff_element_t b);

/*
* Divides two elements of the finite field
*/
void divide(ff_element_t a, ff_element_t b);

/*
* Finds the multiplicative inverse of an element
*/
void inverse(ff_element_t a);



#endif //FINITE_FIELD_H
