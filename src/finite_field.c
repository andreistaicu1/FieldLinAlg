//
// Created by Andrei Staicu on 12/15/22.
//

#include "../headers/finite_field.h"
#include "../headers/min_polynomial.h"


typedef struct finite_field {
    size_t p_char;
    size_t degree;
    size_t *min_polynomial;
}finite_field_t;

typedef struct ff_element {
    finite_field_t *ff;
    size_t *polynomial;
}ff_element_t;

/*
* Creates finite field
*/
finite_field_t *init_finite_field(size_t degree, size_t p_char){
    finite_field_t *ff = malloc(sizeof(finite_field_t));
    ff->degree = degree;
    ff->p_char = p_char;
    ff->min_polynomial = find_min_polynomial(p_char, degree);
    return ff;
}

/*
* Frees finite field
*/
void free_finite_field(finite_field_t * finite_field){
    free(finite_field->min_polynomial);
    free(finite_field);
}

/*
* Creates finite field element
*/
void init_ff_element(finite_field_t *finite_field){
    return;
}


bool ff_equal(finite_field_t *ff_1, finite_field_t *ff_2){
    return ff_1 == ff_2;
}


ff_element_t add(ff_element_t a, ff_element_t b){

    assert(ff_equal(a.ff, b.ff));

    finite_field_t *ff = a.ff;
    size_t new_poly[a.ff->degree];
    size_t carry;

    carry = 0;
    for (size_t i = 0; i < ff->degree; i++){
        size_t val = a.polynomial[i] + b.polynomial[i] + carry;
        carry = floor(val / ff->p_char);
        new_poly[i] = val;
    }

    return (ff_element_t) {.ff = ff, .polynomial = new_poly};
}

// Subtracts two elements of the finite field
//      a - b
void subtract(ff_element_t a, ff_element_t b){
    return;
}

// Multiplies two elements of the finite field
//      a * b
void multiply(ff_element_t a, ff_element_t b){
    return;
}

// Divides two elements of the finite field
//      a / b
void divide(ff_element_t a, ff_element_t b){
    return;
}

// Finds the multiplicative inverse of an element
//      Finds b s.t. a * b = 1
void inverse(ff_element_t a){
    return;
}