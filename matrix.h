#ifndef MATRIX_H
#define MATRIX_H

#include <stdint.h>

/* utility functions */

uint32_t fast_rand(void);

void set_seed(uint32_t value);
void set_nthreads(ssize_t count);
void set_dimensions(ssize_t width);
void *thr_worker(void *arg); 
void wait_flip(void);

void display(const uint32_t* matrix);
void display_row(const uint32_t* matrix, ssize_t row);
void display_column(const uint32_t* matrix, ssize_t column);
void display_element(const uint32_t* matrix, ssize_t row, ssize_t column);

/* matrix operations */

uint32_t* new_matrix(void);
uint32_t* identity_matrix(void);
uint32_t* random_matrix(uint32_t seed);
uint32_t* uniform_matrix(uint32_t value);
uint32_t* sequence_matrix(uint32_t start, uint32_t step);

uint32_t* cloned(const uint32_t* matrix);
uint32_t* reversed(const uint32_t* matrix);
uint32_t* transposed(const uint32_t* matrix);

uint32_t* scalar_add(const uint32_t* matrix, uint32_t scalar);
uint32_t* scalar_mul(const uint32_t* matrix, uint32_t scalar);
uint32_t* matrix_pow(const uint32_t* matrix, uint32_t exponent);
uint32_t* matrix_add(const uint32_t* matrix_a, const uint32_t* matrix_b);
uint32_t* matrix_mul(const uint32_t* matrix_a, const uint32_t* matrix_b);

void scalar_add_worker(const size_t id);
void scalar_mul_worker(const size_t id);
void matrix_add_worker(const size_t id);
void strass_add_worker(const size_t id);
void strass_sub_worker(const size_t id);
void strass_mul_worker(const size_t id);

void matrix_mul_nomem(const uint32_t* matrix_a, const uint32_t* matrix_b, uint32_t* result);

void strass_add(const uint32_t* a, ssize_t aw, ssize_t ah, ssize_t as,
                     const uint32_t* b, ssize_t bw, ssize_t bh, ssize_t bs,
                     uint32_t* c, ssize_t cs);
void strass_sub(const uint32_t* a, ssize_t aw, ssize_t ah, ssize_t as,
                     const uint32_t* b, ssize_t bw, ssize_t bh, ssize_t bs,
                     uint32_t* c, ssize_t cs);
void strass_mul(const uint32_t* a, ssize_t aw, ssize_t ah, ssize_t as,
                     const uint32_t* b, ssize_t bw, ssize_t bh, ssize_t bs,
                     uint32_t* c, ssize_t cs);
void strassen(const uint32_t* a, ssize_t aw, ssize_t ah, ssize_t as,
                   const uint32_t* b, ssize_t bw, ssize_t bh, ssize_t bs,
                   uint32_t* c, ssize_t cs);

/* compute operations */

uint32_t get_sum(const uint32_t* matrix);
uint32_t get_trace(const uint32_t* matrix);
uint32_t get_minimum(const uint32_t* matrix);
uint32_t get_maximum(const uint32_t* matrix);
uint32_t get_frequency(const uint32_t* matrix, uint32_t value);

#endif
