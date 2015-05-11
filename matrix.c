#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <stdbool.h>
#include <pthread.h>
#include <inttypes.h>

#include <string.h>

#include "matrix.h"

#define STRASS_THRESH 200

static uint32_t g_seed = 0;

static ssize_t g_width = 0;
static ssize_t g_height = 0;
static ssize_t g_elements = 0;

static ssize_t g_nthreads = 1;

////////////////////////////////
///     UTILITY FUNCTIONS    ///
////////////////////////////////

/**
 * Returns pseudorandom number determined by the seed
 */
uint32_t fast_rand(void) {

    g_seed = (214013 * g_seed + 2531011);
    return (g_seed >> 16) & 0x7FFF;
}

/**
 * Sets the seed used when generating pseudorandom numbers
 */
void set_seed(uint32_t seed) {

    g_seed = seed;
}

/**
 * Sets the number of threads available
 */
void set_nthreads(ssize_t count) {

    g_nthreads = count;
}

/**
 * Sets the dimensions of the matrix
 */
void set_dimensions(ssize_t order) {

    g_width = order;
    g_height = order;

    g_elements = g_width * g_height;
}

/**
 * Displays given matrix
 */
void display(const uint32_t* matrix) {

    for (ssize_t y = 0; y < g_height; y++) {
        for (ssize_t x = 0; x < g_width; x++) {
            if (x > 0) printf(" ");
            printf("%" PRIu32 " ", matrix[y * g_width + x]);
        }

        printf("\n");
    }
}

/**
 * Displays given matrix row
 */
void display_row(const uint32_t* matrix, ssize_t row) {

    for (ssize_t x = 0; x < g_width; x++) {
        if (x > 0) printf(" ");
        printf("%" PRIu32, matrix[row * g_width + x]);
    }

    printf("\n");
}

/**
 * Displays given matrix column
 */
void display_column(const uint32_t* matrix, ssize_t column) {

    for (ssize_t i = 0; i < g_height; i++) {
        printf("%" PRIu32 "\n", matrix[i * g_width + column]);
    }
}

/**
 * Displays the value stored at the given element index
 */
void display_element(const uint32_t* matrix, ssize_t row, ssize_t column) {

    printf("%" PRIu32 "\n", matrix[row * g_width + column]);
}

////////////////////////////////
///   MATRIX INITALISATIONS  ///
////////////////////////////////

/**
 * Returns new matrix with all elements set to zero
 */
uint32_t* new_matrix(void) {

    return calloc(g_elements, sizeof(uint32_t));
}

/**
 * Returns new identity matrix
 */
uint32_t* identity_matrix(void) {

    uint32_t* matrix = new_matrix();

    for (ssize_t i = 0; i < g_width; i++) {
        matrix[i * g_width + i] = 1;
    }

    return matrix;
}

/**
 * Returns new matrix with elements generated at random using given seed
 */
uint32_t* random_matrix(uint32_t seed) {

    uint32_t* matrix = new_matrix();

    set_seed(seed);

    for (ssize_t i = 0; i < g_elements; i++) {
        matrix[i] = fast_rand();
    }

    return matrix;
}

/**
 * Returns new matrix with all elements set to given value
 */
uint32_t* uniform_matrix(uint32_t value) {

    uint32_t* matrix = new_matrix();

    for (ssize_t i = 0; i < g_elements; i++) {
        matrix[i] = value;
    }

    return matrix;
}

uint32_t* uniform_matrix_nomem(uint32_t value, uint32_t* matrix) {

    for (ssize_t i = 0; i < g_elements; i++) {
        matrix[i] = value;
    }

    return matrix;
}

/**
 * Returns new matrix with elements in sequence from given start and step
 */
uint32_t* sequence_matrix(uint32_t start, uint32_t step) {

    uint32_t* matrix = new_matrix();
    uint32_t current = start;

    for (ssize_t i = 0; i < g_elements; i++) {
        matrix[i] = current;
        current += step;
    }

    return matrix;
}

////////////////////////////////
///     MATRIX OPERATIONS    ///
////////////////////////////////

/**
 * Returns new matrix with elements cloned from given matrix
 */
uint32_t* cloned(const uint32_t* matrix) {

    uint32_t* result = new_matrix();

    for (ssize_t i = 0; i < g_elements; i++) {
        result[i] = matrix[i];
    }

    return result;
}

/**
 * Returns new matrix with elements ordered in reverse
 */
uint32_t* reversed(const uint32_t* matrix) {

    uint32_t* result = new_matrix();

    for (ssize_t i = 0; i < g_elements; i++) {
        result[i] = matrix[g_elements - 1 - i];
    }

    return result;
}

/**
 * Returns new transposed matrix
 */
uint32_t* transposed(const uint32_t* matrix) {

    uint32_t* result = new_matrix();

    for (ssize_t y = 0; y < g_height; y++) {
        for (ssize_t x = 0; x < g_width; x++) {
            result[x * g_width + y] = matrix[y * g_width + x];
        }
    }

    return result;
}

/**
 * Returns new matrix with scalar added to each element
 */
uint32_t* scalar_add(const uint32_t* matrix, uint32_t scalar) {

    uint32_t* result = new_matrix();
    
    for (ssize_t i = 0; i < g_elements; ++i) {
        result[i] = matrix[i] + scalar;
    }

    return result;
}

/**
 * Returns new matrix with scalar multiplied to each element
 */
uint32_t* scalar_mul(const uint32_t* matrix, uint32_t scalar) {

    uint32_t* result = new_matrix();

    for (ssize_t i = 0; i < g_elements; ++i) {
        result[i] = matrix[i] * scalar;
    }
    
    return result;
}

/**
 * Returns new matrix with elements added at the same index
 */
uint32_t* matrix_add(const uint32_t* matrix_a, const uint32_t* matrix_b) {

    uint32_t* result = new_matrix();

    for (ssize_t i = 0; i < g_elements; ++i) {
        result[i] = matrix_a[i] + matrix_b[i];
    }
    
    return result;
}

/**
 * Assumes result is already zeroed.
 */
void matrix_mul_nomem(const uint32_t* matrix_a, const uint32_t* matrix_b, uint32_t* result) {

    for (ssize_t y = 0; y < g_height; y++) {
        for (ssize_t k = 0; k < g_width; k++) {
            for (ssize_t x = 0; x < g_width; x++) {
                result[y * g_width + x] += matrix_a[y * g_width + k] * matrix_b[k * g_width + x];
            }
        }
    }
}

/**
 * Returns new matrix, multiplying the two matrices together
 */
uint32_t* matrix_mul(const uint32_t* matrix_a, const uint32_t* matrix_b) {

    uint32_t* result = new_matrix();
    
    //TODO: another matrix_mul_nomem that doesn't 0 elements as it goes.
    matrix_mul_nomem(matrix_a, matrix_b, result);  

    return result;
}

void strass_add(const uint32_t* a, ssize_t aw, ssize_t ah, 
                 const uint32_t* b, ssize_t bw, ssize_t bh, 
                 uint32_t* c) {
    ssize_t max_w = (aw > bw ? aw : bw);
    ssize_t max_h = (ah > bh ? ah : bh);
    
    for (ssize_t i = 0; i < max_h - 1; ++i) {
        for (ssize_t j = 0; j < max_w - 1; ++j) {
           c[max_w*i + j] = a[max_w*i + j] + b[max_w*i + j];
        }
    }

    if (aw > bw) {
        for (ssize_t i = 1; i < max_h - 1; ++i) {
            c[max_w*i - 1] = a[max_w*i - 1];
        }
    } else if (aw < bw) {
        for (ssize_t i = 1; i < max_h - 1; ++i) {
            c[max_w*i - 1] = b[max_w*i - 1];
        }
    } else {
        for (ssize_t i = 1; i < max_h - 1; ++i) {
            c[max_w*i - 1] = a[max_w*i - 1] + b[max_w*i - 1];
        }
    }

    if (ah > bh) {
        for (ssize_t j = 0; j < max_w - 1; ++j) {
            c[max_h - 1 + j] = a[max_h - 1 + j];
        }
    } else if (ah < bh) {
        for (ssize_t j = 0; j < max_w - 1; ++j) {
            c[max_h - 1 + j] = b[max_h - 1 + j];
        }
    } else {
        for (ssize_t j = 1; j < max_w - 1; ++j) {
            c[max_h - 1 + j] = a[max_h - 1 + j] + b[max_h - 1 + j];
        }
    }

    if (ah == bh && aw == bw) {
        c[max_w*max_h - 1] = a[max_w*max_h - 1] + b[max_w*max_h - 1];
    } 
    
}

void strass_sub(const uint32_t* a, ssize_t aw, ssize_t ah, 
                 const uint32_t* b, ssize_t bw, ssize_t bh, 
                 uint32_t* c) {
    ssize_t max_w = (aw > bw ? aw : bw);
    ssize_t max_h = (ah > bh ? ah : bh);
    
    for (ssize_t i = 0; i < max_h - 1; ++i) {
        for (ssize_t j = 0; j < max_w - 1; ++j) {
           c[max_w*i + j] = a[max_w*i + j] - b[max_w*i + j];
        }
    }

    if (aw > bw) {
        for (ssize_t i = 1; i < max_h - 1; ++i) {
            c[max_w*i - 1] = a[max_w*i - 1];
        }
    } else if (aw < bw) {
        for (ssize_t i = 1; i < max_h - 1; ++i) {
            c[max_w*i - 1] = 0 - b[max_w*i - 1];
        }
    } else {
        for (ssize_t i = 1; i < max_h - 1; ++i) {
            c[max_w*i - 1] = a[max_w*i - 1] - b[max_w*i - 1];
        }
    }

    if (ah > bh) {
        for (ssize_t j = 0; j < max_w - 1; ++j) {
            c[max_h - 1 + j] = a[max_h - 1 + j];
        }
    } else if (ah < bh) {
        for (ssize_t j = 0; j < max_w - 1; ++j) {
            c[max_h - 1 + j] = 0 - b[max_h - 1 + j];
        }
    } else {
        for (ssize_t j = 1; j < max_w - 1; ++j) {
            c[max_h - 1 + j] = a[max_h - 1 + j] - b[max_h - 1 + j];
        }
    }

    if (ah == bh && aw == bw) {
        c[max_w*max_h - 1] = a[max_w*max_h - 1] - b[max_w*max_h - 1];
    } 
    
}

void strass_mul(const uint32_t* a, ssize_t aw, ssize_t ah, 
                 const uint32_t* b, ssize_t bw, ssize_t bh, 
                 uint32_t* c) {

        ssize_t min_p = (aw < bh ? aw : bh);
        ssize_t max_w = (aw < bw ? aw : bw);

        for (ssize_t y = 0; y < ah; y++) {
            for (ssize_t k = 0; k < min_p; k++) {
                for (ssize_t x = 0; x < bw; x++) {
                    c[y * max_w + x] += a[y * g_width + k] * b[k * g_width + x];
                }
            }
        }
}

void strassen(const uint32_t* a, ssize_t aw, ssize_t ah, 
              const uint32_t* b, ssize_t bw, ssize_t bh, 
              uint32_t* c) {
    
    
    // Get largest dimension
    ssize_t n = (aw > bw ? aw : bw);
    ssize_t m = (ah > bh ? ah : bh);
    n = (n > m ? n : m);
    
    if (n <= STRASS_THRESH) {
        strass_mul(a, aw, ah, b, bw, bh, c);
        return;
    } else {
    
        m = (n&1 ? (n+1)/2 : n/2); 
        n = m*m;

        // Pointers to each quadrant of input matrices
        
        ssize_t a11w = aw&1 ? (aw+1)/2 : aw/2;
        ssize_t a11h = ah&1 ? (ah+1)/2 : ah/2;
        ssize_t a22w = aw-a11w;
        ssize_t a22h = ah-a11h;
        ssize_t b11w = bw&1 ? (bw+1)/2 : bw/2;
        ssize_t b11h = bh&1 ? (bh+1)/2 : bh/2;
        ssize_t b22w = bw-b11w;
        ssize-t b22h = bh-b11h;

        uint32_t* a11 = a;
        uint32_t* a12 = a + a11w;
        uint32_t* a21 = a + g_width*a11h;
        uint32_t* a22 = a21 + a11w;
        uint32_t* b11 = b;
        uint32_t* b12 = b + b11w;
        uint32_t* b21 = b + g_width*b11h;
        uint32_t* b22 = b21 + b11w;
        
        //TODO: consider zeroing fewer of these, when multiplications happen, instead of calloc
        uint32_t* S = (uint32_t*)calloc(11, sizeof(uint32_t) * n);
        strass_sub(b12, b22w, b11h, b22, b22w, b22h, S + n);
        strass_add(a11, a11w, a11h, a12, a22w, a11h, S + 2*n);
        strass_sub(b21, b11w, b22h, b11, b11w, b11h, S + 3*n);
        strass_add(a11, a11w, a11h, a22, a22w, a22h, S + 4*n);
        strass_add(b11, b11w, b11h, b22, b22w, b22h, S + 5*n);
        strass_sub(a12, a22w, a11h, a22, a22w, a22h, S + 6*n);
        strass_add(b21, b11w, b22h, b22, b22w, b22h, S + 7*n);
        strass_sub(a11, a11w, a11h, a21, a11w, a22h, S + 8*n);
        strass_add(b11, b11w, b11h, b12, b22w, b11h, S + 9*n);
        strass_add(a21, a11w, a22h, a22, a22w, a22h, S + 10*n);
        
        //P1 = A11*S1
        strassen(a11, a11w, a11h, S + n, m, m, S);
        //P2 = S2*B22
        strassen(S + 2*n, m, m, b22, b22w, b22h, S + n);
        //P3 = S3*B11
        strassen(S + 3*n, m, m, b11, b11w, b11h, S + 2*n);
        //P4 = A22*S4
        strassen(a22, a22w, a22h, S + 4*n, m, m, S + 3*n);
        //P5 = S5*S6
        strassen(S + 5*n, m, m, S + 6*n, m, m, S + 4*n);
        //P6 = S7*S8
        strassen(S + 7*n, m, m, S + 8*n, m, m, S + 5*n);
        //P7 = S9*S10
        strassen(S + 9*n, m, m, S + 10*n, m, m, S + 6*n);
        
        

        free(S);
    }
}


/**
 * Returns new matrix, powering the matrix to the exponent
 */
uint32_t* matrix_pow(const uint32_t* matrix, uint32_t exponent) {
    
    
    uint32_t* result = identity_matrix();
    
    if (exponent == 0) return result;
    int width = 32 -__builtin_clz(exponent);

    uint32_t* mpowers = (uint32_t *)calloc(width, sizeof(uint32_t)*g_elements);

    for (int i = 0; i < width; ++i) {
        if (i == 0) {
            memcpy(mpowers, matrix, g_elements*sizeof(uint32_t));
        } else {
            matrix_mul_nomem((mpowers + g_elements*(i-1)), (mpowers + g_elements*(i-1)), (mpowers + g_elements*i));
        }
    }
    
    uint32_t* temp = new_matrix();

    for (uint32_t i = 0; i < width; ++i) {
        if (exponent & (1 << i)) {
            matrix_mul_nomem(result, (mpowers + g_elements*i), temp);
            
            uint32_t* swap = result;
            result = temp;
            temp = swap;
        }
    }
    free(mpowers);
    free(temp);
    

    // TODO: diagonalisation

    return result;
}

////////////////////////////////
///       COMPUTATIONS       ///
////////////////////////////////

/**
 * Returns the sum of all elements
 */
uint32_t get_sum(const uint32_t* matrix) {

    /*
        to do

        1 2
        2 1 => 6

        1 1
        1 1 => 4
    */
    uint32_t sum = 0;

    for (uint32_t i = 0; i < g_elements; ++i) {
        sum += matrix[i];
    }


    return sum;
}

/**
 * Returns the trace of the matrix
 */
uint32_t get_trace(const uint32_t* matrix) {

    /*
        to do

        1 0
        0 1 => 2

        2 1
        1 2 => 4
    */

    uint32_t trace = 0;

    for (uint32_t i = 0; i < g_width; ++i) {
        trace += matrix[i * g_width + i];
    }

    return trace;
}

/**
 * Returns the smallest value in the matrix
 */
uint32_t get_minimum(const uint32_t* matrix) {

    uint32_t min = matrix[0];

    for (uint32_t i = 1; i < g_elements; ++i) {
        if (matrix[i] < min) {
            min = matrix[i];
        }
    } 
    
    /*
        to do

        1 2
        3 4 => 1

        4 3
        2 1 => 1
    */

    return min;
}

/**
 * Returns the largest value in the matrix
 */
uint32_t get_maximum(const uint32_t* matrix) {

    uint32_t max = matrix[0];

    for (uint32_t i = 1; i < g_elements; ++i) {
        if (matrix[i] > max) {
            max = matrix[i];
        }
    }
    /*
        to do

        1 2
        3 4 => 4

        4 3
        2 1 => 4
    */

    return max;
}

/**
 * Returns the frequency of the value in the matrix
 */
uint32_t get_frequency(const uint32_t* matrix, uint32_t value) {
    
    uint32_t freq = 0;

    for (uint32_t i = 0; i < g_elements; ++i) {
        if (matrix[i] == value) ++freq;
    }

    /*
        to do

        1 1
        1 1 :: 1 => 4

        1 0
        0 1 :: 2 => 0
    */

    return freq;
}
