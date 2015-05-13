#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <stdbool.h>
#include <pthread.h>
#include <inttypes.h>

#include <string.h>

#include "matrix.h"

#define STRASS_THRESH 100

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
    
    //matrix_mul_nomem(matrix_a, matrix_b, result);  
    
    strassen(matrix_a, g_width, g_height, g_width, matrix_b, g_width, g_height, g_width, result);

    return result;
}

void strass_add(const uint32_t* a, ssize_t aw, ssize_t ah, ssize_t as,
                 const uint32_t* b, ssize_t bw, ssize_t bh, ssize_t bs,
                 uint32_t* c, ssize_t cs) {
    ssize_t max_w = (aw > bw ? aw : bw);
    ssize_t max_h = (ah > bh ? ah : bh);
    
    if (cs == -1) cs = max_w;

    for (ssize_t i = 0; i < max_h - 1; ++i) {
        for (ssize_t j = 0; j < max_w - 1; ++j) {
           c[cs*i + j] = a[as*i + j] + b[bs*i + j];
        }
    }

    if (aw > bw) {
        for (ssize_t i = 1; i < max_h - 1; ++i) {
            c[cs*i - 1] = a[as*i - 1];
        }
    } else if (aw < bw) {
        for (ssize_t i = 1; i < max_h - 1; ++i) {
            c[cs*i - 1] = b[bs*i - 1];
        }
    } else {
        for (ssize_t i = 1; i < max_h - 1; ++i) {
            c[cs*i - 1] = a[as*i - 1] + b[bs*i - 1];
        }
    }

    if (ah > bh) {
        for (ssize_t j = 0; j < max_w - 1; ++j) {
            c[cs*(max_h - 1) + j] = a[as*(max_h - 1) + j];
        }
    } else if (ah < bh) {
        for (ssize_t j = 0; j < max_w - 1; ++j) {
            c[cs*(max_h - 1) + j] = b[bs*(max_h - 1) + j];
        }
    } else {
        for (ssize_t j = 1; j < max_w - 1; ++j) {
            c[cs*(max_h - 1) + j] = a[as*(max_h - 1) + j] + b[bs*(max_h - 1) + j];
        }
    }

    if (ah == bh && aw == bw) {
        c[cs*max_h - 1] = a[max_w*max_h - 1] + b[max_w*max_h - 1];
    } 
    
}

void strass_sub(const uint32_t* a, ssize_t aw, ssize_t ah, ssize_t as,
                 const uint32_t* b, ssize_t bw, ssize_t bh, ssize_t bs,
                 uint32_t* c, ssize_t cs) {
    ssize_t max_w = (aw > bw ? aw : bw);
    ssize_t max_h = (ah > bh ? ah : bh);
    
    if (cs == -1) cs = max_w;

    for (ssize_t i = 0; i < max_h - 1; ++i) {
        for (ssize_t j = 0; j < max_w - 1; ++j) {
           c[cs*i + j] = a[as*i + j] - b[bs*i + j];
        }
    }

    if (aw > bw) {
        for (ssize_t i = 1; i < max_h - 1; ++i) {
            c[cs*i - 1] = 0 - a[as*i - 1];
        }
    } else if (aw < bw) {
        for (ssize_t i = 1; i < max_h - 1; ++i) {
            c[cs*i - 1] = 0 - b[bs*i - 1];
        }
    } else {
        for (ssize_t i = 1; i < max_h - 1; ++i) {
            c[cs*i - 1] = a[as*i - 1] - b[bs*i - 1];
        }
    }

    if (ah > bh) {
        for (ssize_t j = 0; j < max_w - 1; ++j) {
            c[cs*(max_h - 1) + j] = 0 - a[as*(max_h - 1) + j];
        }
    } else if (ah < bh) {
        for (ssize_t j = 0; j < max_w - 1; ++j) {
            c[cs*(max_h - 1) + j] = 0 - b[bs*(max_h - 1) + j];
        }
    } else {
        for (ssize_t j = 1; j < max_w - 1; ++j) {
            c[cs*(max_h - 1) + j] = a[as*(max_h - 1) + j] - b[bs*(max_h - 1) + j];
        }
    }

    if (ah == bh && aw == bw) {
        c[cs*max_h - 1] = a[max_w*max_h - 1] - b[max_w*max_h - 1];
    } 
    
}

void strass_mul(const uint32_t* a, ssize_t aw, ssize_t ah, ssize_t as,
                 const uint32_t* b, ssize_t bw, ssize_t bh, ssize_t bs,
                 uint32_t* c) {

    ssize_t min_p = (aw < bh ? aw : bh);
    ssize_t max_w = (aw < bw ? aw : bw);

    for (ssize_t y = 0; y < ah; y++) {
        for (ssize_t k = 0; k < min_p; k++) {
            for (ssize_t x = 0; x < bw; x++) {
                c[y * max_w + x] += a[y * as + k] * b[k * bs + x];
            }
        }
    }
}

void strassen(const uint32_t* a, ssize_t aw, ssize_t ah, ssize_t as,
              const uint32_t* b, ssize_t bw, ssize_t bh, ssize_t bs,
              uint32_t* c) {
    
    // Get largest dimension
    ssize_t max_w = (aw > bw ? aw : bw);
    ssize_t max_h = (ah > bh ? ah : bh);
    ssize_t n = (max_w > max_h ? max_w : max_h);
    
    if (n <= STRASS_THRESH) {
        strass_mul(a, aw, ah, as, b, bw, bh, bs, c);
        return;
    } else {
    
        ssize_t m = (n&1 ? (n+1)/2 : n/2); 
        n = m*m;

        // Pointers to each quadrant of input/output matrices
        
        ssize_t a11w = aw&1 ? (aw+1)/2 : aw/2;
        ssize_t a11h = ah&1 ? (ah+1)/2 : ah/2;
        ssize_t a22w = aw-a11w;
        ssize_t a22h = ah-a11h;
        ssize_t b11w = bw&1 ? (bw+1)/2 : bw/2;
        ssize_t b11h = bh&1 ? (bh+1)/2 : bh/2;
        ssize_t b22w = bw-b11w;
        ssize_t b22h = bh-b11h;

        const uint32_t* a11 = a;
        const uint32_t* a12 = a + a11w;
        const uint32_t* a21 = a + as*a11h;
        const uint32_t* a22 = a21 + a11w;
        
        const uint32_t* b11 = b;
        const uint32_t* b12 = b + b11w;
        const uint32_t* b21 = b + bs*b11h;
        const uint32_t* b22 = b21 + b11w;

        uint32_t* c11 = c;
        uint32_t* c12 = c + m;
        uint32_t* c21 = c + max_w*m;
        uint32_t* c22 = c21 + m;

        
        //TODO: consider zeroing fewer of these, when multiplications happen, instead of calloc
        uint32_t* S = (uint32_t*)calloc(11, sizeof(uint32_t) * n);
        strass_sub(b12, b22w, b11h, bs, b22, b22w, b22h, bs, S + n, -1);    // S1
        strass_add(a11, a11w, a11h, as, a12, a22w, a11h, as, S + 2*n, -1);  // S2
        strass_sub(b21, b11w, b22h, bs, b11, b11w, b11h, bs, S + 3*n, -1);  // S3
        strass_add(a11, a11w, a11h, as, a22, a22w, a22h, as, S + 4*n, -1);  // S4
        strass_add(b11, b11w, b11h, bs, b22, b22w, b22h, bs, S + 5*n, -1);  // S5
        strass_sub(a12, a22w, a11h, as, a22, a22w, a22h, as, S + 6*n, -1);  // S6
        strass_add(b21, b11w, b22h, bs, b22, b22w, b22h, bs, S + 7*n, -1);  // S7
        strass_sub(a11, a11w, a11h, as, a21, a11w, a22h, as, S + 8*n, -1);  // S8
        strass_add(b11, b11w, b11h, bs, b12, b22w, b11h, bs, S + 9*n, -1);  // S9
        strass_add(a21, a11w, a22h, as, a22, a22w, a22h, as, S + 10*n, -1); // S10
        
        //P1 = A11*S1
        strassen(a11, a11w, a11h, as, S + n, m, m, m, S);
        //P2 = S2*B22
        strassen(S + 2*n, m, m, m, b22, b22w, b22h, bs, S + n);
        //P3 = S3*B11
        strassen(S + 3*n, m, m, m, b11, b11w, b11h, bs, S + 2*n);
        //P4 = A22*S4
        strassen(a22, a22w, a22h, as, S + 4*n, m, m, m, S + 3*n);
        //P5 = S5*S6
        strassen(S + 5*n, m, m, m, S + 6*n, m, m, m, S + 4*n);
        //P6 = S7*S8
        strassen(S + 7*n, m, m, m, S + 8*n, m, m, m, S + 5*n);
        //P7 = S9*S10
        strassen(S + 9*n, m, m, m, S + 10*n, m, m, m, S + 6*n);
        
        // C11 = P5 + P4 - P2 + P6
        strass_add(S + 4*n, m, m, m, S + 3*n, m, m, m, S + 10*n, -1);
        strass_sub(S + 5*n, m, m, m, S + n, m, m, m, S + 9*n, -1);
        strass_add(S + 10*n, m, m, m, S + 9*n, m, m, m, c11, max_w);
        // C12 = P1 + P2
        strass_add(S, m, m, m, S + n, m, m, m, c12, max_w);
        // C21 = P3 + P4
        strass_add(S + 2*n, m, m, m, S + 3*n, m, m, m, c21, max_w);
        // C22 = P1 + P5 - P3 - P7
        strass_add(S, m, m, m, S + 4*n, m, m, m, S + n, -1);
        strass_add(S + 2*n, m, m, m, S + 6*n, m, m, m, S + 3*n, -1);
        strass_sub(S + n, m, m, m, S + 3*n, m, m, m, c22, max_w);
        
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
            //matrix_mul_nomem((mpowers + g_elements*(i-1)), (mpowers + g_elements*(i-1)), (mpowers + g_elements*i));
            strassen((mpowers + g_elements*(i-1)), g_width, g_height, g_width,
                     (mpowers + g_elements*(i-1)), g_width, g_height, g_width,
                     (mpowers + g_elements*i));
        }
    }
    
    uint32_t* temp = new_matrix();

    for (uint32_t i = 0; i < width; ++i) {
        if (exponent & (1 << i)) {
            //matrix_mul_nomem(result, (mpowers + g_elements*i), temp);
            strassen(result, g_width, g_height, g_width,
                     (mpowers + g_elements*i), g_width, g_height, g_width,
                     temp);

            uint32_t* swap = result;
            result = temp;
            temp = swap;
            uniform_matrix_nomem(0, temp);
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

    return freq;
}
