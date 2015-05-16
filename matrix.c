#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <stdbool.h>
#include <pthread.h>
#include <inttypes.h>

#include <string.h>
//#include <assert.h>

#include "matrix.h"

#define STRASS_THRESH 100
#define STRASS_MIN 5

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
            printf("%" PRIu32, matrix[y * g_width + x]);
        }

        printf("\n");
    }
}

void strass_display(const uint32_t* matrix, ssize_t mw, ssize_t mh, ssize_t ms) {

    for (ssize_t y = 0; y < mh; y++) {
        for (ssize_t x = 0; x < mw; x++) {
            if (x > 0) printf(" ");
            printf("%" PRIu32 " ", matrix[y * ms + x]);
        }

        printf("\n");
    }
    printf("\n");
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
    
    strassen(matrix_a, g_width, g_height, g_width, matrix_b, g_width, g_height, g_width, result, g_width);

    return result;
}

/**
 * Adds a and b, placing the result in c.
 * Input matrix dimensions may differ by at most 2 per side.
 * In the case where one matrix is smaller than the other in one dimension,
 * that matrix is implicitly padded with zeroes on that edge.
 *
 * c is assumed to be exactly large enough to hold the resulting sum.
 *
 * a, b, c may all be submatrices, being part of a larger parent in memory.
 * Hence, while aw, ah, bw, bh define the dimensions of a and b;
 * as, bs, cs are the stride required to get to the corresponding next row.
 * That is, it is the width of the parent matrix.
 */
void strass_add(const uint32_t* a, ssize_t aw, ssize_t ah, ssize_t as,
                 const uint32_t* b, ssize_t bw, ssize_t bh, ssize_t bs,
                 uint32_t* c, ssize_t cs) {
    
    for (ssize_t y = 0; y < ah; ++y) {
        for (ssize_t x = 0; x < aw; ++x) {
            c[cs*y + x] = a[as*y + x];
        }
    }

    for (ssize_t y = 0; y < bh; ++y) {
        for (ssize_t x = 0; x < bw; ++x) {
            c[cs*y + x] += b[bs*y + x];
        }
    }
    /*ssize_t min_w = (aw < bw ? aw : bw);
    ssize_t min_h = (ah < bh ? ah : bh);
    ssize_t max_w = (aw > bw ? aw : bw);
    ssize_t max_h = (ah > bh ? ah : bh);
    ssize_t w_off = max_w - 1;
    ssize_t h_off = max_h - 1;
    
    // Add the parts which definitely exist, the common largest submatrices.
    for (ssize_t i = 0; i < min_h; ++i) {
        for (ssize_t j = 0; j < min_w; ++j) {
           c[cs*i + j] = a[as*i + j] + b[bs*i + j];
        }
    }
    
    // Add the last columns.
    if (aw > bw) {
        for (ssize_t i = 0; i < ah; ++i) {
            c[cs*i + w_off] = a[as*i + w_off];
        }
        
        if ((max_w - min_w) > 1) {
            ssize_t w_off2 = w_off - 1;
            for (ssize_t i = 0; i < ah; ++i) {
                c[cs*i + w_off2] = a[as*i + w_off2];
            }   
        }
    } else if (aw < bw) {
        for (ssize_t i = 0; i < bh; ++i) {
            c[cs*i + w_off] = b[bs*i + w_off];
        }
        
        if ((max_w - min_w) > 1) {
            ssize_t w_off2 = w_off - 1;
            for (ssize_t i = 0; i < bh; ++i) {
                c[cs*i + w_off2] = b[bs*i + w_off2];
            }   
        }
    }
    
    // Add the last row.
    if (ah > bh) {
        ssize_t cv = cs*h_off;
        ssize_t av = as*h_off;
        for (ssize_t j = 0; j < aw; ++j) {
            c[cv + j] = a[av + j];
        }
        
        if ((max_h - min_h) > 1) {
            ssize_t cv2 = cv - cs;
            ssize_t av2 = av - as;
            for (ssize_t j = 0; j < aw; ++j) {
                c[cv2 + j] = a[av2 + j];
            }
        }
    } else if (ah < bh) {
        ssize_t cv = cs*h_off;
        ssize_t bv = bs*h_off;
        for (ssize_t j = 0; j < bw; ++j) {
            c[cv + j] = b[bv + j];
        }

        if ((max_h - min_h) > 1) {
            ssize_t cv2 = cv - cs;
            ssize_t bv2 = bv - bs;
            for (ssize_t j = 0; j < bw; ++j) {
                c[cv2 + j] = b[bv2 + j];
            }
        }
    }*/
}

/**
 * Subtracts b from a, placing the result in c.
 * Input matrix dimensions may differ by at most 1 per side.
 *
 * Details are as strass_add.
 */ 
void strass_sub(const uint32_t* a, ssize_t aw, ssize_t ah, ssize_t as,
                 const uint32_t* b, ssize_t bw, ssize_t bh, ssize_t bs,
                 uint32_t* c, ssize_t cs) {
    
    for (ssize_t y = 0; y < ah; ++y) {
        for (ssize_t x = 0; x < aw; ++x) {
            c[cs*y + x] = a[as*y + x];
        }
    }

    for (ssize_t y = 0; y < bh; ++y) {
        for (ssize_t x = 0; x < bw; ++x) {
            c[cs*y + x] -= b[bs*y + x];
        }
    }
    
    
    /*ssize_t min_w = (aw < bw ? aw : bw);
    ssize_t min_h = (ah < bh ? ah : bh);
    ssize_t max_w = (aw > bw ? aw : bw);
    ssize_t max_h = (ah > bh ? ah : bh);
    ssize_t w_off = max_w - 1;
    ssize_t h_off = max_h - 1;
    
    // Subtract the parts which definitely exist, the common largest submatrices.
    for (ssize_t i = 0; i < min_h; ++i) {
        for (ssize_t j = 0; j < min_w; ++j) {
           c[cs*i + j] = a[as*i + j] - b[bs*i + j];
        }
    }
    
    // Handle the last column.
    if (aw > bw) {
        for (ssize_t i = 0; i < ah; ++i) {
            c[cs*i + w_off] = a[as*i + w_off];
        }
        
        if ((max_w - min_w) > 1) {
            ssize_t w_off2 = w_off - 1;
            for (ssize_t i = 0; i < ah; ++i) {
                c[cs*i + w_off2] = a[as*i + w_off2];
            }   
        }
    } else if (aw < bw) {
        for (ssize_t i = 0; i < bh; ++i) {
            c[cs*i + w_off] = 0 - b[bs*i + w_off];
        }

        if ((max_w - min_w) > 1) {
            ssize_t w_off2 = w_off - 1;
            for (ssize_t i = 0; i < bh; ++i) {
                c[cs*i + w_off2] = 0 - b[bs*i + w_off];
            }   
        }
    }
    
    // Handle the last row.
    if (ah > bh) {
        ssize_t cv = cs*h_off;
        ssize_t av = as*h_off;
        for (ssize_t j = 0; j < aw; ++j) {
            c[cv + j] = a[av + j];
        }

        if ((max_h - min_h) > 1) {
            ssize_t cv2 = cv - cs;
            ssize_t av2 = av - as;
            for (ssize_t j = 0; j < aw; ++j) {
                c[cv2 + j] = a[av2 + j];
            }
        }
    } else if (ah < bh) {
        ssize_t cv = cs*h_off;
        ssize_t bv = bs*h_off;
        for (ssize_t j = 0; j < bw; ++j) {
            c[cv + j] = 0 - b[bv + j];
        }

        if ((max_h - min_h) > 1) {
            ssize_t cv2 = cv - cs;
            ssize_t bv2 = bv - bs;
            for (ssize_t j = 0; j < bw; ++j) {
                c[cv2 + j] = 0 - b[bv2 + j];
            }
        }
    }*/
}

/**
 * Base case multiplication,
 * Multiplies a by b, placing the result in c.
 * Input matrix dimensions may differ by at most 1 per side.
 *
 * Details are as strass_add.
 */
void strass_mul(const uint32_t* a, ssize_t aw, ssize_t ah, ssize_t as,
                 const uint32_t* b, ssize_t bw, ssize_t bh, ssize_t bs,
                 uint32_t* c, ssize_t cs) {

    ssize_t min_p = (aw < bh ? aw : bh);

    for (ssize_t y = 0; y < ah; y++) {
        for (ssize_t k = 0; k < min_p; k++) {
            for (ssize_t x = 0; x < bw; x++) {
                c[y*cs + x] += a[y*as + k] * b[k*bs + x];
            }
        }
    }
}

void strassen(const uint32_t* a, ssize_t aw, ssize_t ah, ssize_t as,
              const uint32_t* b, ssize_t bw, ssize_t bh, ssize_t bs,
              uint32_t* c, ssize_t cs) {
    //assert(abs(aw - bw) <= 1);
    //assert(abs(ah - bh) <= 1);
    
    if (aw == 0 || ah == 0 || bw == 0 || bh == 0) {
        return;
    }

    // Get largest dimension
    ssize_t max_w = (aw > bw ? aw : bw);
    ssize_t max_h = (ah > bh ? ah : bh);
    ssize_t m = (max_w > max_h ? max_w : max_h);
    //printf("m : %zd", m);
    
    if (m <= STRASS_THRESH || aw <= STRASS_MIN || bw <= STRASS_MIN) {
        
        //strass_display(a, aw, ah, as);
        //strass_display(b, bw, bh, bs);

        strass_mul(a, aw, ah, as, b, bw, bh, bs, c, cs);
        return;
    } else {
         
        //printf("aw: %zd ah: %zd  bw: %zd bh: %zd\n", aw, ah, bw, bh);
        // Pointers to each quadrant of input/output matrices
        /*if (aw == 3 && bw == 1) {
            strass_display(a, aw, ah, as);
            strass_display(b, bw, bh, bs);
        }*/

        ssize_t a11w = aw&1 ? (aw+1)/2 : aw/2;
        ssize_t a11h = ah&1 ? (ah+1)/2 : ah/2;
        ssize_t a22w = aw-a11w;
        ssize_t a22h = ah-a11h;
        
        ssize_t b11w = bw&1 ? (bw+1)/2 : bw/2;
        ssize_t b11h = bh&1 ? (bh+1)/2 : bh/2;
        ssize_t b22w = bw-b11w;
        ssize_t b22h = bh-b11h;
        
        ssize_t c11w = (max_w&1 ? (max_w+1)/2 : max_w/2);
        ssize_t c11h = (max_h&1 ? (max_h+1)/2 : max_h/2); 
        ssize_t c22w = max_w-c11w;
        ssize_t c22h = max_h-c11h;


        if (aw > bw) {
            b11w = a11w;
            b22w = bw-b11w;    
        } else if (aw < bw) {
            a11w = b11w;
            a22w = aw-a11w;
        } else {

        }

        if (ah > bh) {
            b11h = a11h;
            b22h = bh-b11h;
        } else if (ah < bh) {
            a11h = b11h;
            a22h = ah-a11h;
        } else {

        }

        const uint32_t* a11 = a;
        const uint32_t* a12 = a + a11w;
        const uint32_t* a21 = a + as*a11h;
        const uint32_t* a22 = a21 + a11w;

        const uint32_t* b11 = b;
        const uint32_t* b12 = b + b11w;
        const uint32_t* b21 = b + bs*b11h;
        const uint32_t* b22 = b21 + b11w;

        uint32_t* c11 = c;
        uint32_t* c12 = c + c11w;
        uint32_t* c21 = c + cs*c11h;
        uint32_t* c22 = c21 + c11w;

        //strass_display(a11, a11w, a11h, as);
        //strass_display(a12, a22w, a11h, as);
        //strass_display(a21, a11w, a22h, as);
        //strass_display(a22, a22w, a22h, as);

        //strass_display(b11, b11w, b11h, bs);
        //strass_display(b12, b22w, b11h, bs);
        //strass_display(b21, b11w, b22h, bs);
        //strass_display(b22, b22w, b22h, bs);

        m = (c11w > c11h ? c11w : c11h);
        ssize_t n = m*m;
        uint32_t* S = (uint32_t*)calloc(11, sizeof(uint32_t) * n);
        
        strass_sub(b12, b22w, b11h, bs, b22, b22w, b22h, bs, S + n, m);    // S1  = b12 - b22
        strass_add(a11, a11w, a11h, as, a12, a22w, a11h, as, S + 2*n, m);  // S2  = a11 + a12 
        strass_add(a21, a11w, a22h, as, a22, a22w, a22h, as, S + 3*n, m);  // S3  = a21 + a22
        strass_sub(b21, b11w, b22h, bs, b11, b11w, b11h, bs, S + 4*n, m);  // S4  = b21 - b11
        strass_add(a11, a11w, a11h, as, a22, a22w, a22h, as, S + 5*n, m);  // S5  = a11 + a22
        strass_add(b11, b11w, b11h, bs, b22, b22w, b22h, bs, S + 6*n, m);  // S6  = b11 + b22
        strass_sub(a12, a22w, a11h, as, a22, a22w, a22h, as, S + 7*n, m);  // S7  = a12 - a22
        strass_add(b21, b11w, b22h, bs, b22, b22w, b22h, bs, S + 8*n, m);  // S8  = b21 + b22
        strass_sub(a11, a11w, a11h, as, a21, a11w, a22h, as, S + 9*n, m);  // S9  = a11 - a21
        strass_add(b11, b11w, b11h, bs, b12, b22w, b11h, bs, S + 10*n, m); // S10 = b11 + b12
        


        //printf("\nEsses\n");
        //strass_display(S+n, m, m, m);
        //strass_display(S+2*n, m, m, m);
        //strass_display(S+3*n, m, m, m);
        //strass_display(S+4*n, m, m, m);
        //strass_display(S+5*n, m, m, m);
        //strass_display(S+6*n, m, m, m);
        //strass_display(S+7*n, m, m, m);
        //strass_display(S+8*n, m, m, m);
        //strass_display(S+9*n, m, m, m);
        //strass_display(S+10*n, m, m, m);


        //printf("\nPees\n");


        //P1 = A11*S1
        memset(S, 0, n*sizeof(uint32_t));
        strassen(a11, a11w, a11h, as, S + n, m, m, m, S, m);
        //P2 = S2*B22
        memset(S + n, 0, n*sizeof(uint32_t));
        strassen(S + 2*n, m, m, m, b22, b22w, b22h, bs, S + n, m);
        //P3 = S3*B11
        memset(S + 2*n, 0, n*sizeof(uint32_t));
        strassen(S + 3*n, m, m, m, b11, b11w, b11h, bs, S + 2*n, m);
        //P4 = A22*S4
        memset(S + 3*n, 0, n*sizeof(uint32_t));
        strassen(a22, a22w, a22h, as, S + 4*n, m, m, m, S + 3*n, m);
        //P5 = S5*S6
        memset(S + 4*n, 0, n*sizeof(uint32_t));
        strassen(S + 5*n, m, m, m, S + 6*n, m, m, m, S + 4*n, m);
        //P6 = S7*S8
        memset(S + 5*n, 0, 2*n*sizeof(uint32_t));
        strassen(S + 7*n, m, m, m, S + 8*n, m, m, m, S + 5*n, m);
        //P7 = S9*S10
        strassen(S + 9*n, m, m, m, S + 10*n, m, m, m, S + 6*n, m);
        
        //strass_display(S, m, m, m);
        //strass_display(S+n, m, m, m);
        //strass_display(S+2*n, m, m, m);
        //strass_display(S+3*n, m, m, m);
        //strass_display(S+4*n, m, m, m);
        //strass_display(S+5*n, m, m, m);
        //strass_display(S+6*n, m, m, m);
        
        //printf("successful for submults.\n");
        // C11 = P5 + P4 - P2 + P6
        //memset(S + 9*n, 0, 2*n*sizeof(uint32_t));
        strass_add(S + 4*n, m, m, m, S + 3*n, m, m, m, S + 10*n, m);
        strass_sub(S + 5*n, m, m, m, S + n, m, m, m, S + 9*n, m);
        strass_add(S + 10*n, c11w, c11h, m, S + 9*n, c11w, c11h, m, c11, cs);
        //printf("C11\n");
        // C12 = P1 + P2
        strass_add(S, c22w, c11h, m, S + n, c22w, c11h, m, c12, cs);
        //printf("C12\n");
        // C21 = P3 + P4
        strass_add(S + 2*n, c11w, c22h, m, S + 3*n, c11w, c22h, m, c21, cs);
        //printf("C21\n");
        // C22 = P1 + P5 - P3 - P7
        //memset(S + n, 0, n*sizeof(uint32_t));
        //memset(S + 3*n, 0, n*sizeof(uint32_t));
        strass_add(S, m, m, m, S + 4*n, m, m, m, S + n, m);
        strass_add(S + 2*n, m, m, m, S + 6*n, m, m, m, S + 3*n, m);
        strass_sub(S + n, c22w, c22h, m, S + 3*n, c22w, c22h, m, c22, cs);
        //printf("C22\n");
        
        //printf("\nCees\n"); 
        //strass_display(c11, c11w, c11h, cs);
        //strass_display(c12, c22w, c11h, cs);
        //strass_display(c21, c11w, c22h, cs);
        //strass_display(c22, c22w, c22h, cs);

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
                     (mpowers + g_elements*i), g_width);
        }
    }
    
    uint32_t* temp = new_matrix();

    for (uint32_t i = 0; i < width; ++i) {
        if (exponent & (1 << i)) {
            //matrix_mul_nomem(result, (mpowers + g_elements*i), temp);
            strassen(result, g_width, g_height, g_width,
                     (mpowers + g_elements*i), g_width, g_height, g_width,
                     temp, g_width);

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
