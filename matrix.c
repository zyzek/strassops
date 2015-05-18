/* Parallelisation strategy:
 *
 * 1. Spool up N threads, 
 * 2. Have global current-op
 * 3. Have global struct for threads, which contains:
 *    i. required operands,
 *    ii. pointer to current function to execute*/
/*
 * worker(threadstruct) {
 *   Spinlock until operation
 *   switch on current op,
 *   get relevant operands,
 *   call relevant function, with fetched operands
 * }
 *
 * threadstruct :
 *  id,
 *  pointer to global struct,
 *
 * global struct :
 *  void pointer, structure depends on current operation 
 */



/*
 * Check misses, cachegrind
 */


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

//typedef uint32_t v4si __attribute__ ((vector_size(128)));

static uint32_t g_seed = 0;

static ssize_t g_width = 0;
static ssize_t g_height = 0;
static ssize_t g_elements = 0;

static ssize_t g_nthreads = 1;

//enum operation {NONE, SCADD, SCMUL, MADD, MMUL, STADD, STSUB, STMUL, GSUM, GTRACE, GMIN, GMAX, GFREQ};
//static enum operation curr_op = NONE;

struct freq_arg {
    const uint32_t* matrix;
    size_t id;
    uint32_t value;
    uint32_t result;
};
typedef struct freq_arg freq_arg;

struct sum_arg {
    const uint32_t* matrix;
    size_t id;
    uint32_t result;
};
typedef struct sum_arg sum_arg;

struct scalar_arg {
    const uint32_t* matrix;
    size_t id;
    uint32_t scalar;
    uint32_t* result;
};
typedef struct scalar_arg scalar_arg;

struct matrix_arg {
    uint32_t* result;
    size_t id;
    const uint32_t* a;
    const uint32_t* b;
};
typedef struct matrix_arg matrix_arg;

struct strass_arg {
    size_t id;
    const uint32_t* a;
    size_t aw;
    size_t ah;
    size_t as;
    const uint32_t* b;
    size_t bw;
    size_t bh;
    size_t bs;
    uint32_t* c;
    size_t cs;
};
typedef struct strass_arg strass_arg;


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

/**
 * Displays a matrix whose stride differs from its width.
 */ 
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
    //uint32_t* matrix = (uint32_t*)aligned_alloc(16, g_elements*sizeof(uint32_t));
    //memset(matrix, 0, g_elements*sizeof(uint32_t));
    //return matrix;
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
uint32_t* old_scalar_add(const uint32_t* matrix, uint32_t scalar) {

    uint32_t* result = new_matrix();
    
    for (ssize_t i = 0; i < g_elements; ++i) {
        result[i] = matrix[i] + scalar;
    }

    return result;
}

static void *scal_add_worker(void *arg) {
    scalar_arg argument = *((scalar_arg *)arg);

    for (size_t y = argument.id; y < g_height; y += g_nthreads) {
        for (size_t x = 0; x < g_width; ++x) {
            argument.result[y*g_height + x] = argument.matrix[y*g_height + x] + argument.scalar;
        }
    }
    
    return NULL;
}

uint32_t* scalar_add(const uint32_t* matrix, uint32_t scalar) {
    uint32_t* result = new_matrix();

    size_t temp = g_nthreads;
    g_nthreads = g_nthreads < g_height ? g_nthreads : g_height;
    scalar_arg args[g_nthreads];

    for (size_t i = 0; i < g_nthreads; ++i) {
        args[i] = (scalar_arg) {
            .matrix = matrix,
            .id = i,
            .scalar = scalar,
            .result = result
        };
    }
    
    pthread_t thread_ids[g_nthreads];

    for (size_t i = 0; i < g_nthreads; ++i) {
        pthread_create(thread_ids + i, NULL, scal_add_worker, args + i);
    }
    
    for (size_t i = 0; i < g_nthreads; ++i) {
        pthread_join(thread_ids[i], NULL);
    }

    g_nthreads = temp;

    return result;
}

/**
 * Returns new matrix with scalar multiplied to each element
 */
uint32_t* old_scalar_mul(const uint32_t* matrix, uint32_t scalar) {

    uint32_t* result = new_matrix();

    for (ssize_t i = 0; i < g_elements; ++i) {
        result[i] = matrix[i] * scalar;
    }
    
    return result;
}

static void *scal_mul_worker(void *arg) {
    scalar_arg argument = *((scalar_arg *)arg);

    for (size_t y = argument.id; y < g_height; y += g_nthreads) {
        for (size_t x = 0; x < g_width; ++x) {
            argument.result[y*g_height + x] = argument.matrix[y*g_height + x] * argument.scalar;
        }
    }
    
    return NULL;
}

uint32_t* scalar_mul(const uint32_t* matrix, uint32_t scalar) {
    uint32_t* result = new_matrix();

    size_t temp = g_nthreads;
    g_nthreads = g_nthreads < g_height ? g_nthreads : g_height;
    scalar_arg args[g_nthreads];

    for (size_t i = 0; i < g_nthreads; ++i) {
        args[i] = (scalar_arg) {
            .matrix = matrix,
            .id = i,
            .scalar = scalar,
            .result = result
        };
    }
    
    pthread_t thread_ids[g_nthreads];

    for (size_t i = 0; i < g_nthreads; ++i) {
        pthread_create(thread_ids + i, NULL, scal_mul_worker, args + i);
    }
    
    for (size_t i = 0; i < g_nthreads; ++i) {
        pthread_join(thread_ids[i], NULL);
    }

    g_nthreads = temp;

    return result;
}

/**
 * Returns new matrix with elements added at the same index
 */
uint32_t* old_matrix_add(const uint32_t* matrix_a, const uint32_t* matrix_b) {

    uint32_t* result = new_matrix();

    for (ssize_t i = 0; i < g_elements; ++i) {
        result[i] = matrix_a[i] + matrix_b[i];
    }

    return result;
}

static void *add_worker(void *arg) {
    matrix_arg argument = *((matrix_arg *)arg);

    for (size_t y = argument.id; y < g_height; y += g_nthreads) {
        for (size_t x = 0; x < g_width; ++x) {
            argument.result[y*g_height + x] = argument.a[y*g_height + x] + argument.b[y*g_height + x];
        }
    }
    
    return NULL;
}

uint32_t* matrix_add(const uint32_t* matrix_a, const uint32_t* matrix_b) {
    uint32_t* result = new_matrix();

    size_t temp = g_nthreads;
    g_nthreads = g_nthreads < g_height ? g_nthreads : g_height;
    matrix_arg args[g_nthreads];

    for (size_t i = 0; i < g_nthreads; ++i) {
        args[i] = (matrix_arg) {
            .result = result,
            .id = i,
            .a = matrix_a,
            .b = matrix_b
        };
    }
    
    pthread_t thread_ids[g_nthreads];

    for (size_t i = 0; i < g_nthreads; ++i) {
        pthread_create(thread_ids + i, NULL, add_worker, args + i);
    }
    
    for (size_t i = 0; i < g_nthreads; ++i) {
        pthread_join(thread_ids[i], NULL);
    }

    g_nthreads = temp;
    
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
    
    strassen(matrix_a, g_width, g_height, g_width, matrix_b, g_width, g_height, g_width, result, g_width);

    return result;
}

/**
 * Adds a and b, placing the result in c.
 * Details are as strassen().
 */
void old_strass_add(const uint32_t* a, ssize_t aw, ssize_t ah, ssize_t as,
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
}

// TODO: see if removing these dereferences is faster.
static void *strass_add_worker(void *arg) {
    strass_arg argument = *((strass_arg *)arg);

    for (ssize_t y = argument.id; y < argument.ah; y += g_nthreads) {
        for (ssize_t x = 0; x < argument.aw; ++x) {
            argument.c[argument.cs*y + x] = argument.a[argument.as*y + x];
        }
    }

    for (ssize_t y = argument.id; y < argument.bh; y += g_nthreads) {
        for (ssize_t x = 0; x < argument.bw; ++x) {
            argument.c[argument.cs*y + x] += argument.b[argument.bs*y + x];
        }
    }
    
    return NULL;
}

void strass_add(const uint32_t* a, ssize_t aw, ssize_t ah, ssize_t as,
                 const uint32_t* b, ssize_t bw, ssize_t bh, ssize_t bs,
                 uint32_t* c, ssize_t cs) {
    size_t temp = g_nthreads;
    g_nthreads = g_nthreads < g_height ? g_nthreads : g_height;
    strass_arg args[g_nthreads];

    for (size_t i = 0; i < g_nthreads; ++i) {
        args[i] = (strass_arg) {
            .id = i,
            .a = a,
            .aw = aw,
            .ah = ah,
            .as = as,
            .b = b,
            .bw = bw,
            .bh = bh,
            .bs = bs,
            .c = c,
            .cs = cs
        };
    }
    
    pthread_t thread_ids[g_nthreads];

    for (size_t i = 0; i < g_nthreads; ++i) {
        pthread_create(thread_ids + i, NULL, strass_add_worker, args + i);
    }
    
    for (size_t i = 0; i < g_nthreads; ++i) {
        pthread_join(thread_ids[i], NULL);
    }

    g_nthreads = temp;
}

/**
 * Subtracts b from a, placing the result in c.
 * Details are as strassen().
 */ 
void old_strass_sub(const uint32_t* a, ssize_t aw, ssize_t ah, ssize_t as,
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
}

// TODO: see if removing these dereferences is faster.
static void *strass_sub_worker(void *arg) {
    strass_arg argument = *((strass_arg *)arg);

    for (ssize_t y = argument.id; y < argument.ah; y += g_nthreads) {
        for (ssize_t x = 0; x < argument.aw; ++x) {
            argument.c[argument.cs*y + x] = argument.a[argument.as*y + x];
        }
    }

    for (ssize_t y = argument.id; y < argument.bh; y += g_nthreads) {
        for (ssize_t x = 0; x < argument.bw; ++x) {
            argument.c[argument.cs*y + x] -= argument.b[argument.bs*y + x];
        }
    }
    
    return NULL;
}

void strass_sub(const uint32_t* a, ssize_t aw, ssize_t ah, ssize_t as,
                 const uint32_t* b, ssize_t bw, ssize_t bh, ssize_t bs,
                 uint32_t* c, ssize_t cs) {
    size_t temp = g_nthreads;
    g_nthreads = g_nthreads < g_height ? g_nthreads : g_height;
    strass_arg args[g_nthreads];

    for (size_t i = 0; i < g_nthreads; ++i) {
        args[i] = (strass_arg) {
            .id = i,
            .a = a,
            .aw = aw,
            .ah = ah,
            .as = as,
            .b = b,
            .bw = bw,
            .bh = bh,
            .bs = bs,
            .c = c,
            .cs = cs
        };
    }
    
    pthread_t thread_ids[g_nthreads];

    for (size_t i = 0; i < g_nthreads; ++i) {
        pthread_create(thread_ids + i, NULL, strass_sub_worker, args + i);
    }
    
    for (size_t i = 0; i < g_nthreads; ++i) {
        pthread_join(thread_ids[i], NULL);
    }

    g_nthreads = temp;
}

static void *strass_mul_worker(void *arg) {
    strass_arg argument = *((strass_arg *)arg);
    
    ssize_t min_p = (argument.aw < argument.bh ? argument.aw : argument.bh);

    for (ssize_t y = argument.id; y < argument.ah; y += g_nthreads) {
        for (ssize_t k = 0; k < min_p; k++) {
            for (ssize_t x = 0; x < argument.bw; x++) {
                argument.c[y*argument.cs + x] += argument.a[y*argument.as + k] * argument.b[k*argument.bs + x];
            }
        }
    }
   
    return NULL;
}

void strass_mul(const uint32_t* a, ssize_t aw, ssize_t ah, ssize_t as,
                 const uint32_t* b, ssize_t bw, ssize_t bh, ssize_t bs,
                 uint32_t* c, ssize_t cs) {

    size_t temp = g_nthreads;
    g_nthreads = g_nthreads < g_height ? g_nthreads : g_height;
    strass_arg args[g_nthreads];

    for (size_t i = 0; i < g_nthreads; ++i) {
        args[i] = (strass_arg) {
            .id = i,
            .a = a,
            .aw = aw,
            .ah = ah,
            .as = as,
            .b = b,
            .bw = bw,
            .bh = bh,
            .bs = bs,
            .c = c,
            .cs = cs
        };
    }
    
    pthread_t thread_ids[g_nthreads];

    for (size_t i = 0; i < g_nthreads; ++i) {
        pthread_create(thread_ids + i, NULL, strass_mul_worker, args + i);
    }
    
    for (size_t i = 0; i < g_nthreads; ++i) {
        pthread_join(thread_ids[i], NULL);
    }

    g_nthreads = temp;
}

/**
 * Base case multiplication, 
 * reverts to this for small-enough matrices.
 * Details are as strassen().
 */
void old_strass_mul(const uint32_t* a, ssize_t aw, ssize_t ah, ssize_t as,
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


/**
 * Multiplies a by b, placing the result in c.
 *
 * In the case where one matrix is smaller than the other,
 * that matrix is implicitly padded with zeroes on that edge.
 *
 * c is assumed to be exactly large enough to hold the resulting product.
 *
 * a, b, c may all be submatrices, being part of a larger parent in memory.
 * Hence, while aw, ah, bw, bh define the dimensions of a and b;
 * as, bs, cs are the stride required to get to the corresponding next row.
 * That is, it is the width of the parent matrix.
 */
void strassen(const uint32_t* a, ssize_t aw, ssize_t ah, ssize_t as,
              const uint32_t* b, ssize_t bw, ssize_t bh, ssize_t bs,
              uint32_t* c, ssize_t cs) {
    
    // All matrices passed in are square; no need to check both dimensions here.
    if (aw <= STRASS_THRESH || bw <= STRASS_THRESH) {
        strass_mul(a, aw, ah, as, b, bw, bh, bs, c, cs);
        return;
    } else {
         
        // Sizes of, and pointers to, each quadrant of input/output matrices.
        
        ssize_t a11w, a11h, a22w, a22h,
                b11w, b11h, b22w, b22h,
                c11w, c11h, c22w, c22h;
        
        if (aw > bw) {
            a11w = aw&1 ? (aw+1)/2 : aw/2;
            b11w = a11w;
        } else {
            b11w = bw&1 ? (bw+1)/2 : bw/2;
            a11w = b11w;
        }

        if (ah > bh) {
            a11h = ah&1 ? (ah+1)/2 : ah/2;
            b11h = a11h;
        } else {
            b11h = bh&1 ? (bh+1)/2 : bh/2;
            a11h = b11h;
        }

        a22w = aw-a11w;
        a22h = ah-a11h;
        b22w = bw-b11w;
        b22h = bh-b11h;
        c11w = a11w;
        c11h = a11h; 
        c22w = (aw > bw ? aw : bw)-c11w;
        c22h = (ah > bh ? ah : bh)-c11h;

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
        

        // Allocate room for intermediate results.
        
        ssize_t m = c11w;
        ssize_t n = m*m;
        uint32_t* S = (uint32_t*)calloc(11, sizeof(uint32_t) * n);
      

        // Perform actual strassen algorithm.
        
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

        strassen(a11, a11w, a11h, as, S + n, m, m, m, S, m); //P1 = A11*S1
        memset(S + n, 0, n*sizeof(uint32_t));
        strassen(S + 2*n, m, m, m, b22, b22w, b22h, bs, S + n, m); //P2 = S2*B22
        memset(S + 2*n, 0, n*sizeof(uint32_t));
        strassen(S + 3*n, m, m, m, b11, b11w, b11h, bs, S + 2*n, m); //P3 = S3*B11
        memset(S + 3*n, 0, n*sizeof(uint32_t));
        strassen(a22, a22w, a22h, as, S + 4*n, m, m, m, S + 3*n, m); //P4 = A22*S4
        memset(S + 4*n, 0, n*sizeof(uint32_t));
        strassen(S + 5*n, m, m, m, S + 6*n, m, m, m, S + 4*n, m); //P5 = S5*S6
        memset(S + 5*n, 0, 2*n*sizeof(uint32_t));
        strassen(S + 7*n, m, m, m, S + 8*n, m, m, m, S + 5*n, m); //P6 = S7*S8
        strassen(S + 9*n, m, m, m, S + 10*n, m, m, m, S + 6*n, m); //P7 = S9*S10
        
        // No need to zero these, since entries of S are the same size, c is already zeroed.
        strass_add(S + 4*n, m, m, m, S + 3*n, m, m, m, S + 10*n, m);
        strass_sub(S + 5*n, m, m, m, S + n, m, m, m, S + 9*n, m);
        strass_add(S + 10*n, c11w, c11h, m, S + 9*n, c11w, c11h, m, c11, cs); // C11 = P5 + P4 - P2 + P6
        strass_add(S, c22w, c11h, m, S + n, c22w, c11h, m, c12, cs); // C12 = P1 + P2
        strass_add(S + 2*n, c11w, c22h, m, S + 3*n, c11w, c22h, m, c21, cs); // C21 = P3 + P4
        strass_add(S, m, m, m, S + 4*n, m, m, m, S + n, m);
        strass_add(S + 2*n, m, m, m, S + 6*n, m, m, m, S + 3*n, m);
        strass_sub(S + n, c22w, c22h, m, S + 3*n, c22w, c22h, m, c22, cs); // C22 = P1 + P5 - P3 - P7
        
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
            strassen((mpowers + g_elements*(i-1)), g_width, g_height, g_width,
                     (mpowers + g_elements*(i-1)), g_width, g_height, g_width,
                     (mpowers + g_elements*i), g_width);
        }
    }
    
    uint32_t* temp = new_matrix();

    for (uint32_t i = 0; i < width; ++i) {
        if (exponent & (1 << i)) {
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
    
    return result;
}

////////////////////////////////
///       COMPUTATIONS       ///
////////////////////////////////

static void *sum_worker(void *arg) {
    sum_arg argument = *((sum_arg *)arg);
    
    uint32_t sum = 0;

    for (size_t y = argument.id; y < g_height; y += g_nthreads) {
        for (size_t x = 0; x < g_width; ++x) {
            sum += argument.matrix[y*g_height + x];
        }
    }
    
    ((sum_arg*)arg)->result += sum;

    return NULL;
}

uint32_t get_sum(const uint32_t* matrix) {

    uint32_t sum = 0;
    size_t temp = g_nthreads;
    g_nthreads = g_nthreads < g_height ? g_nthreads : g_height;
    sum_arg args[g_nthreads];

    for (size_t i = 0; i < g_nthreads; ++i) {
        args[i] = (sum_arg) {
            .matrix = matrix,
            .id = i,
            .result = 0
        };
    }
    
    pthread_t thread_ids[g_nthreads];

    for (size_t i = 0; i < g_nthreads; ++i) {
        pthread_create(thread_ids + i, NULL, sum_worker, args + i);
    }
    
    for (size_t i = 0; i < g_nthreads; ++i) {
        pthread_join(thread_ids[i], NULL);
    }

    for (size_t i = 0; i < g_nthreads; ++i) {
        sum += args[i].result;
    }
    
    g_nthreads = temp;

    return sum;
}

/**
 * Returns the sum of all elements
 */
uint32_t old_get_sum(const uint32_t* matrix) {

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

static void *min_worker(void *arg) {
    sum_arg argument = *((sum_arg *)arg);
    
    uint32_t min = UINT32_MAX;

    for (size_t y = argument.id; y < g_height; y += g_nthreads) {
        for (size_t x = 0; x < g_width; ++x) {
            if (argument.matrix[y*g_height + x] < min) min = argument.matrix[y*g_height + x];
        }
    }
    
    ((sum_arg*)arg)->result = min;

    return NULL;
}

uint32_t get_minimum(const uint32_t* matrix) {

    uint32_t min = matrix[0];
    size_t temp = g_nthreads;
    g_nthreads = g_nthreads < g_height ? g_nthreads : g_height;
    sum_arg args[g_nthreads];

    for (size_t i = 0; i < g_nthreads; ++i) {
        args[i] = (sum_arg) {
            .matrix = matrix,
            .id = i,
            .result = 0
        };
    }
    
    pthread_t thread_ids[g_nthreads];

    for (size_t i = 0; i < g_nthreads; ++i) {
        pthread_create(thread_ids + i, NULL, min_worker, args + i);
    }
    
    for (size_t i = 0; i < g_nthreads; ++i) {
        pthread_join(thread_ids[i], NULL);
    }

    for (size_t i = 0; i < g_nthreads; ++i) {
        if (args[i].result < min) min = args[i].result;
    }
    
    g_nthreads = temp;

    return min;
}

/**
 * Returns the smallest value in the matrix
 */
uint32_t old_get_minimum(const uint32_t* matrix) {

    uint32_t min = matrix[0];

    for (uint32_t i = 1; i < g_elements; ++i) {
        if (matrix[i] < min) {
            min = matrix[i];
        }
    } 
    
    return min;
}

static void *max_worker(void *arg) {
    sum_arg argument = *((sum_arg *)arg);
    
    uint32_t max = 0;

    for (size_t y = argument.id; y < g_height; y += g_nthreads) {
        for (size_t x = 0; x < g_width; ++x) {
            if (argument.matrix[y*g_height + x] > max) max = argument.matrix[y*g_height + x];
        }
    }
    
    ((sum_arg*)arg)->result = max;

    return NULL;
}

uint32_t get_maximum(const uint32_t* matrix) {

    uint32_t max = matrix[0];
    size_t temp = g_nthreads;
    g_nthreads = g_nthreads < g_height ? g_nthreads : g_height;
    sum_arg args[g_nthreads];

    for (size_t i = 0; i < g_nthreads; ++i) {
        args[i] = (sum_arg) {
            .matrix = matrix,
            .id = i,
            .result = 0
        };
    }
    
    pthread_t thread_ids[g_nthreads];

    for (size_t i = 0; i < g_nthreads; ++i) {
        pthread_create(thread_ids + i, NULL, max_worker, args + i);
    }
    
    for (size_t i = 0; i < g_nthreads; ++i) {
        pthread_join(thread_ids[i], NULL);
    }

    for (size_t i = 0; i < g_nthreads; ++i) {
        if (args[i].result > max) max = args[i].result;
    }
    
    g_nthreads = temp;

    return max;
}

/**
 * Returns the largest value in the matrix
 */
uint32_t old_get_maximum(const uint32_t* matrix) {

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
uint32_t old_get_frequency(const uint32_t* matrix, uint32_t value) {
    
    uint32_t freq = 0;

    for (uint32_t i = 0; i < g_elements; ++i) {
        if (matrix[i] == value) ++freq;
    }

    return freq;
}

static void *freq_worker(void *arg) {
    freq_arg argument = *((freq_arg *)arg);
    
    uint32_t freq = 0;

    for (size_t y = argument.id; y < g_height; y += g_nthreads) {
        //printf("thread %zd computing row %zd\n", argument->id, y);
        for (size_t x = 0; x < g_width; ++x) {
            if (argument.matrix[y*g_height + x] == argument.value) ++freq;
        }
    }
    
    ((freq_arg*)arg)->result += freq;

    return NULL;
}

/**
 * Returns the frequency of the value in the matrix
 */
uint32_t get_frequency(const uint32_t* matrix, uint32_t value) {

    uint32_t freq = 0;
    size_t temp = g_nthreads;
    g_nthreads = g_nthreads < g_height ? g_nthreads : g_height;
    freq_arg args[g_nthreads];

    for (size_t i = 0; i < g_nthreads; ++i) {
        args[i] = (freq_arg) {
            .matrix = matrix,
            .id = i,
            .value = value,
            .result = 0
        };
    }
    
    pthread_t thread_ids[g_nthreads];

    for (size_t i = 0; i < g_nthreads; ++i) {
        pthread_create(thread_ids + i, NULL, freq_worker, args + i);
    }
    
    for (size_t i = 0; i < g_nthreads; ++i) {
        pthread_join(thread_ids[i], NULL);
    }

    for (size_t i = 0; i < g_nthreads; ++i) {
        freq += args[i].result;
    }
    
    g_nthreads = temp;

    return freq;
}
