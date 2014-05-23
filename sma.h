#include "plutocore.h"
#include <execinfo.h>
#include <signal.h>

#ifndef __SMA__
#define __SMA__


#define MAXVAL 0xFFFFFFFF

#define sub_nt(seq, pos, base)  ((seq & ~(3 << pos)) + ((base & 3) << pos))
#define ins_nt(seq, pos, base)  ((seq & (MAXVAL << (pos + 2))) + ((seq >> 2) & ~(MAXVAL << pos)) + ((base & ~(MAXVAL << 2)) << pos))
#define del_nt(seq, pos, extra) ((seq & (MAXVAL << (pos + 2))) + ((seq << 2) & ~(MAXVAL << (pos + 2))) + (extra & 3))
#define get_suffix(seq, len)    (seq & ~(MAXVAL << pos))
#define ins_check(seq, pos)     (get_suffix(seq, pos - 2) != get_suffix(seq >> 2, pos - 2))
#define del_check(seq, pos, extra) (get_suffix(seq, pos - 2) != get_suffix(del_nt(seq, pos, extra), pos - 2))

typedef unsigned int seq_t; // Longest seq: 16 nt (32 bits)
typedef struct mstack_t mstack_t;
typedef struct param_t  param_t;
typedef struct mismatch_t mismatch_t;

struct mismatch_t {
   char  offset;
   seq_t seq;
};

struct mstack_t {
   seq_t        seq;
   seq_t        pos;
   seq_t        lim;
   mismatch_t   m[];
};

struct param_t {
   char        slen;
   char      * extra;
   mstack_t ** mstack;
};

void sma (mstack_t **, seq_t, seq_t, char *, int);
void sma_gen (seq_t, int, char, char, param_t, char);
void sma_merge (mstack_t **, int);
void add_mismatch (mstack_t **, seq_t, char);
int  mcomp (const void *, const void *);
void print_stack (mstack_t *, int);

#endif
