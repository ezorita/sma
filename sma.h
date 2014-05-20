


#ifndef __SMA_
#define __SMA__

#define MAXVAL 0xFFFFFFFF
#define sub_nt(seq, pos, base) ((seq & ~(3 << pos)) + (base & (3 << pos))
#define ins_nt(seq, pos, base) ((seq & (MAXVAL << pos)) + ((seq & ~(MAXVAL << pos)) >> 2) + (base & (3 << pos)))

typedef unsigned int seq_t; // Longest seq: 16 nt (32 bits)
typedef struct ustack_t ustack_t;
typedef struct param_t  param_t;

struct seqstack_t {
   seq_t pos;
   seq_t lim;
   seq_t seq[];
};

struct param_t {
   seq_t         sign;
   seq_t         extra;
   int           extracnt;
   seqstack_t ** seqstack;
   // sign:
   //   bit 0: subseq sign (0 neg, 1 pos)
   //   bit 1: 2-subseq sign (0 neg, 1 pos)
   // extra: (used for deletions)
   //   Known nucleotides (count and seq) that follow the sequence.
};

void fill_params(param_t *params, seq_t seq, seq_t extra, int extracnt);

#endif
