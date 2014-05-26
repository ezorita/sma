#include "sma.h"

void
print_stack
(
 mstack_t * stack,
 int slen
 )
{
   for (int k = 0; k < stack->pos; k++) {
      char * seq = idtoseq(stack->m[k].seq, slen);
      fprintf(stdout,"%s\n", seq);
      free(seq);
   }
}

void
sma
(
 mstack_t ** mstack,
 seq_t       seq,
 seq_t       slen,
 char      * extra,
 int         tau
)

{
   if (mstack[tau]->seq == seq) return;

   // Clear stack.
   mstack[tau]->pos = 0;
   if (tau == 0) {
      mstack[tau]->seq = seq;
      add_mismatch(mstack, seq, 0);
      return;
   }

   // Check whether the previous has been computed.
   if (mstack[tau-1]->seq != seq) sma(mstack, seq, slen, extra, tau-1);

   // Run recursive generator.
   struct param_t params = {
      .slen = slen,
      .extra = extra,
      .mstack = mstack + tau
   };
   mstack[tau]->seq = seq;
   sma_gen(seq, 0, tau, 0, params, 0);

   // Sort.
   mstack[tau] = realloc(mstack[tau], 3*sizeof(seq_t) + mstack[tau]->pos*sizeof(mismatch_t));
   mstack[tau]->lim = mstack[1]->pos;
   qsort(mstack[tau]->m, mstack[tau]->pos, sizeof(mismatch_t), mcomp);

   sma_merge(mstack, tau);
}


void
sma_merge
(
 mstack_t ** mstackp,
 int         tau
)
// SYNOPSIS:                                                              
//   Merges a sorted stack containing mismatches of distance leq than tau with all the
//   previously-merged stacks of distance tau-1, tau-2, ... 0 from the original sequence.
//   As a result of the merge, the stack contains only sequences that are exactly at
//   distance tau from the original seq. For the algorithm to work, all the merged
//   stacks of distance leq than the specified tau must be present at the stack array.
//   Merging a previously merged stack returns exactly the same result.
//
// PARAMETERS:                                                            
//   mstackp: a pointer to the beginning of the stack array.
//   tau:     the stack that will be merged.
//                                                                        
// RETURN:                                                                
//   The function does not return anything. The merged list of mismatches is stored in the
//   same original stack. The final stack is realloc'ed to fit exactly the final size of the
//   set.
//                                                                        
// SIDE EFFECTS:                                                          
//   The tau-th stack of the array may be realloc'ed.
//   
{
   // Get the pointer.
   mstack_t   * dst_stack = mstackp[tau];
   seq_t        dst_sz    = dst_stack->pos;
   mismatch_t * dst_lst   = dst_stack->m;

   // Remove duplicates first.
   int p = 0;
   for (int i = 0; i < dst_sz - 1; i++)
      if (dst_lst[i].seq != dst_lst[i+1].seq) dst_lst[p++] = dst_lst[i];

   dst_lst[p++] = dst_lst[dst_sz - 1];
   // Save new size.   
   dst_stack->pos = p;

   // Merge the current stack with all the previous.
   for (int i = 0; i < tau; i++) {
      mstack_t   * ref_stack = mstackp[i];
      seq_t        ref_sz    = ref_stack->pos;
      mismatch_t * ref_lst   = ref_stack->m;

      int j = 0, k = 0, l = 0;
      while (j < ref_sz && k < dst_sz) {
         if (ref_lst[j].seq < dst_lst[k].seq)
            j++;
         else if (ref_lst[j].seq > dst_lst[k].seq) 
            dst_lst[l++] = dst_lst[k++];
         else {
            j++; k++;
         }
      }

      // Complete copy.
      if (k < dst_sz && k != l) 
         memmove(dst_lst + l, dst_lst+ k, (dst_sz - k)*sizeof(mismatch_t));

      // Save the final size.
      dst_sz -= k-l;
   }
   
   dst_stack->pos = dst_sz;

   // Job done, realloc destination stack.
   dst_stack->lim = dst_stack->pos;
   mstackp[tau] = realloc(dst_stack, 3*sizeof(seq_t) + dst_stack->pos*sizeof(mismatch_t));
   if (mstackp[tau] == NULL) {
      fprintf(stderr, "Error (realloc) in sma_merge.\n");
      exit(EXIT_FAILURE);
   }
}


void
sma_gen
(
 seq_t   seq,
 int     start,
 char    tau,
 char    offset,
 param_t params,
 char    last
)
// SYNOPSIS:                                                              
//   Recursively generates all INS/SUB/DEL possible combinations up to tau mismatches.
//   The resulting set may contain sequences within distance < tau wrt the original seq,
//   due to redundant mismatches.
//   Only straightforward redundant mismatches are avoided:
//     - Inserting/deleting within repeated nucleotides.
//     - Del+Ins and Ins+Del.
//     - Inserting deleting with a constant suffix.
//     - delete then substitute by the previous (same as deleting the next).
//
// PARAMETERS:                                                            
//   seq:    the sequence id.
//   start:  start level (nucleotide).
//   tau:    remaining num of mismatches.
//   offset: offset of the current sequence. (effective num of dels or ins).
//   params: parameter struct.
//   last:   informs about the last action, set to 0. (used in recursive generation).
//                                                                        
// RETURN:                                                                
//   The set of sequences is returned in the mismatch stack provided in params.
//                                                                        
// SIDE EFFECTS:                                                          
//   The mismatch stack may be realloc'ed.
//   
{
   // Tau exhausted, save seq and return.
   if (tau == 0) {
      // Remember to clear the upper bits!
      add_mismatch(params.mstack, seq & ~(MAXVAL << (2*params.slen)), offset);
      return;
   }

   // End of chain
   if (start == params.slen) return;

   // Create a mask for the current nucleotide.
   unsigned char pos = 2*(params.slen - start - 1);

   // Count.
   char base = (seq >> pos) & 3;
   char next = (seq >> (pos - 2)) & 3;

   // SUBSTITUTIONS:
   char deleted = (last < 0 ? -last-1 : -1);
   for (int i = 0; i < 4; i++) {
      if (i == base || i == deleted) continue;
      sma_gen(sub_nt(seq, (seq_t)pos, (seq_t)i), start+1, tau-1, offset, params, 0);
   }

   if (start < params.slen - 1) {
      // DELETIONS: (Only if q[i] != q[i+1] and previous was not an insertion)
      if ((base != next) && del_check(seq, pos, params.extra[offset]) && last <= 0) {
         if (params.extra[offset] < 0) {
            for (int i = 0; i < 4; i++) {
               sma_gen(del_nt(seq, (seq_t)pos, (seq_t)i), start, tau-1, offset-1, params, -base-1);
            }
         }
         else {
            sma_gen(del_nt(seq, (seq_t)pos, (seq_t)params.extra[offset]), start, tau-1, offset-1, params, -base-1);
         }
      }
   
      // INSERTIONS:
      if (ins_check(seq, pos) && (last >= 0)) {
         for (int i = 0; i < 4; i++) {
            if (i == base) continue;
            sma_gen(ins_nt(seq,(seq_t)pos, (seq_t)i), start+1, tau-1, offset+1, params, base + 1);
         }
      }
   }

   // Continue w/ next nucleotide.
   sma_gen(seq, start+1, tau, offset, params, 0);
}

void
add_mismatch
(
 mstack_t ** stackp,
 seq_t       seq,
 char        offset
)
{
   mstack_t * mstack = *stackp;
   
   if (mstack->pos >= mstack->lim) {
      size_t newsize = mstack->lim * 2;
      mstack_t * p = realloc(mstack, 3*sizeof(seq_t) + newsize * sizeof(mismatch_t));
      if (p == NULL) {
         fprintf(stderr, "mismatch realloc failed: %s", strerror(errno));
         exit(EXIT_FAILURE);
      }
      *stackp = mstack = p;
      mstack->lim = newsize;
   }

   struct mismatch_t mismatch = {
      .offset = offset,
      .seq    = seq
   };

   mstack->m[mstack->pos++] = mismatch;
}


int
mcomp
(
 const void * a,
 const void * b
 )
{
   mismatch_t * ma = (mismatch_t *) a;
   mismatch_t * mb = (mismatch_t *) b;
   if (ma->seq > mb->seq)
      return 1;
   else if (ma->seq < mb->seq)
      return -1;
   else
      return 0;
}
