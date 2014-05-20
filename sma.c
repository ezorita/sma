#include "sma.h"

// TODO:
// - When increasing index, check whether the num of remaining nucleotides is less than the available tau,
//   in this case do not continue beacuse it's not possible to create seqs of distance tau.
// - Decrease the sequence len when inserting unknown nucleotides after a deletion!! Otherwise there'll be
//   repeated sequences! (a param to count the accumulated 'unkn' may be necessary).

ustack_t * sma (uint seq, uint slen, unsigned char tau) {
   
}


void fill_params (param_t * params, seq_t seq, seq_t extra, seq_t extracnt) {
   
}

void rec_sma (seq_t seq, int slen, int start, int tau, int offset, param_t params, seq_t *count) {

   // Write sequence if tau is exhausted.
   if (tau == 0) {
      stack_add(seq, params.seqstack + offset);
      return;
   } 

   // Return when it's not possible to create seqs of distance tau anymore.
   if (slen - start < tau) return;
   
   // Create a mask for the current nucleotide.
   unsigned char pos = 2*(slen - start - 1);
   seq_t ntmask = 3 << pos;
   seq_t inc    = 1 << pos;

   // Iterate:
   do {
      // Go one level deeper instead of substituting for the same nucleotide.
      if ((seq & ntmask) == ((*count) & ntmask)) {
         rec_sma(seq, slen, start+1, tau, offset, params, count);
      }
      // The new nucleotide is not the same as the original.
      else {
         // If this is the last nt, only substitute.
         if (start == slen - 1) {
            stack_add(sub_nt(seq, pos, (*count)), params.seqstack + offset);
         }
         // The new nucleotide is different than the next one as well: Just insert/substitute.
         else if (((seq << 2) & ntmask) != ((*count) & ntmask)) {
            // If subsequence is positive: INS -> SUB.
            if (params.sign & inc) {
               // INS. (offset + 1 to count the insertion)
               seq_t subcount = 0;
               rec_sma(ins_nt(seq, pos, (*count)), slen, start + 1, tau - 1, offset + 1, params, &subcount);
               // SUB.
               subcount = 0;
               rec_sma(sub_nt(seq, pos, (*count)), slen, start + 1, tau - 1, offset, params, &subcount);
            }
            // If subsequence is negative: SUB -> INS.
            else {
               // SUB.
               seq_t subcount = 0;
               rec_sma(sub_nt(seq, pos, (*count)), slen, start + 1, tau - 1, offset, params, &subcount);
               // INS. (offset + 1 to count the insertion)
               subcount = 0;
               rec_sma(ins_nt(seq, pos, (*count)), slen, start + 1, tau - 1, offset + 1, params, &subcount);
            }
         }
         // The new nucleotide is the same as the next one: Time for deletions!
         else {
            
         }
      }
      (*count) += inc;
    
   } while (((*count) & ntmask) < ntmask);
   
}
