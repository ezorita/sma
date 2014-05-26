#include "sma.h"

// TODO:
// - When increasing index, check whether the num of remaining nucleotides is less than the available tau,
//   in this case do not continue beacuse it's not possible to create seqs of distance tau.
// - Decrease the sequence len when inserting unknown nucleotides after a deletion!! Otherwise there'll be
//   repeated sequences! (a param to count the accumulated 'unkn' may be necessary).

void SIGSEGV_handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(EXIT_FAILURE);
}


int main(int argc, char *argv[])
{
   if (argc != 4) {
      fprintf(stderr, "usage: sma <seq> <seqlen> <tau>\n");
      exit(EXIT_FAILURE);
   }

   signal(SIGSEGV, SIGSEGV_handler);
   
   char  * sequence = argv[1];
   int slen = atoi(argv[2]);
   int tau  = atoi(argv[3]);

   int nids;
   seq_t * seqs = seqtoid(sequence, &nids, slen);
   seq_t seq = seqs[0];
   free(seqs);

   if (nids > 1)
      fprintf(stderr, "WARNING: ignoring 'N' (%s will be used)\n", idtoseq(seq, slen));

   // Prepare extras.   
   char * extra = malloc((2*tau - 1)*sizeof(char));
   extra += tau - 1;

   // Fill the seq part.
   int i = 1;
   for (; i < tau; i++)
      extra[i] = (seq >> (2*(i-1))) & 3;

   // Fill the extra part.
   int nextras = strlen(sequence + slen);
   if (nextras > 0)
      seqs = seqtoid(sequence + slen, &nids, nextras);

   seq_t remainder = seqs[0];
   free(seqs);

   i = 0;
   for (; i < min(nextras, tau); i++)
      extra[-i] = (remainder >> 2*(nextras - 1 - i)) & 3;
   for(; i < tau; i++)
      extra[-i] = -1;

   // Print seqs and extras.
   fprintf(stdout, "sequence: %#010X\nextras:\n", seq);
   for (i = -tau+1; i < tau; i++) fprintf(stdout,"[%d] %d\n",i,extra[i]);

   // Allocate stack.
   mstack_t ** mstack = malloc((tau+1)*sizeof(mstack_t *));
   for (int i = 0; i <= tau; i++) {
      int ssize = 1 << (4*i);
      mstack[i] = malloc(3*sizeof(seq_t) + ssize*sizeof(mismatch_t));
      mstack[i]->lim = ssize;
      mstack[i]->pos = 0;
      mstack[i]->seq = 0xFFFFFFFF;
   }

   sma(mstack, seq, slen, extra, tau);
   fprintf(stderr, "Number of seqs %d.\n", mstack[tau]->pos);

   return 0;
}

