/*-------------------------------------------------------------

Copyright (C) 2000 Peter Clote. 
All Rights Reserved.

Permission to use, copy, modify, and distribute this
software and its documentation for NON-COMMERCIAL purposes
and without fee is hereby granted provided that this
copyright notice appears in all copies.


THE AUTHOR AND PUBLISHER MAKE NO REPRESENTATIONS OR
WARRANTIES ABOUT THE SUITABILITY OF THE SOFTWARE, EITHER
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
PARTICULAR PURPOSE, OR NON-INFRINGEMENT. THE AUTHORS
AND PUBLISHER SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED
BY LICENSEE AS A RESULT OF USING, MODIFYING OR DISTRIBUTING
THIS SOFTWARE OR ITS DERIVATIVES.

-------------------------------------------------------------*/


/*************************************************
	Original Program: smithWaterman.c
	Peter Clote, 11 Oct 2000

Program for local sequence alignment, using the Smith-Waterman
algorithm and assuming a LINEAR gap penalty.
A traceback is used to determine the alignment, and
is determined by choosing that direction, for
which S(i-1,j-1)+sigma(a_i,b_j), S(i-1,j)+Delta and 
S(i,j-1)+Delta is maximum, i.e.  for which 

                    _
                   |
                   | H(i-1,j-1) + sigma(a_i,b_j)  (diagonal)
H(i,j) =  MAX of   | H(i-1,j)   + delta           (up)
                   | H(i,j-1)   + delta           (left)
                   | 0
                   | _


is a maximum.

*************************************************/

/* Maximum size of the sequences: they should have equal size */

/*
 FIX: SOME BUGS HAS BEEN SOLVED
      SOME CODE HAS BEEN IMPROVED 
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>		// assertion handling
#include <ctype.h>		// character handling
#include <stdlib.h>		// def of RAND_MAX

   /* Note:                      	 */
   /* N must be the size of the arrays   */
   /* Here it is assumed that the two	 */
   /* sequences have the same size       */

/* PARALLEL CODE */
#include <mpi.h>		// Use MPI

/* The following C construct must be understood as:
 *    condition ? code_when_true : code when false 
 * is a shorthand for
 *     if (condition) { code_when_true  } 
 *        else        { code_when_false } */
#define min(a,b) (((a)<(b)) ? (a) : (b)) /* Complete the value returned when
					    the condition evaluates to false. */

#define MAX_SEQ 50

/* Macro to ease checking for errors when allocating memory dynamically. */
#define CHECK_NULL(_check) {\
   if ((_check)==NULL) \
      fprintf(stderr, "Null Pointer allocating memory\n");\
   }

#include <sys/time.h>
double getusec_() {
        struct timeval time;
        gettimeofday(&time, NULL);
        return ((double)time.tv_sec * (double)1e6 + (double)time.tv_usec);
}

#define START_COUNT_TIME stamp = getusec_();
#define STOP_COUNT_TIME(_m) stamp = getusec_() - stamp;\
                        stamp = stamp/1e6;\
                        printf ("%s: %0.6f\n",(_m), stamp);


#define AA 20			// Maximum dimension for the similarity matrix
#define MAX2(x,y)     ((x)<(y) ? (y) : (x))
#define MAX3(x,y,z)   (MAX2(x,y)<(z) ? (z) : MAX2(x,y))

// function prototypes
void error (char *);		/** error handling */

int char2AAmem[256];
int AA2charmem[AA];
void initChar2AATranslation (void);


/* PARALLEL CODE */
int getRowCount(int rowsTotal, int mpiRank, int mpiSize) {
  /* Adjust slack of rows in case rowsTotal is not exactly divisible */
  return (rowsTotal / mpiSize) + (rowsTotal % mpiSize > mpiRank);
}


int
main (int argc, char *argv[])
{

  // variable declarations
  FILE *in1, *in2, *sm;
  char ch;
  int temp;
  int i, j, diag, down, right, DELTA;
  //int tempi, tempj, x, y;
  //int topskip, bottomskip;
  //char *s1out, *s2out;  /* Pointers to arrays that will store the two results */
  //int Aend, Bend, Abegin, Bbegin;
  int max, Max, xMax, yMax;
  // Max is first found maximum in similarity matrix H
  // max is auxilliary to determine largest of
  // diag,down,right, xMax,yMax are h-coord of Max
  short *s1, *s2;  /* Pointers to arrays that will store the two sequences */
  int *hptr;			// Scoring matrix
  int **h;			// accessed as h[i][j]
  int sim[AA][AA];		// Similarity matrix
  //short **xTraceback, **yTraceback;
  //short *xTracebackptr, *yTracebackptr;
  int N; /* Lenght of sequence to analyze */
  int dim1; /* Lenght of 1st sequence to analyze */
  int dim2; /* Lenght of 2nd sequence to analyze */
  int nc; /* Number of characters read with fscanf */
  int BS; /* Block size */

  /* PARALLEL CODE */
  /* Local matrices */
  short *s1l;  /* Pointers to local arrays that will store the 1st sequence */
  // Following variables are not needed unless traceback is implemented:
  //int *hlptr;			// Scoring matrix
  //int **hl;			// accessed as hl[i][j]

  int rank, nprocs, nrows;
  MPI_Status status;

  double stamp;
  START_COUNT_TIME;

  /* PARALLEL CODE */
  MPI_Init (&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &nprocs );


  if (!rank ) { /* Use process 0 as the master process and have it process the input */
   /**** Error handling for input file ****/
   if (argc < 6) 
    {
      fprintf (stderr, "Usage:\n");
      fprintf (stderr, " mpirun.mpich -np P %s Seq1 Seq2 SM gapPenalty BS [N]\n", argv[0]);
      fprintf (stderr, "   where:\n");
      fprintf (stderr, "     P: Number of MPI processes\n");
      fprintf (stderr, "     Seq1: File containing 1st sequence\n");
      fprintf (stderr, "     Seq2: File containing 2nd sequence\n");
      fprintf (stderr, "     SM: File containing Similarity Matrix\n");
      fprintf (stderr, "     gapPenalty: Penalty of a gap\n");
      fprintf (stderr, "     BS: Block size for the parallel implementation\n");
      fprintf (stderr, "     [N]: Max length of sequences to use (optional)\n");
      MPI_Abort(MPI_COMM_WORLD, 911);
      exit (1);
    }


   /***** Initialization of input file and pair array **/
   in1 = fopen (argv[1], "r");
   in2 = fopen (argv[2], "r");
   sm = fopen (argv[3], "r");
   DELTA = atoi (argv[4]);
   BS = atoi (argv[5]);
   N = MAX_SEQ;
   if (argc >= 7)
     N = atoi (argv[6]);
   if (argc > 7)
     printf ("WARNING: Additional parameters after argv[6]=%s not used.\n",
	     argv[6]);

   assert(in1);
   assert(in2);
   assert(sm);
   assert(N > 0);
   assert((0 < BS) && (BS <= N));

   printf("Gap Penalty=%d\n", DELTA);
   printf("Block Size used BS=%d\n", BS);

   CHECK_NULL ((s1 = (short *) malloc (sizeof (short) * (N + 1))));
   CHECK_NULL ((s2 = (short *) malloc (sizeof (short) * (N + 1))));

   initChar2AATranslation ();

   /** Read similarity matrix **/
   if (fscanf (sm, "%*s") == EOF) /* Initial string in file not used.
				     Instead, the coding is hard coded via
				     routine initChar2AATranslation (void) */
     {
       fprintf (stderr, "Similarity Matrix file empty\n");
       fclose (sm);
       exit (1);
     } 
   for (i = 0; i < AA; i++)  /* Use initial AA rows in the matrix in the file */
     for (j = 0; j <= i; j++)
       {
	 if (fscanf (sm, "%d ", &temp) == EOF)
	   {
	     fprintf (stderr, "Similarity Matrix file empty\n");
	     fclose (sm);
	     exit (1);
	   }
	 sim[i][j] = temp;
       }
   fclose (sm);
   for (i = 0; i < AA; i++)
     for (j = i + 1; j < AA; j++)
       sim[i][j] = sim[j][i];	// symmetrify



   /** Read first file into array "s1" **/
 
   i = 0;
   do
     {
       nc = fscanf (in1, "%c", &ch);
       if (nc > 0 && char2AAmem[(int) ch] >= 0)
	 {
	   s1[++i] = char2AAmem[(int) ch];
	 }
     }
   while (nc > 0 && (i < N));

   dim1 = i;	/* We store the length of the 1st sequence analyzed in dim1 */
   fclose (in1);
   //assert(dim1 >= N);
   if (dim1<N)
     printf("WARNING: Sequence 1 read from file %s has %d values (<N=%d)\n",
	    argv[1], dim1, N);
   printf("Length of sequence 1 analyzed=%d\n", dim1);
  

   /** Read second file into array "s2" **/
   i = 0;
   do
    {
      nc = fscanf (in2, "%c", &ch);
      if (nc > 0 && char2AAmem[(int) ch] >= 0)
	{
	  s2[++i] = char2AAmem[(int) ch];
	}
    }
   while (nc > 0 && (i < N));

   dim2 = i;	/* We store the length of the 2nd sequence analyzed in dim2 */
   fclose (in2);
   //assert(s2[0] >= N);
   if (dim2<N)
     printf("WARNING: Sequence 2 read from file %s has %d values (<N=%d)\n",
	    argv[2], dim2, N);

   printf("Length of sequence 2 analyzed=%d\n", dim2);

  } /* End of work specific for the master. */

  /* PARALLEL CODE HERE */
  MPI_Bcast( &dim1 , 1    , MPI_INT, 0, MPI_COMM_WORLD ); /* Broadcast dim1 */
  MPI_Bcast( &dim2 , 1    , MPI_INT, 0, MPI_COMM_WORLD ); /* Broadcast dim2 */
  MPI_Bcast( &BS   , 1    , MPI_INT, 0, MPI_COMM_WORLD ); /* Broadcast BS */
  MPI_Bcast( &DELTA, 1    , MPI_INT, 0, MPI_COMM_WORLD ); /* Broadcast DELTA */
  MPI_Bcast( &sim  , AA*AA, MPI_INT, 0, MPI_COMM_WORLD  ); /* Broadcast matrix sim */

 /*int rowsTotal, int mpiRank, int mpiSize*/
  nrows = getRowCount( nrows, rank, size  );
  printf ("rank:%d  nrows=%d \n", rank, nrows );
  /* Allocate the local vector  s1l used to keep the 1st sequence:
     it must have nrows+1 elements of data type short. */
  CHECK_NULL ((s1l = (short *) malloc ( ... ) ));

  if (rank) {  /* Have process other than the master allocate space */
   //CHECK_NULL ((s1 = (short *) malloc (sizeof (short) * (dim1 + 1))));
   CHECK_NULL ((s2 = (short *) malloc (sizeof (short) * (dim2 + 1))));
  }

  //MPI_Bcast( s1+1, dim1, MPI_SHORT, 0, MPI_COMM_WORLD );
  MPI_Scatter( s1+1, nrows, ... ); 
  /* s1+1: use pointer arithmetic to point to the second element in the array
   *       since the 1st position is not used.
   *        ___________________      ____
   *       |      |     |                | 
   *       |s1[0] |s1[1]|       ....     | 
   *       |______|_____|_____       ____| 
   *        ^      ^
   *        |      |
   *        s1    s1+1
   *
   *       */
  MPI_Bcast( s2+1, dim2, ... );


#if 0
  CHECK_NULL ((s1out = (char *) malloc (sizeof (char) * 2 * dim1)));
  CHECK_NULL ((s2out = (char *) malloc (sizeof (char) * 2 * dim2)));
#endif

  int dim1sz;
  /* Allocate full matrix with dim1 rows in the root process, in preparation for traceback,
   * but only a part (nrows+1 rows) for all the other processes. */
  dim1sz = (rank ? nrows : ...); /* condition ? true_case : false_case */
  CHECK_NULL ((hptr = (int *) malloc (sizeof (int) * (dim1sz + 1) * (dim2 + 1))));
  CHECK_NULL ((h = (int **) malloc (sizeof (int *) * (dim1sz + 1))));
  printf ("rank:%d  dim1sz=%d \n", rank, dim1sz );

  /* Mount h[dim1sz][dim2] */
  for (i = 0; i <= dim1sz; i++)
    h[i] = hptr + i * (dim2 + 1); /* Let h[i] point to the beginnin of row i */

#if 0
  CHECK_NULL ((xTracebackptr =
	       (short *) malloc (sizeof (short) * (dim1sz + 1) * (dim2 + 1))));
  CHECK_NULL ((xTraceback = (short **) malloc (sizeof (short *) * (dim1sz + 1))));
  /* Mount xTraceback[dim1sz][dim2] */
  for (i = 0; i <= dim1sz; i++)
    xTraceback[i] = xTracebackptr + i * (dim2 + 1);

  CHECK_NULL ((yTracebackptr =
	       (short *) malloc (sizeof (short) * (dim1sz + 1) * (dim2 + 1))));
  CHECK_NULL ((yTraceback = (short **) malloc (sizeof (short *) * (dim1sz + 1))));
  /* Mount yTraceback[dim1sz][dim2] */
  for (i = 0; i <= dim1sz; i++)
    yTraceback[i] = yTracebackptr + i * (dim2 + 1);
#endif

  

  /* You may want to delete yTraceback and xTraceback updates on the following      */
  /* process since we are not interested on it. See comments below. It is up to you. */

  /** initialize traceback array **/
  Max = xMax = yMax = 0;
#if 0
  for (i = 0; i <= dim1sz; i++)        /* TODO: Check if nrows is ok for all XXXXXXXXXXX */
    for (j = 0; j <= dim2; j++)
      {
	xTraceback[i][j] = -1;
	yTraceback[i][j] = -1;
      }

#endif
  STOP_COUNT_TIME("Initialization time in seconds");

  START_COUNT_TIME;

  /** compute "h" local similarity array **/
  for (i = 0; i <= nrows; i++)
    h[i][0] = 0;
  for (j = 0; j <= dim2; j++)
    h[0][j] = 0;

  /* PARALLEL CODE */
  for (int jj = 1; jj <= dim2; jj += ...) {   /* Strip mining: Define blocks of size BS in the columns */

   if (rank != 0) MPI_Recv( ... );

   for (i = 1; i <= ... ; i++)
    for (j = jj; j < min(...,...); j++) /* Strip mining: Traverse BS columns within the current block */
      {
	diag  = h[i - 1][j - 1] + sim[ s1l[i] ][ s2[j] ];
	down  = h[i - 1][j    ] + DELTA;
	right = h[i    ][j - 1] + DELTA;
	max = MAX3 (diag, down, right);
	if (max <= 0)
	  {
	    h[i][j] = 0;
	    //xTraceback[i][j] = -1;
	    //yTraceback[i][j] = -1;
	    // these values already -1
	  }
	else if (max == diag)
	  {
	    h[i][j] = diag;
	    //xTraceback[i][j] = i - 1;
	    //yTraceback[i][j] = j - 1;
	  }
	else if (max == down)
	  {
	    h[i][j] = down;
	    //xTraceback[i][j] = i - 1;
	    //yTraceback[i][j] = j;
	  }
	else
	  {
	    h[i][j] = right;
	    //xTraceback[i][j] = i;
	    //yTraceback[i][j] = j - 1;
	  }
	if (max > Max)
	  {
	    Max = max;
	    xMax = i;
	    yMax = j;
	  }
      }	// end for loop
   if (rank != nprocs-1 ) MPI_Send( ... );
   /* Here we could use a non-blocking send. */
  }

#if 1
  /* PARALLEL CODE */
    int localres[2];
    int globalres[2];
    localres[0] = Max;
    localres[1] = rank;

    MPI_Allreduce(localres, ... );

    if (rank == 0) {
        printf("Rank %d has maximum value of %d\n",
               globalres[1], globalres[0]);
    }


#endif


  STOP_COUNT_TIME("Computation of scoring matrix time in seconds");
  /* PARALLEL CODE */
  int rem   = dim1 % nprocs; /* If nprocs does not divide dim1 exactly one extra row is given to the first rem processes
                                i.e. those having rank < rem. */
  int extra = (rank < rem) ? rank : rem; /* How many extra rows do we need to account for? */
  int grow  = (dim1/nprocs)*rank + extra +xMax; /* Global index between 1 and dim1 for xMax. */
  printf ("rank:%d  Max=%d xMax=%d (global xMax=%d) yMax=%d\n", rank, Max, xMax, grow, yMax);
  if (rank==nprocs-1)
  printf ("Last element in matrix: rank=%d\t(local index)\t(global index)\n\t\t\t\th[%d][%d] =\th[%d][%d] = %d\n",
          rank, nrows, dim2, grow, dim2, h[nrows][dim2]);

#if 0
  START_COUNT_TIME;

  /* Parallelization STOPS here. We are not interested on the parallelization  
     of the traceback process, and the generation of the match result. */



  // Reset to max point to do alignment
  i = xMax;
  j = yMax;
  x = y = 0;
  topskip = bottomskip = 1;
  while (i > 0 && j > 0 && h[i][j] > 0)
    {
      if (topskip && bottomskip)
	{
	  s1out[x++] = AA2charmem[s1[i]];
	  s2out[y++] = AA2charmem[s2[j]];
	}
      else if (topskip)
	{
	  s1out[x++] = '-';
	  s2out[y++] = AA2charmem[s2[j]];
	}
      else if (bottomskip)
	{
	  s1out[x++] = AA2charmem[s1[i]];
	  s2out[y++] = '-';
	}
      topskip = (j > yTraceback[i][j]);
      bottomskip = (i > xTraceback[i][j]);
      tempi = i;
      tempj = j;
      i = xTraceback[tempi][tempj];
      j = yTraceback[tempi][tempj];
    }



  // print alignment
  printf ("\n");
  printf ("\n");
  for (i = x - 1; i >= 0; i--)
    printf ("%c", s1out[i]);
  printf ("\n");
  for (j = y - 1; j >= 0; j--)
    printf ("%c", s2out[j]);
  printf ("\n");
  printf ("\n");

  STOP_COUNT_TIME("Traceback and print results time in seconds");
#endif

  if ( ... ) free(s1); /* Only the master process allocated it */
  free(s2); free(h); free(hptr); 
  free(s1l);
  //free(s1out); free(s2out);
  //free(xTracebackptr); free(yTracebackptr); 
  //free(xTraceback); free(yTraceback); 

  /* PARALLEL CODE */
  //i=
  MPI_Finalize ();
  //printf ("rank:%d  err=%d \n", rank, i );

}

void
error (char *s)
{
  fprintf (stderr, "%s\n", s);
  exit (1);
}


void
initChar2AATranslation (void)
{
  int i;
  for (i = 0; i < 256; i++)
    char2AAmem[i] = -1;
  char2AAmem['c'] = char2AAmem['C'] = 0;
  AA2charmem[0] = 'c';
  char2AAmem['g'] = char2AAmem['G'] = 1;
  AA2charmem[1] = 'g';
  char2AAmem['p'] = char2AAmem['P'] = 2;
  AA2charmem[2] = 'p';
  char2AAmem['s'] = char2AAmem['S'] = 3;
  AA2charmem[3] = 's';
  char2AAmem['a'] = char2AAmem['A'] = 4;
  AA2charmem[4] = 'a';
  char2AAmem['t'] = char2AAmem['T'] = 5;
  AA2charmem[5] = 't';
  char2AAmem['d'] = char2AAmem['D'] = 6;
  AA2charmem[6] = 'd';
  char2AAmem['e'] = char2AAmem['E'] = 7;
  AA2charmem[7] = 'e';
  char2AAmem['n'] = char2AAmem['N'] = 8;
  AA2charmem[8] = 'n';
  char2AAmem['q'] = char2AAmem['Q'] = 9;
  AA2charmem[9] = 'q';
  char2AAmem['h'] = char2AAmem['H'] = 10;
  AA2charmem[10] = 'h';
  char2AAmem['k'] = char2AAmem['K'] = 11;
  AA2charmem[11] = 'k';
  char2AAmem['r'] = char2AAmem['R'] = 12;
  AA2charmem[12] = 'r';
  char2AAmem['v'] = char2AAmem['V'] = 13;
  AA2charmem[13] = 'v';
  char2AAmem['m'] = char2AAmem['M'] = 14;
  AA2charmem[14] = 'm';
  char2AAmem['i'] = char2AAmem['I'] = 15;
  AA2charmem[15] = 'i';
  char2AAmem['l'] = char2AAmem['L'] = 16;
  AA2charmem[16] = 'l';
  char2AAmem['f'] = char2AAmem['F'] = 17;
  AA2charmem[17] = 'L';
  char2AAmem['y'] = char2AAmem['Y'] = 18;
  AA2charmem[18] = 'y';
  char2AAmem['w'] = char2AAmem['W'] = 19;
  AA2charmem[19] = 'w';
}

