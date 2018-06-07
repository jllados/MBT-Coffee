#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"


/*******************************************************************************/
/*                myers and Miller                                             */
/*                                                                             */
/*	makes DP between the the ns[0] sequences and the ns[1]                 */
/*                                                                             */
/*	for MODE, see the function  get_dp_cost                                */
/*******************************************************************************/

#define gap(k)  ((k) <= 0 ? 0 : CL->gop*SCORE_K+CL->gep*SCORE_K*(k))    /* k-symbol indel cost */
static int *sapp;                /* Current script append ptr */
static int  last;                /* Last script op appended */
static int  printPtr;
                       /* Append "Delete k" op */

static int *sns;
static int **sls;

void progDel(int k)
{

    if (last < 0)
    {
        last = sapp[printPtr - 1] -= k;
    }
    else
    {
        last = sapp[printPtr++] =  - (k);
    }
}

void progAdd(int k)
{

    if (last < 0)
    {
        sapp[printPtr - 1] = k;
        sapp[printPtr++] = last;
    }
    else
    {
        last = sapp[printPtr++] = k;
    }
}

void progAlign()
{

    sapp[printPtr++] = last = 0;
}


int openPenalty1(int i, int j, Constraint_list *CL) {
	int gapCol = 32;
    int g;

    if (i == 0 || i == CL->lenprfsub){
        return (0);
    }

    g = CL->prfstd[j][gapCol] +  CL->prfsub[i][gapCol];
		
    return (g);
}

int openPenalty2(int i, int j, Constraint_list *CL) {
	int gapCol = 32;
    int g;

    if (j == 0 || j == CL->lenprfstd) {
        return (0);
    }

    g = CL->prfsub[i][gapCol] + CL->prfstd[j][gapCol];

    return (g);
}

int extPenalty1(int i, int j, Constraint_list *CL) {
	int lenCol = 33;
    int h;

    if (i == 0 || i == CL->lenprfsub) {
        return (0);
    }

    h = CL->prfstd[j][lenCol];

    return (h);
}

int extPenalty2(int i, int j, Constraint_list *CL) {
    int lenCol = 33;
    int h;

    if (j == 0 || j == CL->lenprfstd) {
        return (0);
    }

    h = CL->prfsub[i][lenCol];

    return (h);
}

int gapPenalty1(int i, int j, int k, Constraint_list *CL) {
	int gapCol = 32;
	int lenCol = 33;
    int ix;
    int gp;
    int g, h = 0;

    if (k <= 0)
    {
        return (0);
    }
    if (i == 0 || i == CL->lenprfsub) {
        return (0);
    }

    g = CL->prfstd[j][gapCol] + CL->prfsub[i][gapCol];
    for (ix = 0; ix < k && ix + j < CL->lenprfstd; ix++) {
        h += CL->prfstd[ix + j][lenCol];
    }

    gp = g + h;

    return (gp);
}

int gapPenalty2(int i, int j, int k, Constraint_list *CL) {
    int gapCol = 32;
	int lenCol = 33;
	int ix;
    int gp;
    int g, h = 0;

    if (k <= 0)
    {
        return (0);
    }
    if (j == 0 || j == CL->lenprfstd)
    {
        return (0);
    }

    g = CL->prfsub[i][gapCol] + CL->prfstd[j][gapCol];
    for (ix = 0; ix < k && ix + i < CL->lenprfsub; ix++)
    {
        h += CL->prfsub[ix + i][lenCol];
    }

    gp = g + h;

    return (gp);
}

int progTracepath(int *alnPath1, int *alnPath2) {
    int i, j, k, pos, toDo;
    int alignLen;
    pos = 0;

    toDo = printPtr - 1;

    for (i = 1; i <= toDo; ++i) {
        if (sapp[i] == 0) {
            alnPath1[pos] = 2;
            alnPath2[pos] = 2;
            ++pos;
        } else {
            if ((k = sapp[i]) > 0) {
                for (j = 0; j <= k - 1; ++j) {
                    alnPath2[pos + j] = 2;
                    alnPath1[pos + j] = 1;
                }
                pos += k;
            } else {
                k = (sapp[i] < 0) ? sapp[i] * - 1: sapp[i];
                for (j = 0; j <= k - 1; ++j) {
                    alnPath1[pos + j] = 2;
                    alnPath2[pos + j] = 1;
                }
                pos += k;
            }
        }
    }
	
	alignLen = pos;
    return alignLen;
}

void addGGaps(Alignment *A, int alignmentLength, int *ns, int **ls, int *alnPath1, int *alnPath2) {
    int j;
    int i, ix;
    int len, l1, l2;
	const int extraEndElemNum = 0;
	char ** char_buf;
	
	A=realloc_aln2  ( A,A->max_n_seq,alignmentLength+1);
	char_buf=declare_char (A->max_n_seq,alignmentLength+1);
	
	l1=strlen (A->seq_al[ls[0][0]]);
	l2=strlen (A->seq_al[ls[1][0]]);

    for (j = 0; j < ns[0]; j++) {
        ix = 0;
        for (i = 0; i < alignmentLength; i++) {
            if (alnPath1[i] == 2) {
                if (ix < (l1 - extraEndElemNum)) {
                    char_buf[ls[0][j]][i] = A->seq_al[ls[0][j]][ix];
                } else {
                    char_buf[ls[0][j]][i] = '\0';
                }
                ix++;
            } else if (alnPath1[i] == 1) {
                char_buf[ls[0][j]][i] = '-'; //insertion in first alignment...
            } else
            {
                //cerr << "Error in aln_path\n";
            }
        }
        char_buf[ls[0][j]][i] = '\0';

        len = alignmentLength;
    }

    for (j = 0; j < ns[1]; j++) {
        ix = 0;
        for (i = 0; i < alignmentLength; i++) {
            if (alnPath2[i] == 2) {
                if (ix < (l2 - extraEndElemNum)) {
                    char_buf[ls[1][j]][i] = A->seq_al[ls[1][j]][ix];
                } else {
                    char_buf[ls[1][j]][i] = '\0';
                }
                ix++;
            } else if (alnPath2[i] == 1) {
                char_buf[ls[1][j]][i] = '-'; //insertion in second alignment...
            } else {
                //cerr << "Error in alnPath\n";
            }
        }
        char_buf[ls[1][j]][i] = '\0';

        len = alignmentLength;
    }
	
	A->len_aln=len;
	A->nseq=ns[0]+ns[1];
	
	for ( i=0; i< ns[0]; i++){char_buf[ls[0][i]][len]='\0'; sprintf ( A->seq_al[ls[0][i]], "%s", char_buf[ls[0][i]]);}
	for ( i=0; i< ns[1]; i++){char_buf[ls[1][i]][len]='\0'; sprintf ( A->seq_al[ls[1][i]], "%s", char_buf[ls[1][i]]);}
	
	free_char ( char_buf, -1);
	
	//for ( i=0; i< ns[0]; i++){for (ix = 0; ix < alignmentLength; ix++) {fprintf(stdout, "%c ",A->seq_al[ls[0][i]][ix]);}fprintf(stdout, "\n");fprintf(stdout, "\n");}
	
	//for ( i=0; i< ns[1]; i++){for (ix = 0; ix < alignmentLength; ix++) {fprintf(stdout, "%c ",A->seq_al[ls[1][i]][ix]);}fprintf(stdout, "\n");fprintf(stdout, "\n");}

}

int myers_miller_pair_wise (Alignment *A,int *ns, int **ls,Constraint_list *CL )
	{
	int **pos;
	int a,b, i, j, l,l1, l2, len;
	int *S;
	char ** char_buf;
	int score;
		
	/********Prepare Penalties******/
	//ns2master_ns (ns,ls, &sns,&sls);
	sns=ns;
	sls=ls;

	/********************************/
	
	//fprintf(stdout,"(INFO) Start: myers_miller_pair_wise\n");
	
	pos=aln2pos_simple ( A,-1, ns, ls);
	

	l1=strlen (A->seq_al[ls[0][0]]);
	l2=strlen (A->seq_al[ls[1][0]]);
	S=(int*)vcalloc (l1+l2+1, sizeof (int));
	last=0;
	printPtr = 1;
	sapp=S;
	
	//fprintf(stdout,"(INFO) Start: diff\n");
	
	score=diff (A,ns, ls, 0, l1, 0, l2, 0, 0, CL, pos);	
	diff (NULL,ns, ls, 0, l1, 0, l2, 0, 0, CL, pos);
	
	//fprintf(stdout,"(INFO) End: diff %d\n", score);
	
	/*fprintf(stdout, "---- SAPP -- %d --\n", l1+l2+1);
	for(i=0;i<l1+l2+1;i++)
		fprintf(stdout, "%d ",sapp[i]);
	fprintf(stdout, "\n");
	fprintf(stdout, "---- SAPPEND ----\n");*/
	
	int alnPath1[l1+l2+1], alnPath2[l1+l2+1];
	
	memset(alnPath1, 0, sizeof(alnPath1));
	memset(alnPath2, 0, sizeof(alnPath2));
	
	int alignmentLength = progTracepath(alnPath1, alnPath2);

	/*fprintf(stdout, "---- alnPath1 ----\n");
	for(i=0;i<l1+l2+1;i++)
		fprintf(stdout, "%d ",alnPath1[i]);
	fprintf(stdout, "\n");
	fprintf(stdout, "---- alnPath1END ----\n");
	
	fprintf(stdout, "---- alnPath2 ----\n");
	for(i=0;i<l1+l2+1;i++)
		fprintf(stdout, "%d ",alnPath2[i]);
	fprintf(stdout, "\n");
	fprintf(stdout, "---- alnPath2END ----\n");*/
	

	addGGaps(A, alignmentLength, ns, ls, alnPath1, alnPath2);
	
	//fprintf(stdout, "---- addGGapsEND ----\n");
	
	vfree (S);
	l1=strlen (A->seq_al[ls[0][0]]);
	l2=strlen (A->seq_al[ls[1][0]]);
	
	if ( l1!=l2) exit(1);
		
	free_int (pos, -1);
	
	//fprintf(stdout, "------------ END FUNC %d ------------\n\n\n\n\n", score);
	
	return score;
	}


int diff (Alignment *A, int *ns, int **ls, int s1, int M,int s2, int N , int go1, int go2, Constraint_list *CL, int **pos) {
	static int *HH;
	static int *DD;
	/* Forward cost-only vectors */
	static int *RR;
	static int *SS;
	static int *gS;
    /* Reverse cost-only vectors */
    int midi, midj, type;
    /* Midpoint, type, and cost */
    int midh;

	/*TREATMENT OF THE TERMINAL GAP PENALTIES*/
	/*TG_MODE=0---> gop and gep*/
	/*TG_MODE=1---> ---     gep*/
	/*TG_MODE=2---> ---     ---*/

	 if ( !HH) {
		int L;
	    L=M+N+1;

	    HH=(int*)vcalloc (L, sizeof (int));
	    DD=(int*)vcalloc (L, sizeof (int));
	    RR=(int*)vcalloc (L, sizeof (int));
	    SS=(int*)vcalloc (L, sizeof (int));
	    gS=(int*)vcalloc (L, sizeof (int));
	}

	 if ( A==NULL) {
	    vfree(HH);
	    vfree(DD);
	    vfree(RR);
	    vfree(SS);
	    vfree(gS);
	    HH=DD=RR=SS=gS=NULL;
	    return 0;
	}
	 
{
	int i, j;
	int hh, f, e, s;
    int t, tl, g, h;
	
	if (N <= 0) {
		if (M > 0) 
			progDel(M);
			return ( - gapPenalty1(s1, s2, M, CL));
	}
	
	if (M <= 1) { 
		if (M <= 0) {
			progAdd(N);
		  	return ( - gapPenalty2(s1, s2, N, CL));
		}

	    if (go1 == 0) {
	    	midh =  - gapPenalty1(s1 + 1, s2 + 1, N, CL);
	    } else {
	    	midh =  - gapPenalty2(s1 + 1, s2, 1, CL) - gapPenalty1(s1 + 1, s2 + 1, N, CL);
	    }
	    midj = 0;
    
		for (j = 1; j <= N; j++) { 
			hh = - gapPenalty1(s1, s2 + 1, j - 1, CL) +(CL->get_dp_cost) (A, pos, sns[0], sls[0],s1, pos, sns[1], sls[1],j-1+s2,CL) - gapPenalty1(s1 + 1, s2 + j + 1, N - j, CL);
			if (hh > midh) {
				midh = hh;
			 	midj = j;
		    }
		}
		
	    if (midj == 0) {
			progAdd(N);
	    	progDel(1);
	    } else {
	    	if (midj > 1) 
	    		progAdd(midj-1);
		  	progAlign();
		  	if (midj < N) 
		  		progAdd(N-midj);
		}
	       
	    return midh;
	}

/* Divide: Find optimum midpoint (midi,midj) of cost midh */
 
	midi = M/2;            /* Forward phase:                          */
	HH[0] = 0;            /*   Compute C(M/2,k) & D(M/2,k) for all k */
	
	t =  - openPenalty1(s1, s2 + 1, CL);
    tl =  - extPenalty1(s1, s2 + 1, CL);
	
	for (j = 1; j <= N; j++) {
	    HH[j] = t = t + tl;
        DD[j] = t - openPenalty2(s1 + 1, s2 + j, CL);
	}
	
	if(go1 == 0) {
		t = 0;
	} else {
		t =  - openPenalty2(s1 + 1, s2, CL);
	}
	tl =  - extPenalty2(s1 + 1, s2, CL);
	
	for (i = 1; i <= midi; i++) { 
		s = HH[0];
	    HH[0] = hh = t = t + tl;
	    f = t - openPenalty1(s1 + i, s2 + 1, CL);

		for (j = 1; j <= N; j++) {
			
			g = openPenalty1(s1 + i, s2 + j, CL);
            h = extPenalty1(s1 + i, s2 + j, CL);
			
			if ((hh =   hh - g - h) > (f =   f   - h)) {
				f = hh;
			}
				
			g = openPenalty2(s1 + i, s2 + j, CL);
            h = extPenalty2(s1 + i, s2 + j, CL);
		 	if ((hh = HH[j] - g - h) > (e = DD[j] - h)) { 
		 		e = hh;
		 	}
		 	
		 	hh = s + (CL->get_dp_cost) (A, pos, sns[0], sls[0],i-1+s1, pos, sns[1], sls[1],j-1+s2,CL);
		 	if (f > hh) hh = f;
		 	if (e > hh) hh = e;
		 
		     
		 	s = HH[j];
		 	HH[j] = hh;
		 	DD[j] = e;
		}
	}
	
	DD[0] = HH[0];
	RR[N] = 0;            /* Reverse phase:                          */
	tl = 0;
	 
	for (j = N-1; j >= 0; j--) {
		g =  - openPenalty1(s1 + M, s2 + j + 1, CL);
        tl -= extPenalty1(s1 + M, s2 + j + 1, CL);
        RR[j] = g + tl;
        SS[j] = RR[j] - openPenalty2(s1 + M, s2 + j, CL);
        gS[j] = openPenalty2(s1 + M, s2 + j, CL);
	}
	
	tl = 0;
	for (i = M-1; i >= midi; i--) {
		s = RR[N];
		
		if (go2 == 0) {
            g = 0;
        } else {
            g =  - openPenalty2(s1 + i + 1, s2 + N, CL);
        }
        tl -= extPenalty2(s1 + i + 1, s2 + N, CL);
		
	    RR[N] = hh = g + tl;
	    t = openPenalty1(s1 + i, s2 + N, CL);
	    f = RR[N] - t;
	     
	    for (j = N-1; j >= 0; j--) { 
	    	g = openPenalty1(s1 + i, s2 + j + 1, CL);
            h = extPenalty1(s1 + i, s2 + j + 1, CL);
	    	
		  	if ((hh =   hh  - g - h) > (f = f - h - g + t)) f = hh;
		  	
		  	t = g;
            g = openPenalty2(s1 + i + 1, s2 + j, CL);
            h = extPenalty2(s1 + i + 1, s2 + j, CL);
            hh = RR[j] - g - h;
		  	if (i == (M - 1)) e = SS[j] - h;
			else { e = SS[j] - h - g + openPenalty2(s1 + i + 2, s2 + j, CL); gS[j] = g; }
		  	if (hh > e) {
                e = hh;
            }
		  	
		  	hh = s + (CL->get_dp_cost) (A, pos, sns[0], sls[0],i+s1, pos, sns[1], sls[1],j+s2,CL);
		  	
		  	if (f > hh) hh = f;
		  	if (e > hh) hh = e;
		 
		  	s = RR[j];
		  	RR[j] = hh;
		  	SS[j] = e;
		}
	}
	SS[N] = RR[N];
	gS[N] = openPenalty2(s1 + midi + 1, s2 + N, CL);
	midh = HH[0]+RR[0];        /* Find optimal midpoint */
	midj = 0;
	type = 1;
	
	for (j = 0; j <= N; j++)
		if ((hh = HH[j] + RR[j]) >= midh)
	    	if (hh > midh || (HH[j] != DD[j] && RR[j] == SS[j])) { 
		    	midh = hh;
		    	midj = j;
		    }
	for (j = N; j >= 0; j--)
		if ((hh = DD[j] + SS[j] + gS[j]) > midh) {
	    	midh = hh;
		 	midj = j;
		 	type = 2;
		}
	}	
	    
/* Conquer: recursively around midpoint */

	if (type == 1) { 
    	diff (A,ns, ls, s1,midi, s2, midj, go1, 1, CL, pos); 
    	diff (A,ns, ls, s1+midi,M-midi, s2+midj, N-midj, 1,go2, CL, pos);
    } else { 
      	diff (A,ns, ls, s1,midi-1, s2, midj, go1,0, CL, pos); 
      	progDel(2);
      	diff (A,ns, ls, s1+midi+1, M-midi-1,s2+midj, N-midj,0,go2, CL, pos); 
    }
  	return midh;
  }













/******************************COPYRIGHT NOTICE*******************************/
/*© Centro de Regulacio Genomica */
/*and */
/*Cedric Notredame */
/*2012-07-12 19:05:45. */
/*All rights reserved.*/
/*This file is part of T-COFFEE.*/
/**/
/*    T-COFFEE is free software; you can redistribute it and/or modify*/
/*    it under the terms of the GNU General Public License as published by*/
/*    the Free Software Foundation; either version 2 of the License, or*/
/*    (at your option) any later version.*/
/**/
/*    T-COFFEE is distributed in the hope that it will be useful,*/
/*    but WITHOUT ANY WARRANTY; without even the implied warranty of*/
/*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the*/
/*    GNU General Public License for more details.*/
/**/
/*    You should have received a copy of the GNU General Public License*/
/*    along with Foobar; if not, write to the Free Software*/
/*    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA*/
/*...............................................                                                                                      |*/
/*  If you need some more information*/
/*  cedric.notredame@europe.com*/
/*...............................................                                                                                                                                     |*/
/**/
/**/
/*	*/
/******************************COPYRIGHT NOTICE*******************************/
