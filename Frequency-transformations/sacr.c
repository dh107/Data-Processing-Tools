/*
 *  sacr.c
 *
 *  Revised: 14 August 2007
 *
 *  Purpose: to read SAC binary files and write ASCII data files
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* limit number samples */
#define NM 55000                          /*------------------- Not sure the definite meaning of NM ?*/
int main()
{
   FILE *log, *out;                       /*------------------- log is the parameter file created by the csh*/
   float st[NM];                          /*------------------- st is the value of data of each points*/
   float bt,dt,t;
   int nm = NM;
   int e,i,j,ns,nt;
   char x[26],ur[30],uw[30];              /*------------------- what's x, ur and uw? */
/* determine number sequences */
   log = fopen("log","r");                /*------------------- open the parameter file, log */
   fscanf(log,"%2d\n",&ns);               /*------------------- ns is the number of SAC files that will be converted to txt*/
   if (ns == 0)                           /*------------------- if the number of SAC files is zero, then no need to convert*/
     goto done;
/* input file names */
   for (j = 1; j <= ns; j++)              /*------------------- the outside recycle, for all those files */
   {
     fscanf(log,"%26s.sac\n",x);
     strcpy(ur,x);                        /*------------------- copy the string "x" to "ur" */  
//     strcat(ur,".sac");
     strcpy(uw,x);                        /*------------------- copy the string "x" to "uw" */
     strcat(uw,".txt");                   /*------------------- these five lines are the operation to the file names*/
/*   invoke subroutines from sacio.a library */
     rsac1(ur,st,&nt,&bt,&dt,&nm,&e,30);  /*------------------- the function provided by sac, for change the type of data*/
/*   output variables */
     printf("%30s  %.5d  %6.4f\n",uw,nt,dt);   /*-------------- "uw" is the file name */
                                               /*-------------- "nt" should be the npts in one certain file*/
                                               /*-------------- "dt" is the delta of nearest data points*/
/*   output sequences */
     out = fopen(uw,"w");
     for (i = 0; i < nt; i++)                  /*-------------- one question: how to know the exact value of "nt"? */
     {
       t = bt + i * dt;                        /*-------------- "bt" is the begin-time, thus "t" refers time */
//       fprintf(out,"%13.6E %13.6E\n",t,st[i]);
       fprintf(out,"%f %f\n",t,st[i]);
     }
     close(out);
   }
/* finish */
   done:
   close(log);
   return(0);
}
