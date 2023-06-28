#include<sys/file.h>
#include<unistd.h>
#include<stdlib.h>
#include<stdio.h>

/*convert SAC binary to headerless binary file (fromHelm.c generated) */
/*or Helmberger ascii */
int main(ac,av)
 int ac;
 char **av;
 {
 int i, npts, ihd[40],fd1, fd2, fd3, fd4, fd5;
 float dt, fhd[70], *tr1, *tr2, *tr3, *tr4, *tr5;
 char chd[8][24];

fd1=open("thf",O_RDONLY,0644);
fd2=open("rhf",O_RDONLY,0644);
fd3=open("zhf",O_RDONLY,0644);
fd4=open("rvf",O_RDONLY,0644);
fd5=open("zvf",O_RDONLY,0644);


read(fd1,fhd,70*4);  /*Read Sac Float Field*/
read(fd1,ihd,40*4);  /*Read Sac Int   Field*/
read(fd1,chd,24*8);  /*Read Sac Char. Field*/
npts=ihd[9];
dt=fhd[0];
fprintf(stderr,"npts=%d   dt=%f\n",npts,dt);
tr1=(float *)malloc(sizeof(float)*npts);
read(fd1,tr1,npts*sizeof(float));
close(fd1);

read(fd2,fhd,70*4);  /*Read Sac Float Field*/
read(fd2,ihd,40*4);  /*Read Sac Int   Field*/
read(fd2,chd,24*8);  /*Read Sac Char. Field*/
tr2=(float *)malloc(sizeof(float)*npts);
read(fd2,tr2,npts*sizeof(float));
close(fd2);

read(fd3,fhd,70*4);  /*Read Sac Float Field*/
read(fd3,ihd,40*4);  /*Read Sac Int   Field*/
read(fd3,chd,24*8);  /*Read Sac Char. Field*/
tr3=(float *)malloc(sizeof(float)*npts);
read(fd3,tr3,npts*sizeof(float));
close(fd3);

read(fd4,fhd,70*4);  /*Read Sac Float Field*/
read(fd4,ihd,40*4);  /*Read Sac Int   Field*/
read(fd4,chd,24*8);  /*Read Sac Char. Field*/
tr4=(float *)malloc(sizeof(float)*npts);
read(fd4,tr4,npts*sizeof(float));
close(fd4);

read(fd5,fhd,70*4);  /*Read Sac Float Field*/
read(fd5,ihd,40*4);  /*Read Sac Int   Field*/
read(fd5,chd,24*8);  /*Read Sac Char. Field*/
tr5=(float *)malloc(sizeof(float)*npts);
read(fd5,tr5,npts*sizeof(float));
close(fd5);

fprintf(stdout,"%% %d  %f\n",npts,dt);
for(i=0; i < npts; i++)
   {
   fprintf(stdout,"%+12.6E %+12.6E %+12.6E %+12.6E %+12.6E\n",tr1[i],tr2[i],tr3[i],tr4[i],tr5[i]); 
   }

free(tr1);
free(tr2);
free(tr3);
free(tr4);
free(tr5);
}


