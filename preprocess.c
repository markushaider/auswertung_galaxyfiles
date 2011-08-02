#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>



int main(int argc, char *argv[]) {

  if ( argc != 3 ) {
    printf("Arguments: lpath spath\n");
    return 1;
  }

  FILE *out;
  int nfiles;
  int ngalxs;
  char fname[140]; 

  char * lpath = argv[1];
  char * spath = argv[2];

  /* get the number of timesteps */
  sprintf(fname,"wc -l %s/galaxies",lpath);
  printf("test: %s \n", fname);
  out = popen(fname, "r");
  fscanf(out, "%d", &nfiles);
  printf("Dateien: %d\n", nfiles);
  pclose(out);
  /* get the number of galaxies */
  sprintf(fname,"wc -l %s/gt00001-*.06",lpath);
  out = popen(fname, "r");
  fscanf(out, "%d", &ngalxs);
  printf("Galaxien: %d\n", ngalxs);
  pclose(out);

  /* Counters */
  int gfn = 0, ngal = 0, k = 0, i=0;

  /* Array für Glaxy-Dateinamen und Zeit */
  char **fnames = (char **)malloc(sizeof(char*)*nfiles);
  float *time  = (float *)malloc(sizeof(float)*nfiles);
  float *z  = (float *)malloc(sizeof(float)*nfiles);
  float *a  = (float *)malloc(sizeof(float)*nfiles);


  /* Öffnen von "galaxies" welches Dateinamen und Zeit beinhaltet */
  sprintf(fname,"%s/galaxies",lpath);
  FILE *gfile = fopen(fname,"r");
  if(gfile == NULL)
  {
    printf("Cannot open galaxy file\n");
    exit(-1);
  }

  /* Einlesen der Dateinamen und der Zeit */  

  for (gfn=0; gfn<nfiles; gfn++) {
    fnames[gfn] = (char*)malloc(sizeof(char)*15);
    fscanf(gfile,"%s %f %f %*f %*f %*i\n",fnames[gfn],&time[gfn], &a[gfn]);
    /* printf("gfn %i\n",gfn); */
  }
  fclose(gfile);
  for (gfn =0; gfn<nfiles; gfn++) {
    z[gfn]=1.0/a[gfn]-1.0;
    printf("name %s, time %f, z %f\n",fnames[gfn],time[gfn], z[gfn]);
   }

  /* get number of fileds */
  char zeile[800];
  sprintf(fname,"%s/%s",lpath,fnames[0]);
  gfile = fopen(fname,"r");
  int nfields=0;
  while((zeile[k] = fgetc(gfile))!='\n')
    k++;
  strtok(zeile, " ");
  nfields=1;
  while(strtok(NULL," ")) {
    nfields++;
  }
  printf("Die Datei hat %i Felder\n",nfields);

  float *line  = (float *)malloc(sizeof(float)*nfields);
  float *cluster_stf = (float *)calloc(nfiles,sizeof(float));
  float *cluster_stf_avg = (float *)calloc(nfiles,sizeof(float));
  float *coldmetal_avg = (float *)calloc(nfiles,sizeof(float));
  int *nGalaxies = (int *)calloc(nfiles,sizeof(int));
  float stf_tot=0;

  float **stf_histo = (float **)malloc(nfiles * sizeof(float *));
  for(gfn = 0; gfn < nfiles; gfn++)
      stf_histo[gfn] = (float *)calloc(ngalxs,sizeof(float));
  float **mcold_histo = (float **)malloc(nfiles * sizeof(float *));
  for(gfn = 0; gfn < nfiles; gfn++)
      mcold_histo[gfn] = (float *)calloc(ngalxs, sizeof(float));
  float **rdisk_histo = (float **)malloc(nfiles * sizeof(float *));
  for(gfn = 0; gfn < nfiles; gfn++)
      rdisk_histo[gfn] = (float *)calloc(ngalxs, sizeof(float));
  float **coldmetal_histo = (float **)malloc(nfiles * sizeof(float *));
  for(gfn = 0; gfn < nfiles; gfn++)
      coldmetal_histo[gfn] = (float *)calloc(ngalxs,sizeof(float));
  float **mhalo_histo = (float **)malloc(nfiles * sizeof(float *));
  for(gfn = 0; gfn < nfiles; gfn++)
      mhalo_histo[gfn] = (float *)calloc(ngalxs,sizeof(float));
  float **pos_x = (float **)malloc(nfiles * sizeof(float *));
  for(gfn = 0; gfn < nfiles; gfn++)
      pos_x[gfn] = (float *)calloc(ngalxs, sizeof(float));
  float **pos_y = (float **)malloc(nfiles * sizeof(float *));
  for(gfn = 0; gfn < nfiles; gfn++)
      pos_y[gfn] = (float *)calloc(ngalxs, sizeof(float));
  float **pos_z = (float **)malloc(nfiles * sizeof(float *));
  for(gfn = 0; gfn < nfiles; gfn++)
      pos_z[gfn] = (float *)calloc(ngalxs, sizeof(float));
  float **vel_x = (float **)malloc(nfiles * sizeof(float *));
  for(gfn = 0; gfn < nfiles; gfn++)
      vel_x[gfn] = (float *)calloc(ngalxs, sizeof(float));
  float **vel_y = (float **)malloc(nfiles * sizeof(float *));
  for(gfn = 0; gfn < nfiles; gfn++)
      vel_y[gfn] = (float *)calloc(ngalxs, sizeof(float));
  float **vel_z = (float **)malloc(nfiles * sizeof(float *));
  for(gfn = 0; gfn < nfiles; gfn++)
      vel_z [gfn] = (float *)calloc(ngalxs, sizeof(float ));


  for (gfn = 0; gfn<nfiles; gfn++)
  {
    printf("Datei %i\n",gfn+1);
    printf("fname %s\n",fnames[gfn]);
    sprintf(fname,"%s/%s",lpath,fnames[gfn]);
    gfile = fopen(fname,"r");
    if(gfile == NULL)
    {
      printf("Cannot open file\n");
      exit(-1);
    }


    /* Einlesen der galaxy Dateien */
    for( ngal = 0; ngal<ngalxs; ngal++)
    {

      char * token;
      k=0;
      while(!feof(gfile)&&(zeile[k] = fgetc(gfile))!='\n')
	k++;
      token= strtok(zeile, " ");
      line[0]=atof(token);
      k=1;

      while(token=strtok(NULL," ")) {
	line[k]=atof(token);
	k++;
      }

      /* printf("test Zeile: \n"); */
      /* for(k=0;k<nfields;k++) { */
      /* 	printf("%f ",line[k]); */
      /* } */
      /* printf("\n"); */

      /* for (k = 0;k<nfields; k++) */
      /* { */
      /* 	fscanf(gfile,"%f",&line[k]); */
      /* } */
      /* fscanf(gfile,"\n"); */
      
      /* if(line[0]>0&&line[2]>1E-3) */
      /* if (line[46]!=0) { */
      /* better to ask for field 46 if present */
      if(line[0]>0 && line[2]>1E-3) {

	/* Summe aller 5 Starformation Rates */
	stf_tot=line[33]+line[34]+line[35]+line[36]+line[37];
	cluster_stf[gfn]+=stf_tot;
	coldmetal_avg[gfn]+=line[19];
        nGalaxies[gfn]++;
	/* Histograms */
	stf_histo[gfn][ngal]=stf_tot;
	mcold_histo[gfn][ngal]=line[2];
	rdisk_histo[gfn][ngal]=line[0];
	coldmetal_histo[gfn][ngal]=line[19];
	mhalo_histo[gfn][ngal]=line[8];
	pos_x[gfn][ngal]=line[25];
	pos_y[gfn][ngal]=line[26];
	pos_z[gfn][ngal]=line[27];
	vel_x[gfn][ngal]=line[28];
	vel_y[gfn][ngal]=line[29];
	vel_z[gfn][ngal]=line[30];
      } else {
	stf_histo[gfn][ngal]=NAN;
	mcold_histo[gfn][ngal]=NAN;
	rdisk_histo[gfn][ngal]=NAN;
	coldmetal_histo[gfn][ngal]=NAN;
	mhalo_histo[gfn][ngal]=NAN;
	pos_x[gfn][ngal]=NAN;
	pos_y[gfn][ngal]=NAN;
	pos_z[gfn][ngal]=NAN;
	vel_x[gfn][ngal]=NAN;
	vel_y[gfn][ngal]=NAN;
	vel_z[gfn][ngal]=NAN;
      }

    }
    cluster_stf_avg[gfn]=cluster_stf[gfn]/nGalaxies[gfn];
    coldmetal_avg[gfn]/=nGalaxies[gfn];

    printf("Anzahl der nichtexistenten Galaxien: %i\n",ngalxs-nGalaxies[gfn]);
    printf("Anzahl der existenten Galaxien: %i\n",nGalaxies[gfn]);
    fclose(gfile);
  }

  free(line);

  for(gfn = 0; gfn < nfiles; gfn++)
    free(fnames[gfn]);
  free(fnames);

  FILE *output;
  sprintf(fname,"%s/info.dat",spath);
  printf("Debug: %s\n",fname);
  printf("%i %i\n",nfiles,ngalxs);
  output=fopen(fname,"w");
  fprintf(output,"#Timesteps Galaxies\n");
  fprintf(output,"%i %i\n",nfiles,ngalxs);
  fclose(output);

  /* Binary Output */
  sprintf(fname,"%s/cluster_stf.dat",spath);
  output=fopen(fname,"wb");
  for (gfn = 0; gfn < nfiles; gfn++) {
    fwrite(&time[gfn],sizeof(float),1,output);
    fwrite(&z[gfn],sizeof(float),1,output);
    fwrite(&cluster_stf[gfn],sizeof(float),1,output);
  }
  fclose(output);

  sprintf(fname,"%s/cluster_stf_avg.dat",spath);
  output=fopen(fname,"wb");
  for (gfn = 0; gfn < nfiles; gfn++) {
    fwrite(&time[gfn],sizeof(float),1,output);
    fwrite(&z[gfn],sizeof(float),1,output);
    fwrite(&cluster_stf_avg[gfn],sizeof(float),1,output);
  }
  fclose(output);

  sprintf(fname,"%s/coldmetal_avg.dat",spath);
  output=fopen(fname,"wb");
  for (gfn = 0; gfn < nfiles; gfn++) {
    fwrite(&time[gfn],sizeof(float),1,output);
    fwrite(&z[gfn],sizeof(float),1,output);
    fwrite(&coldmetal_avg[gfn],sizeof(float),1,output);
  }
  fclose(output);

  sprintf(fname,"%s/ngalaxies.dat",spath);
  output=fopen(fname,"wb");
  for (gfn = 0; gfn < nfiles; gfn++) {
    fwrite(&time[gfn],sizeof(float),1,output);
    fwrite(&z[gfn],sizeof(float),1,output);
    fwrite(&nGalaxies[gfn],sizeof(float),1,output);
  }
  fclose(output);

  sprintf(fname,"%s/stf_histo.dat",spath);
  output=fopen(fname,"wb");
  for (gfn = 0; gfn < nfiles; gfn++) {
    fwrite(&time[gfn],sizeof(float),1,output);
    fwrite(&z[gfn],sizeof(float),1,output);
    for(ngal=0; ngal<ngalxs; ngal++) {
      fwrite(&stf_histo[gfn][ngal],sizeof(float),1,output);
    }
  }
  fclose(output);

  sprintf(fname,"%s/mcold_histo.dat",spath);
  output=fopen(fname,"wb");
  for (gfn = 0; gfn < nfiles; gfn++) {
    fwrite(&time[gfn],sizeof(float),1,output);
    fwrite(&z[gfn],sizeof(float),1,output);
    for(ngal=0; ngal<ngalxs; ngal++) {
      fwrite(&mcold_histo[gfn][ngal],sizeof(float),1,output);
    }
  }
  fclose(output);

  sprintf(fname,"%s/rdisk_histo.dat",spath);
  output=fopen(fname,"wb");
  for (gfn = 0; gfn < nfiles; gfn++) {
    fwrite(&time[gfn],sizeof(float),1,output);
    fwrite(&z[gfn],sizeof(float),1,output);
    for(ngal=0; ngal<ngalxs; ngal++) {
      fwrite(&rdisk_histo[gfn][ngal],sizeof(float),1,output);
    }
  }
  fclose(output);

  sprintf(fname,"%s/coldmetal_histo.dat",spath);
  output=fopen(fname,"wb");
  for (gfn = 0; gfn < nfiles; gfn++) {
    fwrite(&time[gfn],sizeof(float),1,output);
    fwrite(&z[gfn],sizeof(float),1,output);
    for(ngal=0; ngal<ngalxs; ngal++) {
      fwrite(&coldmetal_histo[gfn][ngal],sizeof(float),1,output);
    }
  }
  fclose(output);

  sprintf(fname,"%s/mhalo_histo.dat",spath);
  output=fopen(fname,"wb");
  for (gfn = 0; gfn < nfiles; gfn++) {
    fwrite(&time[gfn],sizeof(float),1,output);
    fwrite(&z[gfn],sizeof(float),1,output);
    for(ngal=0; ngal<ngalxs; ngal++) {
      fwrite(&mhalo_histo[gfn][ngal],sizeof(float),1,output);
    }
  }
  fclose(output);

  sprintf(fname,"%s/positions_x.dat",spath);
  output=fopen(fname,"wb");
  for (gfn = 0; gfn < nfiles; gfn++) {
    fwrite(&time[gfn],sizeof(float),1,output);
    fwrite(&z[gfn],sizeof(float),1,output);
    for(ngal=0; ngal<ngalxs; ngal++) {
      fwrite(&pos_x[gfn][ngal],sizeof(float),1,output);
    }
  }
  fclose(output);

  sprintf(fname,"%s/positions_y.dat",spath);
  output=fopen(fname,"wb");
  for (gfn = 0; gfn < nfiles; gfn++) {
    fwrite(&time[gfn],sizeof(float),1,output);
    fwrite(&z[gfn],sizeof(float),1,output);
    for(ngal=0; ngal<ngalxs; ngal++) {
      fwrite(&pos_y[gfn][ngal],sizeof(float),1,output);
    }
  }
  fclose(output);

  sprintf(fname,"%s/positions_z.dat",spath);
  output=fopen(fname,"wb");
  for (gfn = 0; gfn < nfiles; gfn++) {
    fwrite(&time[gfn],sizeof(float),1,output);
    fwrite(&z[gfn],sizeof(float),1,output);
    for(ngal=0; ngal<ngalxs; ngal++) {
      fwrite(&pos_z[gfn][ngal],sizeof(float),1,output);
    }
  }
  fclose(output);

  sprintf(fname,"%s/velocity_x.dat",spath);
  output=fopen(fname,"wb");
  for (gfn = 0; gfn < nfiles; gfn++) {
    fwrite(&time[gfn],sizeof(float),1,output);
    fwrite(&z[gfn],sizeof(float),1,output);
    for(ngal=0; ngal<ngalxs; ngal++) {
      fwrite(&vel_x[gfn][ngal],sizeof(float),1,output);
    }
  }
  fclose(output);

  sprintf(fname,"%s/velocity_y.dat",spath);
  output=fopen(fname,"wb");
  for (gfn = 0; gfn < nfiles; gfn++) {
    fwrite(&time[gfn],sizeof(float),1,output);
    fwrite(&z[gfn],sizeof(float),1,output);
    for(ngal=0; ngal<ngalxs; ngal++) {
      fwrite(&vel_y[gfn][ngal],sizeof(float),1,output);
    }
  }
  fclose(output);

  sprintf(fname,"%s/velocity_z.dat",spath);
  output=fopen(fname,"wb");
  for (gfn = 0; gfn < nfiles; gfn++) {
    fwrite(&time[gfn],sizeof(float),1,output);
    fwrite(&z[gfn],sizeof(float),1,output);
    for(ngal=0; ngal<ngalxs; ngal++) {
      fwrite(&vel_z[gfn][ngal],sizeof(float),1,output);
    }
  }
  fclose(output);


  /* ASCII Output */
  /* Ausschreiben */
  /* FILE *output; */
  /* sprintf(fname,"%s/cluster_stf.dat",spath); */
  /* output=fopen(fname,"w"); */
  /* for (gfn = 0; gfn < nfiles; gfn++) { */
  /*   fprintf(output,"%f %f %f\n",time[gfn], z[gfn], cluster_stf[gfn]); */
  /* } */

  /* fclose(output); */
  /* sprintf(fname,"%s/cluster_stf_avg.dat",spath); */
  /* output=fopen(fname,"w"); */
  /* for (gfn = 0; gfn < nfiles; gfn++) { */
  /*   fprintf(output,"%f %f %f\n",time[gfn], z[gfn], cluster_stf_avg[gfn]); */
  /* } */
  /* fclose(output); */
  /* sprintf(fname,"%s/coldmetal_avg.dat",spath); */
  /* output=fopen(fname,"w"); */
  /* for (gfn = 0; gfn < nfiles; gfn++) { */
  /*   fprintf(output,"%f %f %f\n",time[gfn], z[gfn], coldmetal_avg[gfn]); */
  /* } */
  /* fclose(output); */
  /* sprintf(fname,"%s/ngalaxies.dat",spath); */
  /* output=fopen(fname,"w"); */
  /* for (gfn = 0; gfn < nfiles; gfn++) { */
  /*   fprintf(output,"%f %f %i\n",time[gfn], z[gfn], nGalaxies[gfn]); */
  /* } */
  /* fclose(output); */


  /* sprintf(fname,"%s/stf_histo.dat",spath); */
  /* output=fopen(fname,"w"); */
  /* for (gfn = 0; gfn < nfiles; gfn++) { */
  /*   fprintf(output,"%f %f",time[gfn],z[gfn]); */
  /*   for(ngal=0; ngal<ngalxs; ngal++) { */
  /*     fprintf(output," %f", stf_histo[gfn][ngal]); */
  /*   } */
  /*   fprintf(output,"\n"); */
  /* } */
  /* fclose(output); */

  /* sprintf(fname,"%s/mcold_histo.dat",spath); */
  /* output=fopen(fname,"w"); */
  /* for (gfn = 0; gfn < nfiles; gfn++) { */
  /*   fprintf(output,"%f %f",time[gfn],z[gfn]); */
  /*   for(ngal=0; ngal<ngalxs; ngal++) { */
  /*     fprintf(output," %f", mcold_histo[gfn][ngal]); */
  /*   } */
  /*   fprintf(output,"\n"); */
  /* } */
  /* fclose(output); */

  /* sprintf(fname,"%s/rdisk_histo.dat",spath); */
  /* output=fopen(fname,"w"); */
  /* for (gfn = 0; gfn < nfiles; gfn++) { */
  /*   fprintf(output,"%f %f",time[gfn],z[gfn]); */
  /*   for(ngal=0; ngal<ngalxs; ngal++) { */
  /*     fprintf(output," %f", rdisk_histo[gfn][ngal]); */
  /*   } */
  /*   fprintf(output,"\n"); */
  /* } */
  /* fclose(output); */

  /* sprintf(fname,"%s/coldmetal_histo.dat",spath); */
  /* output=fopen(fname,"w"); */
  /* for (gfn = 0; gfn < nfiles; gfn++) { */
  /*   fprintf(output,"%f %f",time[gfn],z[gfn]); */
  /*   for(ngal=0; ngal<ngalxs; ngal++) { */
  /*     fprintf(output," %f", coldmetal_histo[gfn][ngal]); */
  /*   } */
  /*   fprintf(output,"\n"); */
  /* } */
  /* fclose(output); */

  /* sprintf(fname,"%s/mhalo_histo.dat",spath); */
  /* output=fopen(fname,"w"); */
  /* for (gfn = 0; gfn < nfiles; gfn++) { */
  /*   fprintf(output,"%f %f",time[gfn],z[gfn]); */
  /*   for(ngal=0; ngal<ngalxs; ngal++) { */
  /*     fprintf(output," %f", mhalo_histo[gfn][ngal]); */
  /*   } */
  /*   fprintf(output,"\n"); */
  /* } */
  /* fclose(output); */

  /* sprintf(fname,"%s/positions_ascii_x.dat",spath); */
  /* output=fopen(fname,"w"); */
  /* for (gfn = 0; gfn < nfiles; gfn++) { */
  /*   fprintf(output,"%f %f",time[gfn],z[gfn]); */
  /*   for(ngal=0; ngal<ngalxs; ngal++) { */
  /*     fprintf(output," %f", pos_x[gfn][ngal]); */
  /*   } */
  /*   fprintf(output,"\n"); */
  /* } */
  /* fclose(output); */

  /* sprintf(fname,"%s/positions_y.dat",spath); */
  /* output=fopen(fname,"w"); */
  /* for (gfn = 0; gfn < nfiles; gfn++) { */
  /*   fprintf(output,"%f %f",time[gfn],z[gfn]); */
  /*   for(ngal=0; ngal<ngalxs; ngal++) { */
  /*     fprintf(output," %f", pos_y[gfn][ngal]); */
  /*   } */
  /*   fprintf(output,"\n"); */
  /* } */
  /* fclose(output); */

  /* sprintf(fname,"%s/positions_z.dat",spath); */
  /* output=fopen(fname,"w"); */
  /* for (gfn = 0; gfn < nfiles; gfn++) { */
  /*   fprintf(output,"%f %f",time[gfn],z[gfn]); */
  /*   for(ngal=0; ngal<ngalxs; ngal++) { */
  /*     fprintf(output," %f", pos_z[gfn][ngal]); */
  /*   } */
  /*   fprintf(output,"\n"); */
  /* } */
  /* fclose(output); */

  /* sprintf(fname,"%s/velocity_x.dat",spath); */
  /* output=fopen(fname,"w"); */
  /* for (gfn = 0; gfn < nfiles; gfn++) { */
  /*   fprintf(output,"%f %f",time[gfn],z[gfn]); */
  /*   for(ngal=0; ngal<ngalxs; ngal++) { */
  /*     fprintf(output," %f", vel_x[gfn][ngal]); */
  /*   } */
  /*   fprintf(output,"\n"); */
  /* } */
  /* fclose(output); */

  /* sprintf(fname,"%s/velocity_y.dat",spath); */
  /* output=fopen(fname,"w"); */
  /* for (gfn = 0; gfn < nfiles; gfn++) { */
  /*   fprintf(output,"%f %f",time[gfn],z[gfn]); */
  /*   for(ngal=0; ngal<ngalxs; ngal++) { */
  /*     fprintf(output," %f", vel_y[gfn][ngal]); */
  /*   } */
  /*   fprintf(output,"\n"); */
  /* } */
  /* fclose(output); */

  /* sprintf(fname,"%s/velocity_z.dat",spath); */
  /* output=fopen(fname,"w"); */
  /* for (gfn = 0; gfn < nfiles; gfn++) { */
  /*   fprintf(output,"%f %f",time[gfn],z[gfn]); */
  /*   for(ngal=0; ngal<ngalxs; ngal++) { */
  /*     fprintf(output," %f", vel_z[gfn][ngal]); */
  /*   } */
  /*   fprintf(output,"\n"); */
  /* } */
  /* fclose(output); */

  /* writeout of psi files, which are readable within amira */
/*   /\* position for amira *\/ */
/*   for (gfn = 0; gfn < nfiles; gfn++) { */
/*     sprintf(fname,"%s/pos_amira_all_%04i.psi",spath,gfn+1); */
/*     output=fopen(fname,"w"); */
/*     fprintf(output,"# PSI Format 1.0\n"); */
/*     fprintf(output,"#\n"); */
/*     fprintf(output,"# column[0] = \"x\"\n"); */
/*     fprintf(output,"# column[1] = \"y\"\n"); */
/*     fprintf(output,"# column[2] = \"z\"\n"); */
/*     fprintf(output,"\n"); */
/*     fprintf(output,"%i 1 1\n",nGalaxies[gfn]); */
/*     fprintf(output,"1 1 1\n"); */
/*     fprintf(output,"1 1 1\n"); */
/*     fprintf(output,"1 1 1\n"); */
/*     fprintf(output,"\n"); */
/*     for(ngal=0; ngal<ngalxs; ngal++) { */
/*       if(pos_x[gfn][ngal]!=NAN */
/* ) */
/* 	fprintf(output,"%f %f %f\n", pos_x[gfn][ngal], pos_y[gfn][ngal], pos_z[gfn][ngal]); */
/*     } */
/*     fclose(output); */
/*   }
 */
  /* int counter; */
  /* for (gfn = 0; gfn < nfiles; gfn++) { */
  /*   sprintf(fname,"%s/pos_amira_5000_%04i.psi",spath,gfn+1); */
  /*   output=fopen(fname,"w"); */
  /*   fprintf(output,"# PSI Format 1.0\n"); */
  /*   fprintf(output,"#\n"); */
  /*   fprintf(output,"# column[0] = \"x\"\n"); */
  /*   fprintf(output,"# column[1] = \"y\"\n"); */
  /*   fprintf(output,"# column[2] = \"z\"\n"); */
  /*   fprintf(output,"\n"); */
  /*   counter = 0; */
  /*   for(ngal = 0; ngal < ngalxs; ngal++) { */
  /*     if(pos_x[gfn][ngal] <= 5000 && pos_y[gfn][ngal] <= 5000 && pos_z[gfn][ngal] <= 5000) */
  /* 	counter++; */
  /*   } */
  /*   fprintf(output,"%i 1 1\n",counter); */
  /*   fprintf(output,"1 1 1\n"); */
  /*   fprintf(output,"1 1 1\n"); */
  /*   fprintf(output,"1 1 1\n"); */
  /*   fprintf(output,"\n"); */
  /*   for(ngal=0; ngal<ngalxs; ngal++) { */
  /*     if(-5000 <= pos_x[gfn][ngal] <= 5000 && -5000 <= pos_y[gfn][ngal] <= 5000 && -5000 <= pos_z[gfn][ngal] <= 5000) */
  /* 	fprintf(output,"%f %f %f\n", pos_x[gfn][ngal], pos_y[gfn][ngal], pos_z[gfn][ngal]); */
  /*   } */
  /*   fclose(output); */
  /* } */


  for (gfn = 0; gfn < nfiles; gfn++)
    free(stf_histo[gfn]);
  free(stf_histo);
  for (gfn = 0; gfn < nfiles; gfn++)
    free(mcold_histo[gfn]);
  free(mcold_histo);
  for (gfn = 0; gfn < nfiles; gfn++)
    free(rdisk_histo[gfn]);
  free(rdisk_histo);
  for (gfn = 0; gfn < nfiles; gfn++)
    free(coldmetal_histo[gfn]);
  free(coldmetal_histo);
  for (gfn = 0; gfn < nfiles; gfn++)
    free(mhalo_histo[gfn]);
  free(mhalo_histo);
  for (gfn = 0; gfn < nfiles; gfn++)
    free(pos_x[gfn]);
  free(pos_x);
  for (gfn = 0; gfn < nfiles; gfn++)
    free(pos_y[gfn]);
  free(pos_y);
  for (gfn = 0; gfn < nfiles; gfn++)
    free(pos_z[gfn]);
  free(pos_z);
  for (gfn = 0; gfn < nfiles; gfn++)
    free(vel_x[gfn]);
  free(vel_x);
  for (gfn = 0; gfn < nfiles; gfn++)
    free(vel_y[gfn]);
  free(vel_y);
  for (gfn = 0; gfn < nfiles; gfn++)
    free(vel_z[gfn]);
  free(vel_z);

  free(cluster_stf);
  free(cluster_stf_avg);
  free(coldmetal_avg);
  free(nGalaxies);
  free(time);
  free(z);
  free(a);
  



  return 0;
}
