/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_qeq_self.h"
#include <cmath>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024
#define EV2KCAL 23.02

/* ---------------------------------------------------------------------- */

PairQeqSelf::PairQeqSelf(LAMMPS *lmp) : Pair(lmp)
{
}

/* ---------------------------------------------------------------------- */

PairQeqSelf::~PairQeqSelf()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);

    memory->destroy(chi);
    memory->destroy(eta);
    memory->destroy(gamma);
    memory->destroy(zeta);
    memory->destroy(zcore);
  }
}

/* ---------------------------------------------------------------------- */

void PairQeqSelf::compute(int eflag, int vflag)
  {
  int i,ii,inum;
  int *ilist;
  double en_tmp;

  en_tmp = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;

  for(ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    if (eflag) { 
      en_tmp = EV2KCAL * (chi[type[i]] * q[i] + (eta[type[i]] / 2.) * q[i]*q[i]);
    }

    if (evflag) ev_tally(i,i,nlocal,newton_pair,0.0,en_tmp,0.0,0.0,0.0,0.0);
  }

//  if (vflag_fdotr) virial_fdotr_compute();

}

/* ---------------------------------------------------------------------- */

void PairQeqSelf::read_file(char *file)
{
  int i;
  int params_per_line = 6;
  char **words = new char*[params_per_line+1];

  int ntypes = atom->ntypes;
  int *setflag = new int[ntypes+1];
  for (i=0; i <= ntypes; ++i) setflag[i] = 0;

  memory->create(chi,ntypes+1,"qeq:chi");
  memory->create(eta,ntypes+1,"qeq:eta");
  memory->create(gamma,ntypes+1,"qeq:gamma");
  memory->create(zeta,ntypes+1,"qeq:zeta");
  memory->create(zcore,ntypes+1,"qeq:zcore");

  // open file on proc 0

  FILE *fp;
  if (comm->me == 0) {
    fp = force->open_potential(file);
    if (fp == NULL) {
      char str[128];
      snprintf(str,128,"Cannot open fix qeq parameter file %s",file);
      error->one(FLERR,str);
    }
  }

  int n,nwords,eof,nlo,nhi;
  char line[MAXLINE],*ptr;

  eof = 0;

  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = atom->count_words(line);
    if (nwords == 0) continue;

    if (nwords < 6)
      error->all(FLERR,"Invalid fix qeq parameter file");

    for (n=0, words[n] = strtok(line," \t\n\r\f");
         n < 6;
         words[++n] = strtok(NULL," \t\n\r\f"));

    force->bounds(FLERR,words[0],ntypes,nlo,nhi);
    for (n=nlo; n <=nhi; ++n) {
      chi[n]     = force->numeric(FLERR,words[1]);
      eta[n]     = force->numeric(FLERR,words[2]);
      gamma[n]   = force->numeric(FLERR,words[3]);
      zeta[n]    = force->numeric(FLERR,words[4]);
      zcore[n]   = force->numeric(FLERR,words[5]);
      setflag[n] = 1;
    }
  }

  for (n=1; n <= ntypes; ++n)
    if (setflag[n] == 0)
      error->all(FLERR,"Invalid fix qeq parameter file");

  delete [] words;
  delete [] setflag;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairQeqSelf::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);
  read_file(arg[1]);

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairQeqSelf::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double cut_one = cut_global;
  if (narg == 3) cut_one = force->numeric(FLERR,arg[2]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairQeqSelf::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");

  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 1;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairQeqSelf::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return cut[i][j];
}

