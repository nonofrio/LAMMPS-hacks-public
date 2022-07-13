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

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_morse_shielded_tap.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairMorseShieldedTap::PairMorseShieldedTap(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairMorseShieldedTap::~PairMorseShieldedTap()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(d0);
    memory->destroy(alpha);
    memory->destroy(gamma_w);
    memory->destroy(p_vdW1);
    memory->destroy(r_vdW);
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void PairMorseShieldedTap::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,dr,dexp,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);

        double swa = 0.0;
        double swb = sqrt(cutsq[itype][jtype]);
        double d7 = pow( swb, 7.0 );
        double swa2 = pow( swa, 2.0 );
        double swa3 = pow( swa, 3.0 );
        double swb2 = pow( swb, 2.0 );
        double swb3 = pow( swb, 3.0 );

        double Tap7 =  20.0 / d7;
        double Tap6 = -70.0 * (swa + swb) / d7;
        double Tap5 =  84.0 * (swa2 + 3.0*swa*swb + swb2) / d7;
        double Tap4 = -35.0 * (swa3 + 9.0*swa2*swb + 9.0*swa*swb2 + swb3 ) / d7;
        double Tap3 = 140.0 * (swa3*swb + 3.0*swa2*swb2 + swa*swb3 ) / d7;
        double Tap2 =-210.0 * (swa3*swb2 + swa2*swb3) / d7;
        double Tap1 = 140.0 * swa3 * swb3 / d7;
        double Tap0 = (-35.0*swa3*swb2*swb2 + 21.0*swa2*swb3*swb2 -
                           7.0*swa*swb3*swb3 + swb3*swb3*swb ) / d7;

        double Tap = Tap7 * r + Tap6;
        Tap = Tap * r + Tap5;
        Tap = Tap * r + Tap4;
        Tap = Tap * r + Tap3;
        Tap = Tap * r + Tap2;
        Tap = Tap * r + Tap1;
        Tap = Tap * r + Tap0;

        double dTap = 7*Tap7 * r + 6*Tap6;
        dTap = dTap * r + 5*Tap5;
        dTap = dTap * r + 4*Tap4;
        dTap = dTap * r + 3*Tap3;
        dTap = dTap * r + 2*Tap2;
        dTap += Tap1/r;

        double powr_vdW1 = pow(r, p_vdW1[itype][jtype]);
        double powgi_vdW1 = pow( 1.0 / gamma_w[itype][jtype], p_vdW1[itype][jtype]);
        double fn13 = pow( powr_vdW1 + powgi_vdW1, (1.0 / p_vdW1[itype][jtype]) );
        double exp1 = exp( alpha[itype][jtype] * (1.0 - fn13 / r_vdW[itype][jtype]) );
        double exp2 = exp( 0.5 * alpha[itype][jtype] * (1.0 - fn13 / r_vdW[itype][jtype]) );
        double dfn13 = pow( powr_vdW1 + powgi_vdW1, (1.0 / p_vdW1[itype][jtype]) - 1.0) * pow(r, p_vdW1[itype][jtype] - 2.0);
        double f_vdW = dTap * d0[itype][jtype] * (exp1 - 2.0 * exp2) -
            Tap * d0[itype][jtype] * (alpha[itype][jtype] / r_vdW[itype][jtype]) * (exp1 - exp2) * dfn13;

//        printf("r %f powr_vdW1 %f powgi_vdW1 %f fn13 %f exp1 %f exp2 %f dfn13 %f f_vdW %f\n",r,powr_vdW1,powgi_vdW1,fn13,exp1,exp2,dfn13,f_vdW);

        fpair = factor_lj * f_vdW * (-1.0); 

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          evdwl = Tap * d0[itype][jtype] * (exp1 - 2.0 * exp2);
//            offset[itype][jtype];
          evdwl *= factor_lj;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairMorseShieldedTap::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(d0,n+1,n+1,"pair:d0");
  memory->create(alpha,n+1,n+1,"pair:alpha");
  memory->create(gamma_w,n+1,n+1,"pair:gamma_w");
  memory->create(p_vdW1,n+1,n+1,"pair:p_vdW1");
  memory->create(r_vdW,n+1,n+1,"pair:r_vdW");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMorseShieldedTap::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);

  // reset cutoffs that have been explicitly set

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

void PairMorseShieldedTap::coeff(int narg, char **arg)
{
  if (narg < 7 || narg > 8)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double d0_one = force->numeric(FLERR,arg[2]);
  double alpha_one = force->numeric(FLERR,arg[3]);
  double gamma_one = force->numeric(FLERR,arg[4]);
  double p_one = force->numeric(FLERR,arg[5]);
  double r_one = force->numeric(FLERR,arg[6]);

  double cut_one = cut_global;
  if (narg == 8) cut_one = force->numeric(FLERR,arg[7]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      d0[i][j] = d0_one;
      alpha[i][j] = alpha_one;
      gamma_w[i][j] = gamma_one;
      p_vdW1[i][j] = p_one;
      r_vdW[i][j] = 2.0*r_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairMorseShieldedTap::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  if (offset_flag) {
    double powr_vdW1 = pow(cut[i][j], p_vdW1[i][j]);
    double powgi_vdW1 = pow( 1.0 / gamma_w[i][j], p_vdW1[i][j]);
    double fn13 = pow( powr_vdW1 + powgi_vdW1, (1.0 / p_vdW1[i][j]) );
    double exp1 = exp( alpha[i][j] * (1.0 - fn13 / r_vdW[i][j]) );
    double exp2 = exp( 0.5 * alpha[i][j] * (1.0 - fn13 / r_vdW[i][j]) );
    offset[i][j] = d0[i][j] * (exp1 - 2.0 * exp2);
  } else offset[i][j] = 0.0;

  d0[j][i] = d0[i][j];
  alpha[j][i] = alpha[i][j];
  gamma_w[j][i] = gamma_w[i][j];
  p_vdW1[j][i] = p_vdW1[i][j];
  r_vdW[j][i] = r_vdW[i][j];
  offset[j][i] = offset[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairMorseShieldedTap::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&d0[i][j],sizeof(double),1,fp);
        fwrite(&alpha[i][j],sizeof(double),1,fp);
        fwrite(&gamma_w[i][j],sizeof(double),1,fp);
        fwrite(&p_vdW1[i][j],sizeof(double),1,fp);
        fwrite(&r_vdW[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairMorseShieldedTap::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&d0[i][j],sizeof(double),1,fp);
          fread(&alpha[i][j],sizeof(double),1,fp);
          fread(&gamma_w[i][j],sizeof(double),1,fp);
          fread(&p_vdW1[i][j],sizeof(double),1,fp);
          fread(&r_vdW[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&d0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&alpha[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&gamma_w[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&p_vdW1[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&r_vdW[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairMorseShieldedTap::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairMorseShieldedTap::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairMorseShieldedTap::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g %g\n",i,d0[i][i],alpha[i][i],gamma_w[i][i],p_vdW1[i][i],r_vdW[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairMorseShieldedTap::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g %g\n",
              i,j,d0[i][j],alpha[i][j],gamma_w[i][j],p_vdW1[i][j],r_vdW[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairMorseShieldedTap::single(int i, int j, int itype, int jtype, double rsq,
                         double factor_coul, double factor_lj,
                         double &fforce)
{
  double r,dr,dexp,phi;
  r = sqrt(rsq);

  double swa = 0.0; 
  double swb = sqrt(cutsq[itype][jtype]); 
  double d7 = pow( swb, 7.0 );
  double swa2 = pow( swa, 2.0 );
  double swa3 = pow( swa, 3.0 );
  double swb2 = pow( swb, 2.0 );
  double swb3 = pow( swb, 3.0 );

  double Tap7 =  20.0 / d7;
  double Tap6 = -70.0 * (swa + swb) / d7;
  double Tap5 =  84.0 * (swa2 + 3.0*swa*swb + swb2) / d7;
  double Tap4 = -35.0 * (swa3 + 9.0*swa2*swb + 9.0*swa*swb2 + swb3 ) / d7;
  double Tap3 = 140.0 * (swa3*swb + 3.0*swa2*swb2 + swa*swb3 ) / d7;
  double Tap2 =-210.0 * (swa3*swb2 + swa2*swb3) / d7;
  double Tap1 = 140.0 * swa3 * swb3 / d7;
  double Tap0 = (-35.0*swa3*swb2*swb2 + 21.0*swa2*swb3*swb2 -
                     7.0*swa*swb3*swb3 + swb3*swb3*swb ) / d7;

  double Tap = Tap7 * r + Tap6;
  Tap = Tap * r + Tap5;
  Tap = Tap * r + Tap4;
  Tap = Tap * r + Tap3;
  Tap = Tap * r + Tap2;
  Tap = Tap * r + Tap1;
  Tap = Tap * r + Tap0;

  double dTap = 7*Tap7 * r + 6*Tap6;
  dTap = dTap * r + 5*Tap5;
  dTap = dTap * r + 4*Tap4;
  dTap = dTap * r + 3*Tap3;
  dTap = dTap * r + 2*Tap2;
  dTap += Tap1/r;

  double powr_vdW1 = pow(r, p_vdW1[itype][jtype]);
  double powgi_vdW1 = pow( 1.0 / gamma_w[itype][jtype], p_vdW1[itype][jtype]);
  double fn13 = pow( powr_vdW1 + powgi_vdW1, (1.0 / p_vdW1[itype][jtype]) );
  double exp1 = exp( alpha[itype][jtype] * (1.0 - fn13 / r_vdW[itype][jtype]) );
  double exp2 = exp( 0.5 * alpha[itype][jtype] * (1.0 - fn13 / r_vdW[itype][jtype]) );
  double dfn13 = pow( powr_vdW1 + powgi_vdW1, (1.0 / p_vdW1[itype][jtype]) - 1.0) * pow(r, p_vdW1[itype][jtype] - 2.0);
  phi = Tap * d0[itype][jtype] * (exp1 - 2.0 * exp2); 
  fforce = -1.0 * (dTap * phi - Tap * d0[itype][jtype] * (alpha[itype][jtype] / r_vdW[itype][jtype]) * (exp1 - exp2) * dfn13);
  return factor_lj*phi;
}

/* ---------------------------------------------------------------------- */

void *PairMorseShieldedTap::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"d0") == 0) return (void *) d0;
  if (strcmp(str,"gamma_w") == 0) return (void *) gamma_w;
  if (strcmp(str,"alpha") == 0) return (void *) alpha;
  if (strcmp(str,"p_vdW1") == 0) return (void *) p_vdW1;
  if (strcmp(str,"r_vdW") == 0) return (void *) r_vdW;
  return NULL;
}
