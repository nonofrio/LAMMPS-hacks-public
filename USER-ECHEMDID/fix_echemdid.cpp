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

/* ----------------------------------------------------------------------
   Contributing author: Nicolas Onofrio, The Hong Kong PolyU 
   (nicolas.onofrio@polyu.edu.hk)
------------------------------------------------------------------------- */

#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "fix_echemdid.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "force.h"
#include "kspace.h"
#include "group.h"
#include "pair.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "fix_qeq.h"
#include "modify.h"
#include "citeme.h"

using namespace LAMMPS_NS;
using namespace FixConst;

static const char cite_fix_echemdid[] =
  "fix echemdid command:\n\n"
  "@article{onofrio2015atomic,\n"
  " title={Atomic origin of ultrafast resistance switching in nanoscale electrometallization cells},\n"
  " author={Onofrio, Nicolas and Guzman, David and Strachan, Alejandro},\n"
  " journal={Nature materials},\n"
  " volume={14},\n"
  " number={4},\n"
  " pages={440},\n"
  " year={2015},\n"
  " publisher={Nature Publishing Group}\n"
  "}\n\n"
  "@article{onofrio2015voltage,\n"
  " title={Voltage equilibration for reactive atomistic simulations of electrochemical processes},\n"
  " author={Onofrio, Nicolas and Strachan, Alejandro},\n"
  " journal={The Journal of chemical physics},\n"
  " volume={143},\n"
  " number={5},\n"
  " pages={054109},\n"
  " year={2015},\n"
  " publisher={AIP Publishing}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

FixEChemDID::FixEChemDID(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_echemdid);
  if (narg < 12) error->all(FLERR,"Illegal fix EChemDID command");

// fix ID groupID EChemDID nevery k k_value cut cut_value norm norm_value 
// nelec nelec_value group1 group2 volt volt_value

  nevery = force->inumeric(FLERR,arg[3]);
  k = force->numeric(FLERR,arg[5]);
  rc = force->numeric(FLERR,arg[7]);
  norm = force->numeric(FLERR,arg[9]);
  nelec = force->inumeric(FLERR,arg[11]);
  g1 = group->bitmask[group->find(arg[12])];
  g2 = group->bitmask[group->find(arg[13])];
  volt = force->numeric(FLERR,arg[15]);

  if ((nevery <= 0) || (rc <= 0.0) || (norm <= 0.0) || (nelec <= 0))
    error->all(FLERR,"Illegal fix EChemDID command");

  if ((g1 == -1) || (g2 == -1))
    error->all(FLERR,"Could not find EChemDID boundary group ID");

  comm_forward = 3;

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix  = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  vector_flag = 1;
  size_vector = 2;
  eflag = 0;
  js[0] = js[1] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixEChemDID::~FixEChemDID()
{
}

/* ---------------------------------------------------------------------- */

void FixEChemDID::init()
{
  int flagqeq = 0;
  for (int i = 0; i < modify->nfix; i++) {
    if (strstr(modify->fix[i]->style,"qeq/shielded")) {
      flagqeq = 1;
    }
  }
  if (flagqeq != 1){
    error->all(FLERR,"Fix EChemDID must be called with qeq/shielded");
  }
}

/* ---------------------------------------------------------------------- */

void FixEChemDID::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

int FixEChemDID::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixEChemDID::post_integrate()
{
  int i, j;
  int nlocal = atom -> nlocal;
  double dt = update -> dt;
  int *mask = atom->mask;
  double *locpot, *lap;
  char tmp1[] = "locpot";
  char tmp2[] = "lap";
  get_names(tmp1,locpot);
  get_names(tmp2,lap);

  // Compute locpot 
  for (i = 0; i < nelec; i++) {
    laplacian();
    for (j = 0; j < nlocal; j++) {
      if (mask[j] & groupbit) {
        locpot[j] += k*dt*lap[j]/nelec; 
      }
    }
  }

  comm->forward_comm_fix(this);
}

/* ---------------------------------------------------------------------- */

void FixEChemDID::laplacian()
{
  int inum,jnum,i,j,ii,jj;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,w,n;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *locpot, *lap;
  char tmp1[] = "locpot";
  char tmp2[] = "lap";
  get_names(tmp1,locpot);
  get_names(tmp2,lap);

  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

// Get neighbors
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

// Set boundary conditions for locpot
  eflag = 0;
  js[0] = js[1] = 0.0;
  double dt = update -> dt;
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & g1) {
      js[0] += (0.5*volt-locpot[i])/dt;
      locpot[i] = 0.5*volt;
    }
    if (mask[i] & g2) {
      js[1] += (-0.5*volt-locpot[i])/dt;
      locpot[i] = -0.5*volt;
    }
  }

// Forward locpot to ghost
  comm->forward_comm_fix(this);

// Set laplacian and norm to 0
  for (i = 0; i < nlocal; i++) {
    lap[i] = 0.0;
  }
  n = 0.0;

// Loop over local i
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit) {

      jlist = firstneigh[i];
      jnum = numneigh[i]; 

      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];

// Loop over neighbor j of i
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;
//        printf("%d %d\n",i,j);
        if (mask[j] & groupbit) {

          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx * delx + dely * dely + delz * delz;

// Weighting function and normalization
          if (rsq < rc*rc) w = (1-(rsq)/(rc*rc))*(1-(rsq)/(rc*rc));
          else w = 0.0;

          n += w;
          if (j < nlocal) n += w;

// Compute laplacian 
          lap[i] += w * (locpot[j] - locpot[i]) / rsq;

        } 
      }
    } 
  }

// Compute normalization
//  n = n / (6.0 * nlocal);
//  printf("norm  %f",n);
  for (i = 0; i < nlocal; i++) {
    lap[i] = lap[i] / norm;
  }

}

/* ---------------------------------------------------------------------- */

void FixEChemDID::get_names(char *c,double *&ptr)
{
 int index,flag;
 index = atom->find_custom(c,flag);

 if(index!=-1) ptr = atom->dvector[index];
 else error->all(FLERR,"Fix EChemDID requires fix property/atom d_locpot d_lap command");
}

/* ---------------------------------------------------------------------- */

int FixEChemDID::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i,j,m;
  double *locpot;
  char tmp[] = "locpot";
  get_names(tmp,locpot);

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = locpot[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixEChemDID::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;
  double *locpot;
  char tmp[] = "locpot";
  get_names(tmp,locpot);

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) locpot[i] = buf[m++];
}

/* ---------------------------------------------------------------------- */

double FixEChemDID::compute_vector(int n)
{
  // only sum across procs one time

  if (eflag == 0) {
    MPI_Allreduce(js,js_all,2,MPI_DOUBLE,MPI_SUM,world);
    eflag = 1;
  }
  return js_all[n];
}

/* ------------------------------------------------------------------------- */
