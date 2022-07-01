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

#ifdef FIX_CLASS

FixStyle(echemdid,FixEChemDID)

#else

#ifndef LMP_FIX_ECHEMDID_H
#define LMP_FIX_ECHEMDID_H

#include "fix.h"

namespace LAMMPS_NS {

class FixEChemDID : public Fix {
 public:
  FixEChemDID(class LAMMPS *, int, char **);
  ~FixEChemDID();
  int setmask();
  void init(), post_integrate();
  void init_list(int, class NeighList *);

  virtual int pack_forward_comm(int, int *, double *, int, int *);
  virtual void unpack_forward_comm(int, int, double *);
  void get_names(char *,double *&);
  double compute_vector(int);

 protected:
  int g1,g2,nlocal,nmax;
  double k,nelec,volt,rc,norm;
  double js[2],js_all[2];
  int eflag;
 private:
  void laplacian();
  class NeighList *list;
};

}

#endif
#endif
