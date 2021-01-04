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
   Contributing authors: Ray Shan (Sandia), 
                         Nicolas Onofrio & Udoka Nwankwo (PolyU)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include "fix_qeq_shielded_lr.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "force.h"
#include "group.h"
#include "kspace.h"
#include "respa.h"
#include "memory.h"
#include "error.h"
#include "math_vector.h"
#include <vector>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixQEqShieldedLR::FixQEqShieldedLR(LAMMPS *lmp, int narg, char **arg) :
  FixQEq(lmp, narg, arg) {}

/* ---------------------------------------------------------------------- */

void FixQEqShieldedLR::init()
{
  if (!atom->q_flag)
    error->all(FLERR,"Fix qeq/shielded requires atom attribute q");

  ngroup = group->count(igroup);
  if (ngroup == 0) error->all(FLERR,"Fix qeq/shielded group has no atoms");

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix  = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  int ntypes = atom->ntypes;
  memory->create(shld,ntypes+1,ntypes+1,"qeq:shielding");

  init_shielding();

  int i;
  for (i = 1; i <= ntypes; i++) {
    if (gamma[i] == 0.0)
      error->all(FLERR,"Invalid param file for fix qeq/shielded");
  }

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

}

/* ---------------------------------------------------------------------- */

void FixQEqShieldedLR::init_shielding()
{
  int i,j;
  double d7, swa2, swa3, swb2, swb3;

  int ntypes = atom->ntypes;
  for( i = 1; i <= ntypes; ++i )
    for( j = 1; j <= ntypes; ++j )
      shld[i][j] = pow( gamma[i] * gamma[j], -1.5 );

  if (fabs(swa) > 0.01 && comm->me == 0)
    error->warning(FLERR,"Fix qeq has non-zero lower Taper radius cutoff");
  if (swb < 0)
    error->all(FLERR, "Fix qeq has negative upper Taper radius cutoff");
  else if (swb < 5 && comm->me == 0)
    error->warning(FLERR,"Fix qeq has very low Taper radius cutoff");

  d7 = pow( swb - swa, 7 );
  swa2 = swa*swa;
  swa3 = swa2*swa;
  swb2 = swb*swb;
  swb3 = swb2*swb;

  Tap[7] =  20.0 / d7;
  Tap[6] = -70.0 * (swa + swb) / d7;
  Tap[5] =  84.0 * (swa2 + 3.0*swa*swb + swb2) / d7;
  Tap[4] = -35.0 * (swa3 + 9.0*swa2*swb + 9.0*swa*swb2 + swb3 ) / d7;
  Tap[3] = 140.0 * (swa3*swb + 3.0*swa2*swb2 + swa*swb3 ) / d7;
  Tap[2] =-210.0 * (swa3*swb2 + swa2*swb3) / d7;
  Tap[1] = 140.0 * swa3 * swb3 / d7;
  Tap[0] = (-35.0*swa3*swb2*swb2 + 21.0*swa2*swb3*swb2 -
            7.0*swa*swb3*swb3 + swb3*swb3*swb ) / d7;
}

/* ---------------------------------------------------------------------- */

void FixQEqShieldedLR::pre_force(int vflag)
{
  if (update->ntimestep % nevery) return;

  nlocal = atom->nlocal;

  if( atom->nmax > nmax ) reallocate_storage();

  if( nlocal > n_cap*DANGER_ZONE || m_fill > m_cap*DANGER_ZONE )
    reallocate_matrix();

  reallocate_matrix();

  init_matvec();
  matvecs = CG(b_s, s);         // CG on s - parallel
  matvecs += CG(b_t, t);        // CG on t - parallel
  calculate_Q();
  
  if (force->kspace) force->kspace->qsum_qsq();
}

/* ---------------------------------------------------------------------- */

void FixQEqShieldedLR::init_matvec()
{
  compute_H();

  int inum, ii, i, nat;
  int *ilist;

  int flag_echemdid; // EChemDID
  double *locpot; // EChemDID
  char tmp1[] = "locpot"; // EChemDID
  flag_echemdid = get_names(tmp1,locpot); // EChemDID

  inum = list->inum;
  ilist = list->ilist;
  
  for( ii = 0; ii < inum; ++ii ) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit) {
      Hdia_inv[i] = 1. / eta[ atom->type[i] ];
      if (flag_echemdid == 1) b_s[i] = -chi[atom->type[i]] - locpot[i]; // EChemDID
      else b_s[i] = -chi[atom->type[i]]; //-chi[atom->type[i]]+chi[atom->type[0]];
      b_t[i]      = -1.0;
      t[i] = t_hist[i][2] + 3 * ( t_hist[i][0] - t_hist[i][1] );
      s[i] = 4*(s_hist[i][0]+s_hist[i][2])-(6*s_hist[i][1]+s_hist[i][3]);
      //cout << "b_ =" <<b_s[i]<< endl;
    }
  }
 
  pack_flag = 2;
  comm->forward_comm_fix(this); //Dist_vector( s );
  pack_flag = 3;
  comm->forward_comm_fix(this); //Dist_vector( t );
}

/* ---------------------------------------------------------------------- */

void FixQEqShieldedLR::compute_H()
{
  int inum, jnum, *ilist, *jlist, *numneigh, **firstneigh;
  int i, j, ii, jj;
  double **x;
  double dx, dy, dz, r_sqr, r;
 
  int *type = atom->type;
  x = atom->x;
  int *mask = atom->mask;
  
  /*int nat;
  nat = atom->natoms; 
  std::vector<double> dummyRow(nat,0);
  std::vector<std::vector<double> > AA(nat,dummyRow);
  

  for (i = 0; i < nat; i++){AA[0][i] = 1.;}    // First row of A is all ones	
	
  for (i = 1; i < nat; i++) { for (j = 0; j < nat; j++) {AA[i][j] = calculate_H(i, j) - calculate_H(0, j);}}  */
  

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
//fill in the H matrix
  m_fill = 0;
  r_sqr = 0;
  for( ii = 0; ii < inum; ii++ ) {
    i = ilist[ii];
    if (mask[i] & groupbit) {
      jlist = firstneigh[i];
      jnum = numneigh[i];
      H.firstnbr[i] = m_fill;

      for( jj = 0; jj < jnum; jj++ ) {
        j = jlist[jj];
        j &= NEIGHMASK;

        dx = x[j][0] - x[i][0];
        dy = x[j][1] - x[i][1];
        dz = x[j][2] - x[i][2];
        r_sqr = dx*dx + dy*dy + dz*dz;

        if (r_sqr <= cutoff_sq) {
          H.jlist[m_fill] = j;
          r = sqrt(r_sqr);
          H.val[m_fill] = 0.5 * calculate_H(i, j); //AA[i][j]; //SR_shielded( r, shld[type[i]][type[j]] );
          m_fill++; 
        }
      }
      H.numnbrs[i] = m_fill - H.firstnbr[i];
    }
  }
  
  if (m_fill >= H.m) {
    char str[128];
    sprintf(str,"H matrix size has been exceeded: m_fill=%d H.m=%d\n",m_fill, H.m );
    error->warning(FLERR,str);
    error->all(FLERR,"Fix qeq/shielded has insufficient QEq matrix size");
  }
}

/* ---------------------------------------------------------------------- */

double FixQEqShieldedLR::calculate_H( int i, int j )
{
  double *h = domain->h;
  double lx,ly,lz,xy,xz,yz;
  lx = h[0];
  ly = h[1];
  lz = h[2];
  yz = h[3];
  xz = h[4];
  xy = h[5];

// Compute reciprocal lattice
// Convert lammps -> a, b, c, alpha, beta, gamma
  double a,b,c,angle_alpha,angle_beta,angle_gamma;
  a = lx;
  b = sqrt(ly*ly+xy*xy);
  c = sqrt(lz*lz+xz*xz+yz*yz);
  angle_alpha = acos((xy*xz+ly*yz)/(b*c));
  angle_beta = acos(xz/c);
  angle_gamma = acos(xy/b);

// Convert a, b, c, alpha, beta, gamma -> a1, a2, a3
  double ax,ay,az,bx,by,bz,cx,cy,cz;
  double a1[3], a2[3], a3[3], b1[3], b2[3], b3[3];
  ax = a;
  ay = 0.0;
  az = 0.0;
  bx = b*cos(angle_gamma);
  by = b*sin(angle_gamma);
  bz = 0.0;
  cx = c*cos(angle_beta);
  cy = (b*c*cos(angle_alpha)-bx*cx)/by; //(b*c-bx*cx)/by; // Not sure here
  cz = sqrt(c*c-cx*cx-cy*cy);

  a1[0] = ax;
  a1[1] = ay;
  a1[2] = az;
  a2[0] = bx;
  a2[1] = by;
  a2[2] = bz;
  a3[0] = cx;
  a3[1] = cy;
  a3[2] = cz;

  double vec[3];
  double vf;
  cross(a2,a3,vec);
  vf = 2.0*M_PI/dot(a1,vec);

// b1
  cross(a2,a3,vec);
  for (int ii = 0; ii < 3; ii++) {b1[ii] = 2*M_PI*vec[ii]/dot(a1,vec);};

// b2
  cross(a3,a1,vec);
  for (int ii = 0; ii < 3; ii++) {b2[ii] = 2*M_PI*vec[ii]/dot(a2,vec);};

// b3
  cross(a1,a2,vec);
  for (int ii = 0; ii < 3; ii++) {b3[ii] = 2*M_PI*vec[ii]/dot(a3,vec);};

  int nat, u, v, w;
  nat = atom->natoms;
  int *type = atom->type;
  double **x = atom->x;
  double gamma_shld;
  int Nmax = 6;
  double zeta = 0.27;
  double lmda = 1.2;
  double dx,dy,dz;
  double dhx,dhy,dhz;
  double dr, dh;
  double Jij;
  gamma_shld= shld[type[i]][type[j]];
  
  if (i == j){
     double Astar = 0.0;  
     double Bstar = 0.0;  
     for (u = -Nmax; u <= Nmax; u++){
         for (v = -Nmax; v <= Nmax; v++){
             for (w = -Nmax; w <= Nmax; w++){
                 if (!(u == 0 && v == 0 && w == 0)){ 
                    dx = u*a1[0] + v*a2[0] + w*a3[0];
                    dy = u*a1[1] + v*a2[1] + w*a3[1];
                    dz = u*a1[2] + v*a2[2] + w*a3[2];
                    dhx = u*b1[0] + v*b2[0] + w*b3[0];
                    dhy = u*b1[1] + v*b2[1] + w*b3[1];
                    dhz = u*b1[2] + v*b2[2] + w*b3[2];
                    dr = sqrt(dx*dx + dy*dy + dz*dz); 
                    dh = sqrt(dhx*dhx + dhy*dhy + dhz*dhz);
                    Astar += erfc(dr*sqrt(zeta))/pow((dr*dr*dr+gamma_shld),0.33333333333333);
                    Bstar += (2.0*vf/(dh*dh))*exp((-1.*dh*0.5*dh*0.5)/zeta);  
                 }
             }
         }
     } 
     Jij = eta[type[i]]+0.5*lmda*EV_TO_KCAL_PER_MOL*(Astar+Bstar-2.0*sqrt(zeta/M_PI));
  } else {
          double Alpha = 0.0;
          double Beta = 0.0;
          for (u = -Nmax; u <= Nmax; u++){
              for (v = -Nmax; v <= Nmax; v++){
                  for (w = -Nmax; w <= Nmax; w++){
                      dx = x[j][0] - x[i][0] + u*a1[0] + v*a2[0] + w*a3[0];
                      dy = x[j][1] - x[i][1] + u*a1[1] + v*a2[1] + w*a3[1];
                      dz = x[j][2] - x[i][2] + u*a1[2] + v*a2[2] + w*a3[2];
                      dr = sqrt(dx*dx+dy*dy+dz*dz); 
                      Alpha += erfc(dr*sqrt(zeta))/pow((dr*dr*dr+gamma_shld),0.33333333333333);
                      if (!(u == 0 && v == 0 && w == 0)){ 
                         dhx = u*b1[0] + v*b2[0] + w*b3[0];
                         dhy = u*b1[1] + v*b2[1] + w*b3[1];
                         dhz = u*b1[2] + v*b2[2] + w*b3[2];
                         dx = x[j][0] - x[i][0];
                         dy = x[j][1] - x[i][1];
	                 dz = x[j][2] - x[i][2];
                         dh = sqrt(dhx*dhx + dhy*dhy + dhz*dhz);
                         Beta += (2.0*vf/(dh*dh))*exp(-1.0*(dh*0.5*dh*0.5)/zeta)*cos(dx*dhx + dy*dhy + dz*dhz); // Here careful! 
                      }
                  }
              }
          }
    Jij = 0.5*lmda*EV_TO_KCAL_PER_MOL*(Alpha+Beta);
  }
  return Jij;
}

/* ----------------------------------------------------------------------*/ 

double FixQEqShieldedLR::SR_shielded( double r, double gamma)
{
  double Taper, denom;

  Taper = Tap[7] * r + Tap[6];
  Taper = Taper * r + Tap[5];
  Taper = Taper * r + Tap[4];
  Taper = Taper * r + Tap[3];
  Taper = Taper * r + Tap[2];
  Taper = Taper * r + Tap[1];
  Taper = Taper * r + Tap[0];

  denom = r * r * r + gamma;
  denom = pow(denom,0.3333333333333);

  return Taper * EV_TO_KCAL_PER_MOL / denom;
}

/* ----------------------------------------------------------------------
   return x dot y
------------------------------------------------------------------------- */

double FixQEqShieldedLR::dot(double *x, double *y)
{
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

/* ----------------------------------------------------------------------
   z = x cross y
------------------------------------------------------------------------- */

void FixQEqShieldedLR::cross(double *x, double *y, double *z)
{
  z[0] = x[1]*y[2] - x[2]*y[1];
  z[1] = x[2]*y[0] - x[0]*y[2];
  z[2] = x[0]*y[1] - x[1]*y[0];
}
