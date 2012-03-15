/* This file is part of OOPS.
 *
 * OOPS is (c) 2003-2009 Jacek Gondzio and Andreas Grothey, 
 *                       University of Edinburgh
 *
 * OOPS is distributed in a restricted form in the hope that it will be a useful
 * example of what can be done with SML, however it is NOT released under a free
 * software license.
 *
 * You may only redistribute this version of OOPS with a version of SML. You
 * may not link OOPS with code which is not part of SML licensed under the
 * LGPL v3.
 *
 * You may NOT modify, disassemble, or otherwise reverse engineer OOPS.
 *
 * OOPS is distributed WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 */

#ifndef GLOBALOPT_H
#define GLOBALOPT_H

typedef struct GlobalOpt
{
  int prt;
  int iter_limit;
  double conv_tol;
  int get_diff_dir;
  int n_diff_dir;
  double diff_dir_blocklevel;
  int force_feas;
  int det_num_prob;
  int cure_num_prob;
  double init_regxy;
  double ro_reg;
  double default_pdregth;
  int always_hoc; 
  int weighted_hoc;
  int LA_reorder;
  int LA_use_sparse_schur;
  int WS_do;         
  double WS_gap_adv_cen;
  double WS_target_mu;
  int WS_ret_adv_cen_it;
  int WS_beg_cen_it;
  int WS_target;
  int WS_adjust_to_bnds;
  int WS_adjust_z;
  double WS_balance_xz;
  int WS_unblk;
  int WS_n_unblk_it;
  double WS_unblk_lev;
  int WS_ForwardSens_TenWorst_Analyse;
  int WS_ForwardSens_TenWorst_TakeStep;

}
GlobalOpt;

extern GlobalOpt *glopt;

GlobalOpt *
ReadGlobalOpt(void);

#endif
