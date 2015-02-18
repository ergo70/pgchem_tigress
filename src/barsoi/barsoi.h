/************************************************************************
 * barsoi.h barsoi library functions
 *  
 * Copyright (c) 2004 by Ernst-G. Schmid
 *   
 * This file is part of the xchem::tigress project.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the lesser GNU General Public License as published by
 * the Free Software Foundation version 2.1 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * lesser GNU General Public License for more details.
 ************************************************************************/

#ifdef __cplusplus
extern "C"
{
#endif
  void mm_set_mol (const char *);
  void cm_set_mol (const char *,int);
  void mm_set_current_mol_as_query (void);
  int mm_match();
  void xm_version (char *);
  int mm_get_rings ();
  int mm_get_atom_ring (int);
  void xm_set_ring_perception_algorithm (int);
  void xm_set_strict_typing (int);
  void mm_set_r_s_check (int);
  void mm_set_e_z_check (int);
  void mm_set_chg_check (int);
  void mm_set_rad_check (int);
  void mm_set_iso_check (int);
  void mm_set_exact_match (int);
  void cm_fg_codes (char *);
  void cm_molstat_X (char *);
  void cm_molstat_X_arr (unsigned short *);
  void cm_molstat (char *);
  void cm_molstat_sql_exact (char *);
  void cm_molstat_sql_substruct (char *);
  void cm_tweak_molfile (char *);
#ifdef __cplusplus
}
#endif

#define RPA_DEFAULT 0
#define RPA_SAR 1
#define RPA_SSR 2

#define FEATURE_OFF 0
#define FEATURE_ON 1

#ifndef DLLEXPORT
#ifdef WIN32
#define DLLEXPORT __declspec (dllexport)
#else
#define DLLEXPORT
#endif
#endif
