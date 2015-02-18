/************************************************************************
 * functions.h native chemistry handling functions
 *
 * Copyright (c) 2004,2012 by Ernst-G. Schmid
 *
 * This file is part of the xchem::tigress project.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 2 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * lesser GNU General Public License for more details.
 ************************************************************************/

#include "pgchem_config.h"

Datum pgchem_version (PG_FUNCTION_ARGS);
Datum pgchem_barsoi_version (PG_FUNCTION_ARGS);
Datum pgchem_indigo_version (PG_FUNCTION_ARGS);

//Datum pgchem_is_nostruct (PG_FUNCTION_ARGS);
Datum pgchem_fg_fingerprint (PG_FUNCTION_ARGS);

Datum pgchem_ms_fingerprint_long (PG_FUNCTION_ARGS);
Datum pgchem_ms_fingerprint_short (PG_FUNCTION_ARGS);
Datum pgchem_fgroup_codes (PG_FUNCTION_ARGS);

Datum pgchem_molecule_to_smiles (PG_FUNCTION_ARGS);
Datum pgchem_molecule_to_canonical_smiles (PG_FUNCTION_ARGS);
Datum pgchem_molecule_to_inchi (PG_FUNCTION_ARGS);
Datum pgchem_molecule_to_inchikey (PG_FUNCTION_ARGS);

//Datum pgchem_molecule_to_svg (PG_FUNCTION_ARGS);
//Datum pgchem_molecule_to_png (PG_FUNCTION_ARGS);

Datum pgchem_total_charge (PG_FUNCTION_ARGS);
Datum pgchem_num_atoms (PG_FUNCTION_ARGS);
Datum pgchem_num_heavy_atoms (PG_FUNCTION_ARGS);
Datum pgchem_num_bonds (PG_FUNCTION_ARGS);
Datum pgchem_num_rotatable_bonds (PG_FUNCTION_ARGS);
Datum pgchem_is_chiral (PG_FUNCTION_ARGS);
Datum pgchem_2D (PG_FUNCTION_ARGS);
Datum pgchem_3D (PG_FUNCTION_ARGS);
Datum pgchem_molecule_to_V3000 (PG_FUNCTION_ARGS);
Datum pgchem_molecule_to_V2000 (PG_FUNCTION_ARGS);
//Datum pgchem_VF2 (PG_FUNCTION_ARGS);
Datum pgchem_smartsfilter (PG_FUNCTION_ARGS);
Datum pgchem_smartsfilter_count (PG_FUNCTION_ARGS);


Datum pgchem_molweight (PG_FUNCTION_ARGS);


Datum pgchem_hillformula (PG_FUNCTION_ARGS);


Datum pgchem_exactmass (PG_FUNCTION_ARGS);

/***************************************************************
 *  Accelerated functions from barsoi
 ***************************************************************/

Datum pgchem_ms_fingerprint_long_a (PG_FUNCTION_ARGS);
Datum pgchem_ms_fingerprint_short_a (PG_FUNCTION_ARGS);
Datum pgchem_fgroup_codes_a (PG_FUNCTION_ARGS);
//Datum pgchem_tweak_molecule_a (PG_FUNCTION_ARGS);

Datum pgchem_TPSA (PG_FUNCTION_ARGS);
Datum pgchem_MR (PG_FUNCTION_ARGS);
Datum pgchem_logP (PG_FUNCTION_ARGS);

Datum pgchem_num_H_donors (PG_FUNCTION_ARGS);
Datum pgchem_num_H_acceptors (PG_FUNCTION_ARGS);

Datum pgchem_fp_out (PG_FUNCTION_ARGS);
Datum pgchem_fp_MACCS (PG_FUNCTION_ARGS);
Datum pgchem_describe_fp2 (PG_FUNCTION_ARGS);

Datum pgchem_nbits_set(PG_FUNCTION_ARGS);

Datum pgchem_disconnected(PG_FUNCTION_ARGS);

#ifndef BUILD_WITH_INDIGO
Datum pgchem_mutate_fp (PG_FUNCTION_ARGS);
Datum pgchem_blank_fp (PG_FUNCTION_ARGS);
#endif

Datum pgchem_r_fp_out (PG_FUNCTION_ARGS);
Datum pgchem_r_num_reactants (PG_FUNCTION_ARGS);
Datum pgchem_r_num_products (PG_FUNCTION_ARGS);
Datum pgchem_r_molecule_at (PG_FUNCTION_ARGS);
Datum pgchem_r_reaction_to_smiles (PG_FUNCTION_ARGS);
Datum pgchem_tversky (PG_FUNCTION_ARGS);
Datum pgchem_spectrophore (PG_FUNCTION_ARGS);
Datum pgchem_isotope_pattern (PG_FUNCTION_ARGS);

//Datum reversestring(PG_FUNCTION_ARGS);

//Datum pgchem_msp (PG_FUNCTION_ARGS);




