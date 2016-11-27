/************************************************************************
 * functions.c native chemistry handling functions
 *
 * Copyright (c) 2004,2016 by Ernst-G. Schmid
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


#include <sys/time.h>
#include <sys/resource.h>
#include <postgres.h>
#include <fmgr.h>
#include <funcapi.h>
#include <access/htup_details.h>
#include <utils/varbit.h>
#include <mb/pg_wchar.h>
#include "functions.h"
#include "obwrapper.h"
#include "barsoi/barsoi.h"
#include "reaction/reaction.h"
#include "compression.h"

#ifdef PG_MODULE_MAGIC
PG_MODULE_MAGIC;
#endif


/*
* Returns the version of pgchem.
*/
PG_FUNCTION_INFO_V1 (pgchem_version);

Datum
pgchem_version (PG_FUNCTION_ARGS)
{
    PG_RETURN_CSTRING (PGCHEM_VERSION);
}

/*
* Returns the version of barsoi.
*/
PG_FUNCTION_INFO_V1 (pgchem_barsoi_version);

Datum
pgchem_barsoi_version (PG_FUNCTION_ARGS)
{
    char *barsoi_version_buffer;

    barsoi_version_buffer = (char *) palloc0 (256);

    xm_version (barsoi_version_buffer);

    PG_RETURN_CSTRING (barsoi_version_buffer);
}

/*
* Convert molecule to SMILES.
*/
PG_FUNCTION_INFO_V1 (pgchem_molecule_to_smiles);

Datum
pgchem_molecule_to_smiles (PG_FUNCTION_ARGS)
{
    char *tmpSmiles=NULL;
    text *retval;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
    bool omit_iso_and_chiral_markings = PG_GETARG_BOOL (1);
    int len;
    DECOMPRESSED_DATA *original_data=NULL;

    tmpSmiles =
        ob_smiles_to_smiles (SMIPTR (arg_molecule), omit_iso_and_chiral_markings ? 1 : 0);

    if(tmpSmiles == NULL)
    {
        goto smiles_fail;
    }
    else if (!strlen(tmpSmiles))
    {
        free(tmpSmiles);
smiles_fail:
        original_data=decompress_data(CIPTR(arg_molecule),arg_molecule->compressed_sizeo, arg_molecule->sizeo);
        elog (ERROR, "SMILES generation failed! Offender was :\n %s",original_data->decompressed_data);
    }

    len = strlen (tmpSmiles);

    retval = (text *) palloc0 (len + VARHDRSZ);
    //memset(retval,0x0,len + VARHDRSZ);

    SET_VARSIZE (retval,(len + VARHDRSZ));

    strncpy (VARDATA (retval), tmpSmiles, len);

    free (tmpSmiles);

    PG_RETURN_TEXT_P (retval);
}

/*
* Convert molecule to canonical SMILES.
*/
PG_FUNCTION_INFO_V1 (pgchem_molecule_to_canonical_smiles);

Datum
pgchem_molecule_to_canonical_smiles (PG_FUNCTION_ARGS)
{
    char *tmpSmiles=NULL;
    text *retval;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
    bool omit_iso_and_chiral_markings = PG_GETARG_BOOL (1);
    int len;
    DECOMPRESSED_DATA *original_data=NULL;

    tmpSmiles =
        ob_smiles_to_canonical_smiles (SMIPTR (arg_molecule), omit_iso_and_chiral_markings ? 1 : 0);

    if(tmpSmiles == NULL)
    {
        goto smiles_fail;
    }
    else if (!strlen(tmpSmiles))
    {
        free(tmpSmiles);
smiles_fail:
        original_data = decompress_data(CIPTR(arg_molecule),arg_molecule->compressed_sizeo, arg_molecule->sizeo);
        elog (ERROR, "Canonical SMILES generation failed! Offender was :\n %s",original_data->decompressed_data);
    }

    len = strlen (tmpSmiles);

    retval = (text *) palloc0 (len + VARHDRSZ);
    //memset(retval,0x0,len + VARHDRSZ);

    SET_VARSIZE (retval,(len + VARHDRSZ));

    strncpy (VARDATA (retval), tmpSmiles, len);

    free (tmpSmiles);

    PG_RETURN_TEXT_P (retval);
}

/*
* Convert molecule to InChI.
*/
PG_FUNCTION_INFO_V1 (pgchem_molecule_to_inchi);

Datum
pgchem_molecule_to_inchi (PG_FUNCTION_ARGS)
{
    char *tmpInChI = NULL;
    text *retval;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
    int len;
    DECOMPRESSED_DATA *original_data=NULL;

    tmpInChI = ob_smiles_to_inchi (SMIPTR (arg_molecule));

    if(tmpInChI == NULL)
    {
        goto inchi_fail;
    }
    else if (!strlen(tmpInChI) || strstr (tmpInChI, "InChI=") == NULL)
    {
        free(tmpInChI);
inchi_fail:
        original_data = decompress_data(CIPTR(arg_molecule),arg_molecule->compressed_sizeo, arg_molecule->sizeo);
        elog (ERROR, "InChI generation failed! Offender was :\n %s",original_data->decompressed_data);
    }

    len = strlen (tmpInChI);

    retval = (text *) palloc0 (len + VARHDRSZ);
    //memset(retval,0x0,len + VARHDRSZ);

    SET_VARSIZE (retval,(len + VARHDRSZ));

    strncpy (VARDATA (retval), tmpInChI, len);

    free (tmpInChI);

    PG_RETURN_TEXT_P (retval);
}

/*
* Convert molecule to V3000.
*/
PG_FUNCTION_INFO_V1 (pgchem_molecule_to_V3000);

Datum
pgchem_molecule_to_V3000 (PG_FUNCTION_ARGS)
{
    char *tmpV3000=NULL;
    text *retval;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
    int len;
    DECOMPRESSED_DATA *original_data=NULL;

    if(arg_molecule->original_format==FORMAT_V3000)
    {
        original_data = decompress_data(CIPTR(arg_molecule),arg_molecule->compressed_sizeo, arg_molecule->sizeo);
        tmpV3000 = strdup(original_data->decompressed_data);
    }
    else
    {
        tmpV3000 = ob_smiles_to_V3000 (SMIPTR(arg_molecule));
    }

    if(tmpV3000 == NULL)
    {
        goto v3000_fail;
    }
    else if (!strlen(tmpV3000) || strstr (tmpV3000, "V3000") == NULL)
    {
        free(tmpV3000);
v3000_fail:
        elog (ERROR, "V3000 molfile generation failed! Offender was :\n %s",original_data->decompressed_data);
    }

    len = strlen (tmpV3000);

    retval = (text *) palloc0 (len + VARHDRSZ);
    //memset(retval,0x0,len+VARHDRSZ);

    SET_VARSIZE (retval,(len + VARHDRSZ));

    strncpy (VARDATA (retval), tmpV3000, len);

    free (tmpV3000);

    PG_RETURN_TEXT_P (retval);
}

PG_FUNCTION_INFO_V1 (pgchem_molecule_to_V2000);

Datum
pgchem_molecule_to_V2000 (PG_FUNCTION_ARGS)
{
    char *tmpMolfile=NULL;
    text *retval;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
    int len;
    DECOMPRESSED_DATA *original_data=NULL;

    if(arg_molecule->original_format==FORMAT_V2000)
    {
        original_data = decompress_data(CIPTR(arg_molecule),arg_molecule->compressed_sizeo, arg_molecule->sizeo);
        tmpMolfile = strdup(original_data->decompressed_data);
    }
    else
    {
        tmpMolfile = ob_smiles_to_V2000 (SMIPTR(arg_molecule));
    }

    if(tmpMolfile == NULL)
    {
        goto molfile_fail;
    }
    else if (!strlen(tmpMolfile) || strstr (tmpMolfile, "V2000") == NULL)
    {
        free(tmpMolfile);
molfile_fail:
        elog (ERROR, "V2000 molfile generation failed! Offender was :\n %s",original_data->decompressed_data);
    }

    len = strlen (tmpMolfile);

    retval = (text *) palloc0 (len + VARHDRSZ);
    //memset(retval,0x0,len+VARHDRSZ);

    SET_VARSIZE (retval,(len + VARHDRSZ));

    strncpy (VARDATA (retval), tmpMolfile, len);

    free (tmpMolfile);

    PG_RETURN_TEXT_P (retval);
}

/*
* Calculate molweight of molecule.
*/
PG_FUNCTION_INFO_V1 (pgchem_molweight);

Datum
pgchem_molweight (PG_FUNCTION_ARGS)
{
    float8 retval = 0.0;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

    retval = ob_molweight (SMIPTR (arg_molecule));

    PG_RETURN_FLOAT8 (retval);
}

/*
* Calculate Hill formula of molecule.
*/
PG_FUNCTION_INFO_V1 (pgchem_hillformula);

Datum
pgchem_hillformula (PG_FUNCTION_ARGS)
{
    char *tmpFormula;
    text *retval;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
    int len;

    tmpFormula = ob_hillformula (SMIPTR (arg_molecule));

    len = strlen (tmpFormula);

    retval = (text *) palloc0 (len + VARHDRSZ);
    //memset(retval,0x0,len + VARHDRSZ);

    SET_VARSIZE (retval,(len + VARHDRSZ));

    strncpy (VARDATA (retval), tmpFormula, len);

    free (tmpFormula);

    PG_RETURN_TEXT_P (retval);
}

/*
* Calculate exact mass of molecule.
*/
PG_FUNCTION_INFO_V1 (pgchem_exactmass);

Datum
pgchem_exactmass (PG_FUNCTION_ARGS)
{
    float8 retval = 0.0;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

    retval = ob_exactmass (SMIPTR (arg_molecule));

    PG_RETURN_FLOAT8 (retval);
}

/*
* Calculate total charge of molecule.
*/
PG_FUNCTION_INFO_V1 (pgchem_total_charge);

Datum
pgchem_total_charge (PG_FUNCTION_ARGS)
{
    int32 retval = 0;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

    retval = ob_total_charge (SMIPTR (arg_molecule));

    PG_RETURN_INT32 (retval);
}

/*
* Calculate number of atoms of molecule.
*/
PG_FUNCTION_INFO_V1 (pgchem_num_atoms);

Datum
pgchem_num_atoms (PG_FUNCTION_ARGS)
{
    uint32 retval = 0;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

    retval = ob_num_atoms (SMIPTR (arg_molecule));

    PG_RETURN_UINT32 (retval);
}

/*
* Calculate number of heavy atoms of molecule.
*/
PG_FUNCTION_INFO_V1 (pgchem_num_heavy_atoms);

Datum
pgchem_num_heavy_atoms (PG_FUNCTION_ARGS)
{
    uint32 retval = 0;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

    retval = ob_num_heavy_atoms (SMIPTR (arg_molecule));

    PG_RETURN_UINT32 (retval);
}

/*
* Calculate number of bonds of molecule.
*/
PG_FUNCTION_INFO_V1 (pgchem_num_bonds);

Datum
pgchem_num_bonds (PG_FUNCTION_ARGS)
{
    uint32 retval = 0;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

    retval = ob_num_bonds (SMIPTR (arg_molecule));

    PG_RETURN_UINT32 (retval);
}

/*
* Calculate number of rotatable bonds of molecule.
*/
PG_FUNCTION_INFO_V1 (pgchem_num_rotatable_bonds);

Datum
pgchem_num_rotatable_bonds (PG_FUNCTION_ARGS)
{
    uint32 retval = 0;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

    retval =
        ob_num_rotatable_bonds (SMIPTR (arg_molecule));

    PG_RETURN_UINT32 (retval);
}

/*
* Determine chirality of molecule.
*/
PG_FUNCTION_INFO_V1 (pgchem_is_chiral);

Datum
pgchem_is_chiral (PG_FUNCTION_ARGS)
{
    bool retval = false;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

    retval =
        (ob_is_chiral (SMIPTR(arg_molecule)) ==1) ? true : false;

    PG_RETURN_BOOL (retval);
}

/*
* Molecule has 2D coordinates.
*/
PG_FUNCTION_INFO_V1 (pgchem_2D);

Datum
pgchem_2D (PG_FUNCTION_ARGS)
{
    bool retval = false;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
    DECOMPRESSED_DATA *molfile=NULL;

    if(arg_molecule->original_format == FORMAT_V2000 || arg_molecule->original_format==FORMAT_V3000)
    {
        molfile = decompress_data(CIPTR(arg_molecule),arg_molecule->compressed_sizeo, arg_molecule->sizeo);

        retval =
            (ob_2D (molfile->decompressed_data) == 1) ? true : false;
    }

    PG_RETURN_BOOL (retval);
}

/*
* Molecule has 3D coordinates.
*/
PG_FUNCTION_INFO_V1 (pgchem_3D);

Datum
pgchem_3D (PG_FUNCTION_ARGS)
{
    bool retval = false;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
    DECOMPRESSED_DATA *molfile=NULL;

    if(arg_molecule->original_format == FORMAT_V2000 || arg_molecule->original_format==FORMAT_V3000)
    {
        molfile = decompress_data(CIPTR(arg_molecule),arg_molecule->compressed_sizeo, arg_molecule->sizeo);

        retval =
            (ob_3D (molfile->decompressed_data) == 1) ? true : false;
    }

    PG_RETURN_BOOL (retval);
}

/*
* Calculate molstatistics fingerprint of molecule.
* Long version: n_atoms:6;n_bonds:6...
*/
PG_FUNCTION_INFO_V1 (pgchem_ms_fingerprint_long_a);

Datum
pgchem_ms_fingerprint_long_a (PG_FUNCTION_ARGS)
{
    MOLECULE *arg_molecule;
    text *retval;
    int len;

    bool strict_chg = PG_GETARG_BOOL (1);
    bool strict_iso = PG_GETARG_BOOL (2);
    bool strict_rad = PG_GETARG_BOOL (3);
    //char molfile[MAX_MOL_SIZE + 1];

    const int fpbuffer_len = 1024;
    char fpbuffer[fpbuffer_len];
    char *molfile=NULL;

    arg_molecule = PG_GETARG_MOLECULE_P (0);

    molfile = ob_smiles_to_V2000(SMIPTR(arg_molecule));

    mm_set_chg_check (strict_chg ? FEATURE_ON : FEATURE_OFF);
    mm_set_iso_check (strict_iso ? FEATURE_ON : FEATURE_OFF);
    mm_set_rad_check (strict_rad ? FEATURE_ON : FEATURE_OFF);
    xm_set_ring_perception_algorithm (RPA_SAR);
    cm_set_mol (molfile, FEATURE_ON);
    cm_molstat (fpbuffer);

    free(molfile);

    len = strlen (fpbuffer);

    retval = (text *) palloc0 (len + VARHDRSZ);

    SET_VARSIZE (retval,(len + VARHDRSZ));

    /* copy the return value to a variable of psql type text */
    memcpy (VARDATA (retval), fpbuffer, len);

    /* return as psql type by reference */
    PG_RETURN_TEXT_P (retval);
}

/*
* Calculate molstatistics fingerprint of molecule.
* Short version: 6,6,1,0,0,0,0,6,6,0,0...
*/
PG_FUNCTION_INFO_V1 (pgchem_ms_fingerprint_short_a);

Datum
pgchem_ms_fingerprint_short_a (PG_FUNCTION_ARGS)
{
    MOLECULE *arg_molecule;
    text *retval;
    int len;

    //char molfile[MAX_MOL_SIZE + 1];

    char fpbuffer[1024];

    char *molfile=NULL;

    arg_molecule = PG_GETARG_MOLECULE_P (0);

    molfile = ob_smiles_to_V2000(SMIPTR(arg_molecule));

    memset(fpbuffer,0x0,1024*sizeof(char));

    xm_set_ring_perception_algorithm (RPA_SAR);
    cm_set_mol (molfile, FEATURE_ON);
    cm_molstat_X (fpbuffer);

    free(molfile);

    len = strlen (fpbuffer);

    retval = (text *) palloc0 (len + VARHDRSZ);

    SET_VARSIZE (retval,(len + VARHDRSZ));

    /* copy the return value to a variable of psql type text */
    memcpy (VARDATA (retval), fpbuffer, len);

    /* return as psql type by reference */
    PG_RETURN_TEXT_P (retval);
}
/*
* Calculate functional group codes of molecule.
* 0000A000;...
*/
PG_FUNCTION_INFO_V1 (pgchem_fgroup_codes_a);

Datum
pgchem_fgroup_codes_a (PG_FUNCTION_ARGS)
{
    MOLECULE *arg_molecule;
    text *retval;
    int len;

    //char molfile[MAX_MOL_SIZE + 1];

    const int fpbuffer_len = 1024;
    char fpbuffer[fpbuffer_len];

    char *molfile=NULL;

    arg_molecule = PG_GETARG_MOLECULE_P (0);

    molfile = ob_smiles_to_V2000(SMIPTR(arg_molecule));

    xm_set_ring_perception_algorithm (RPA_SAR);
    cm_set_mol (molfile, FEATURE_ON);
    cm_fg_codes (fpbuffer);

    free(molfile);

    len = strlen (fpbuffer);

    retval = (text *) palloc0 (len + VARHDRSZ);

    SET_VARSIZE (retval,(len + VARHDRSZ));

    /* copy the return value to a variable of psql type text */
    memcpy (VARDATA (retval), fpbuffer, len);

    /* return as psql type by reference */
    PG_RETURN_TEXT_P (retval);
}

/*
* Match molecule against SMARTS pattern, as requested by Mr. Ertl from Novartis. :-)
*/
PG_FUNCTION_INFO_V1 (pgchem_smartsfilter);

Datum
pgchem_smartsfilter (PG_FUNCTION_ARGS)
{
    bool retval = false;
    text *arg_smarts = PG_GETARG_TEXT_P (0);
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (1);
    int smarts_string_len = VARSIZE (arg_smarts) - VARHDRSZ;
    int match = 0;
    char *tmpSMARTS = (char *) palloc0 (smarts_string_len + 1);
    //tmpSMARTS[0] = '\0';

    strncat (tmpSMARTS, VARDATA (arg_smarts), smarts_string_len);

    match = ob_SSS_SMARTS_native (tmpSMARTS, SMIPTR(arg_molecule));

    if(match<0) elog (ERROR, "Invalid SMARTS pattern: %s",tmpSMARTS);

    retval = (match != 0);

    PG_RETURN_BOOL (retval);
}

PG_FUNCTION_INFO_V1 (pgchem_smartsfilter_count);

Datum
pgchem_smartsfilter_count (PG_FUNCTION_ARGS)
{
    //TODO: rewrite for Indigo
    int retval = 0;
    text *arg_smarts = PG_GETARG_TEXT_P (0);
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (1);
    int smarts_string_len = VARSIZE (arg_smarts) - VARHDRSZ;
    char *tmpSMARTS = (char *) palloc0 (smarts_string_len + 1);
    //tmpSMARTS[0] = '\0';

    strncat (tmpSMARTS, VARDATA (arg_smarts), smarts_string_len);

    retval = ob_SSS_SMARTS_native_count (tmpSMARTS, SMIPTR(arg_molecule));

    if(retval<0) elog (ERROR, "Invalid SMARTS pattern: %s",tmpSMARTS);

    PG_RETURN_INT32 (retval);
}

/*
* Predict polar surface area (PSA) of molecule.
*/
PG_FUNCTION_INFO_V1 (pgchem_TPSA);

Datum
pgchem_TPSA (PG_FUNCTION_ARGS)
{
    float8 retval = 0.0;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

    retval = ob_TPSA (SMIPTR (arg_molecule));

    PG_RETURN_FLOAT8 (retval);
}

/*
* Predict molar refractivity (MR) of molecule.
*/
PG_FUNCTION_INFO_V1 (pgchem_MR);

Datum
pgchem_MR (PG_FUNCTION_ARGS)
{
    float8 retval = 0.0;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

    retval = ob_MR (SMIPTR (arg_molecule));

    PG_RETURN_FLOAT8 (retval);
}

/*
* Predict logP of molecule.
*/
PG_FUNCTION_INFO_V1 (pgchem_logP);

Datum
pgchem_logP (PG_FUNCTION_ARGS)
{
    float8 retval = 0.0;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

    retval = ob_logP (SMIPTR (arg_molecule));

    PG_RETURN_FLOAT8 (retval);
}

PG_FUNCTION_INFO_V1 (pgchem_molecule_to_inchikey);

Datum
pgchem_molecule_to_inchikey (PG_FUNCTION_ARGS)
{
    text *retval;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

    retval = (text *) palloc0 (INCHIKEYSZ + VARHDRSZ);

    //memset(retval,0x0,INCHIKEYSZ + VARHDRSZ);

    memcpy (VARDATA (retval), arg_molecule->inchikey, INCHIKEYSZ);

    SET_VARSIZE (retval,(INCHIKEYSZ + VARHDRSZ));

    PG_RETURN_TEXT_P (retval);
}

/*
* Get binary fingerprint as bit varying
*/
PG_FUNCTION_INFO_V1 (pgchem_fp_out);

Datum
pgchem_fp_out (PG_FUNCTION_ARGS)
{
    VarBit *retval;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
    int len;

    len = (OB_FPSIZE2+OB_FPSIZE3)*sizeof(uint32);

    retval = (VarBit *) palloc0 (len + VARBITHDRSZ);

    memcpy(VARBITS(retval),arg_molecule->fp.bytes,len);

    VARBITLEN(retval) = len*8; //8 bit chars

    SET_VARSIZE(retval,(len+VARBITHDRSZ));

    PG_RETURN_VARBIT_P (retval);
}

/*
* Get MACCS binary fingerprint as bit varying
*/
PG_FUNCTION_INFO_V1 (pgchem_fp_MACCS);

Datum
pgchem_fp_MACCS (PG_FUNCTION_ARGS)
{
    VarBit *retval;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
    int len = OB_FPSIZE_MACCS*sizeof(uint32);
    uint32 *tmp_maccs;

    retval = (VarBit *) palloc0 (len + VARBITHDRSZ);

    tmp_maccs = (uint32*) palloc0(len);

    ob_fp_MACCS(SMIPTR(arg_molecule), tmp_maccs);

    memcpy(VARBITS(retval),tmp_maccs,len);

    VARBITLEN(retval) = len*8; //8 bit chars

    SET_VARSIZE(retval,(len+VARBITHDRSZ));

    PG_RETURN_VARBIT_P (retval);
}

/*
* Get reaction binary fingerprint as bit varying
*/
PG_FUNCTION_INFO_V1 (pgchem_r_fp_out);

Datum
pgchem_r_fp_out (PG_FUNCTION_ARGS)
{
    VarBit *retval;
    REACTION *arg_reaction = PG_GETARG_REACTION_P (0);
    int len = 2*OB_FPSIZE2*sizeof(uint32);

    retval = (VarBit *) palloc0 (len + VARBITHDRSZ);

    memcpy(VARBITS(retval),arg_reaction->fp,len);

    VARBITLEN(retval) = len*8; //8 bit chars

    SET_VARSIZE(retval,(len+VARBITHDRSZ));

    PG_RETURN_VARBIT_P (retval);
}

PG_FUNCTION_INFO_V1 (pgchem_describe_fp2);

Datum
pgchem_describe_fp2 (PG_FUNCTION_ARGS)
{
    char *tmpDESC=NULL;
    text *retval;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
    int len;
    DECOMPRESSED_DATA *original_data=NULL;

    tmpDESC = ob_describe_fp2 (SMIPTR (arg_molecule));

    if(tmpDESC == NULL)
    {
        goto desc_fail;
    }
    else if (!strlen(tmpDESC) || strstr (tmpDESC, ">") == NULL)
    {
        free(tmpDESC);
desc_fail:
        original_data = decompress_data(CIPTR(arg_molecule),arg_molecule->compressed_sizeo, arg_molecule->sizeo);
        elog (ERROR, "Description generation failed! Offender was :\n %s",original_data->decompressed_data);
    }

    len = strlen (tmpDESC);

    retval = (text *) palloc0 (len + VARHDRSZ);
    //memset(retval,0x0,len+VARHDRSZ);

    SET_VARSIZE (retval,(len + VARHDRSZ));

    strncpy (VARDATA (retval), tmpDESC, len);

    free (tmpDESC);

    PG_RETURN_TEXT_P (retval);
}


PG_FUNCTION_INFO_V1 (pgchem_num_H_donors);

Datum
pgchem_num_H_donors (PG_FUNCTION_ARGS)
{
    uint32 retval = 0;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

    retval = ob_num_H_donors (SMIPTR (arg_molecule));

    PG_RETURN_UINT32 (retval);
}

PG_FUNCTION_INFO_V1 (pgchem_num_H_acceptors);

Datum
pgchem_num_H_acceptors (PG_FUNCTION_ARGS)
{
    uint32 retval = 0;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

    retval = ob_num_H_acceptors (SMIPTR (arg_molecule));

    PG_RETURN_UINT32 (retval);
}


PG_FUNCTION_INFO_V1 (pgchem_mutate_fp);

Datum
pgchem_mutate_fp (PG_FUNCTION_ARGS)
{
    //TODO: broken
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
    uint32 *offset = arg_molecule->fp.dwords+OB_OFFSET;

    ob_fp3(SMIPTR(arg_molecule), offset);

    PG_RETURN_MOLECULE_P (arg_molecule);
}

PG_FUNCTION_INFO_V1 (pgchem_blank_fp);

Datum
pgchem_blank_fp (PG_FUNCTION_ARGS)
{
    //TODO: broken
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
    uint32 *offset = arg_molecule->fp.dwords+OB_OFFSET;

    memset(offset,0x0,OB_FPSIZE3*(sizeof(uint32)));

    PG_RETURN_MOLECULE_P (arg_molecule);
}


PG_FUNCTION_INFO_V1(pgchem_nbits_set);
Datum
pgchem_nbits_set(PG_FUNCTION_ARGS)
{
    /* how many bits are set in a bitstring? */

    VarBit     *a = PG_GETARG_VARBIT_P(0);
    unsigned char *ap = (unsigned char*) VARBITS(a);

    PG_RETURN_INT32(ob_popcount(ap,VARBITBYTES(a)));
}

PG_FUNCTION_INFO_V1 (pgchem_disconnected);

Datum
pgchem_disconnected (PG_FUNCTION_ARGS)
{
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

    PG_RETURN_BOOL (arg_molecule->disconnected != 0);
}

PG_FUNCTION_INFO_V1 (pgchem_r_num_reactants);

Datum
pgchem_r_num_reactants(PG_FUNCTION_ARGS)
{
    REACTION *arg_reaction = PG_GETARG_REACTION_P (0);

    PG_RETURN_INT32(arg_reaction->num_reactants);
}

PG_FUNCTION_INFO_V1 (pgchem_r_num_products);

Datum
pgchem_r_num_products(PG_FUNCTION_ARGS)
{
    REACTION *arg_reaction = PG_GETARG_REACTION_P (0);

    PG_RETURN_INT32(arg_reaction->num_products);
}

PG_FUNCTION_INFO_V1 (pgchem_r_molecule_at);

Datum
pgchem_r_molecule_at(PG_FUNCTION_ARGS)
{
    REACTION *arg_reaction = PG_GETARG_REACTION_P (0);
    int32 arg_position = PG_GETARG_INT32 (1);
    MOLECULE *retval;
    char *offset = MOLARRAYPTR(arg_reaction);
    int i,len;

    if(arg_position < 1 || arg_position > arg_reaction->num_products+arg_reaction->num_reactants) elog (ERROR, "Molecule index out of bounds: %d",arg_position);

    for(i=1; i<arg_position; i++)
    {
        len = VARSIZE((MOLECULE*)offset)*sizeof(char);
        offset+=len;
    }

    len = VARSIZE((MOLECULE*)offset)*sizeof(char);

    retval = (MOLECULE*) palloc0(len);

    //memset(retval,0x0,len);

    memcpy(retval,(MOLECULE*)offset,len);

    PG_RETURN_MOLECULE_P(retval);
}

PG_FUNCTION_INFO_V1 (pgchem_r_reaction_to_smiles);

Datum
pgchem_r_reaction_to_smiles (PG_FUNCTION_ARGS)
{
    REACTION *reaction = PG_GETARG_REACTION_P (0);
    MOLECULE *tmpMol;
    char* offset = MOLARRAYPTR(reaction);
    text *result;
    char* tmpBuf;
    int i, size=0;

    for(i=0; i<reaction->num_products+reaction->num_reactants; i++)
    {
        size+=((MOLECULE*) offset)->sizesmi*sizeof(char);
        offset+=VARSIZE((MOLECULE*) offset);
    }

    offset = MOLARRAYPTR(reaction);

    tmpBuf = (char *) palloc0 (size+sizeof(char));

    //memset(tmpBuf,0x0,size+sizeof(char));

    for(i=0; i<reaction->num_reactants; i++)
    {
        tmpMol = (MOLECULE*) offset;
        if(strstr(SMIPTR(tmpMol),"\r\n") != NULL)
        {
            strncat(tmpBuf,SMIPTR(tmpMol),tmpMol->sizesmi-3);
        }
        else if(strstr(SMIPTR(tmpMol),"\n") != NULL)
        {
            strncat(tmpBuf,SMIPTR(tmpMol),tmpMol->sizesmi-2);
        }
        if(i<reaction->num_reactants-1) strncat(tmpBuf,".",sizeof(char));
        offset+=VARSIZE(tmpMol)*sizeof(char);
    }

    strncat(tmpBuf,">>",2*sizeof(char));

    for(i=0; i<(reaction->num_products); i++)
    {
        tmpMol = (MOLECULE*) offset;
        if(strstr(SMIPTR(tmpMol),"\r\n") != NULL)
        {
            strncat(tmpBuf,SMIPTR(tmpMol),tmpMol->sizesmi-3);
        }
        else if(strstr(SMIPTR(tmpMol),"\n") != NULL)
        {
            strncat(tmpBuf,SMIPTR(tmpMol),tmpMol->sizesmi-2);
        }
        if(i<reaction->num_products-1) strncat(tmpBuf,".",sizeof(char));
        offset+=VARSIZE(tmpMol)*sizeof(char);
    }

    result = (text *) palloc0(strlen(tmpBuf)+VARHDRSZ);

    //memset(result,0x0,strlen(tmpBuf)+VARHDRSZ);

    memcpy(VARDATA(result),tmpBuf,strlen(tmpBuf));

    SET_VARSIZE(result,strlen(tmpBuf)+VARHDRSZ);

    pfree(tmpBuf);

    PG_RETURN_TEXT_P (result);
}

PG_FUNCTION_INFO_V1 (pgchem_tversky);

Datum pgchem_tversky (PG_FUNCTION_ARGS)
{
    //TODO: migrate to tversky_16
    MOLECULE *mol1_prototype = PG_GETARG_MOLECULE_P (0);
    MOLECULE *mol2_variant = PG_GETARG_MOLECULE_P (1);
    float8 alpha_prototype = PG_GETARG_FLOAT8 (2);
    float8 beta_variant = PG_GETARG_FLOAT8 (3);

    PG_RETURN_FLOAT8 (ob_tversky
                      (mol1_prototype->fp.dwords,  mol2_variant->fp.dwords, OB_FPSIZE2, alpha_prototype, beta_variant));
}

PG_FUNCTION_INFO_V1 (pgchem_spectrophore);

Datum
pgchem_spectrophore (PG_FUNCTION_ARGS)
{
    text *retval;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
    char* spectrophore;
    int len;

    spectrophore = ob_spectrophore(SMIPTR(arg_molecule));

    len = strlen(spectrophore);

    retval = (text *) palloc0 (len + VARHDRSZ);

    //memset(retval,0x0,len + VARHDRSZ);

    memcpy (VARDATA (retval), spectrophore, len);

    SET_VARSIZE (retval,(len + VARHDRSZ));

    PG_RETURN_TEXT_P (retval);
}






