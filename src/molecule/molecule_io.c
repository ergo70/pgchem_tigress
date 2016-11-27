/************************************************************************
 * molecule_io.c molecule input/output support functions
 *
 * Copyright (c) 2007,2016 by Ernst-G. Schmid
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
#include "postgres.h"
#include "libpq/md5.h"
#include "fmgr.h"
#include "libpq/pqformat.h"	/* needed for send/recv functions */
#include "molecule.h"
#include "obwrapper.h"
#include "barsoi/barsoi.h"
#include <assert.h>
#include "compression.h"

//static char bzhash[MOLHASHSZ];
//static bool inited = false;

//Datum pgchem_molecule_to_new_molecule (PG_FUNCTION_ARGS);
Datum molecule_in (PG_FUNCTION_ARGS);
Datum molecule_in_text (PG_FUNCTION_ARGS);
Datum molecule_in_varchar (PG_FUNCTION_ARGS);
Datum molecule_in_bytea (PG_FUNCTION_ARGS);
Datum molecule_out (PG_FUNCTION_ARGS);
Datum molecule_recv (PG_FUNCTION_ARGS);
Datum molecule_send (PG_FUNCTION_ARGS);
//Datum pgchem_inchi_to_molecule (PG_FUNCTION_ARGS);
//Datum pgchem_V3000_to_molecule (PG_FUNCTION_ARGS);
Datum pgchem_strip_salts (PG_FUNCTION_ARGS);
Datum pgchem_add_hydrogens (PG_FUNCTION_ARGS);
Datum pgchem_remove_hydrogens (PG_FUNCTION_ARGS);
Datum pgchem_original_format (PG_FUNCTION_ARGS);


MOLECULE *
new_molecule (char *smiles, char *original_data)
{
    unsigned int sizeo;
    unsigned int sizesmi;
    size_t totalsize;
    MOLECULE *result;
    char *inchikey = NULL;
    COMPRESSED_DATA *compressed_data;

    sizeo = strlen (original_data)+1;

    sizesmi = strlen (smiles)+1;

    compressed_data = compress_data(original_data, sizeo);

    totalsize = CALCDATASZ (compressed_data->compressed_size, sizesmi);

    result = (MOLECULE *) palloc0 (totalsize);

    //memset (result, 0x0, totalsize);

    if (strchr (smiles, '.') != NULL)
        result->disconnected = true;

    result->sizeo = sizeo;
    result->compressed_sizeo=compressed_data->compressed_size;
    result->sizesmi = sizesmi;

    strncpy (SMIPTR(result), smiles, sizesmi);

    //strncpy (MFPTR(result), molfile, sizemf);
    memcpy(CIPTR(result),compressed_data->compressed_data, compressed_data->compressed_size);

    inchikey = ob_smiles_to_inchikey (smiles);

    //printf("%s\n",inchi);
    //printf("%s\n",smiles);
    //printf("%s\n",molfile);

    if(inchikey == NULL)
    {
        goto inchikey_fail;
    }
    else if (strlen(inchikey) != INCHIKEYSZ)
    {
        elog (ERROR, "%s INCHIKEYSIZE %d != %ld\n",inchikey, INCHIKEYSZ, strlen(inchikey));
        free(inchikey);
inchikey_fail:
        elog (ERROR, "Molecule generation failed! Offender was :\n %s",original_data);
    }

    //pg_md5_hash (inchi, strlen (inchi) + 1, result->molhash);
    //calculateDigestFromBuffer(inchi, strlen (inchi), result->molhash);

    memcpy(result->inchikey, inchikey, INCHIKEYSZ);

    //calculateDigestFromBuffer(smiles, sizesmi, result->molhash);

    //result->molhash[16]='\0';

    free (inchikey);

    ob_fp(smiles, result->fp.dwords);

    //result->fp.popcount = ob_popcount(result->fp.dwords,OB_FPSIZE);

    /*if(ancillarydata != NULL)
    {
        //printf("%d %d %d %d\n",ancsize,offset[0],offset[1],offset[2]);
        free(ancillarydata);
    }*/

    //memset(fp3,0x0,FPSIZE3*sizeof(unsigned int));
    //ob_fp2 (molfile, result->fp);

    //offset = result->fp+OFFSET;

    //ob_fp3 (molfile, offset);

    //result->popcount = ob_popcount((uint8 *)result->fp,FPSIZE*sizeof(uint32));

    //printf("popcount: %i\n",result->popcount);

    //merge_fps(result->fp2, fp3);

    SET_VARSIZE (result,totalsize);

    return result;
}


static MOLECULE *make_molecule(const char *raw_input, int size)
{
    MOLECULE *result = NULL;
    char *input = NULL;
    char *smiles = NULL;
    char *endptr = NULL;
    //bool freemolfile = false;
    bool freesmiles = false;
    unsigned int new_len = 0;
    unsigned int input_format = FORMAT_SMILES;

    if(strstr (raw_input, "M  END") != NULL)
    {
        input = palloc0 (size+sizeof(char));
        memcpy (input, raw_input, size);
        endptr = strstr (input, "M  END") + strlen("M  END")*sizeof(char);
        *endptr = 0x0;
        new_len = strlen(input);
        pfree (input);
        input = palloc0 (new_len + 1);
        strncpy(input,raw_input,new_len);
        input[new_len] = 0x0;
    }
    else
    {
        input = palloc0 (size+1);
        memcpy (input, raw_input, size);
        input[size] = 0x0;
    }

    if (strstr (input, "V2000") != NULL || strstr (input, "M  END") != NULL)
    {
        if (strstr (input, "V3000") != NULL)
        {
            input_format = FORMAT_V3000;
        }
        else
        {
            input_format = FORMAT_V2000;
        }
    }
    else if  (strstr (input, "InChI=") != NULL)
    {
        input_format = FORMAT_INCHI;
    }

    switch (input_format)
    {
    case FORMAT_V3000:
        //molfile = ob_V3000_to_molfile (input);
        smiles = ob_molfile_to_canonical_smiles (input,0);
        //freemolfile = true;
        freesmiles = true;
        break;
    case FORMAT_V2000:
        //molfile = input;
        smiles = ob_molfile_to_canonical_smiles (input,0);
        freesmiles = true;
        break;
    case FORMAT_INCHI:
        //molfile = ob_inchi_to_molfile (input);
        smiles = ob_inchi_to_canonical_smiles (input,0);
        //freemolfile = true;
        freesmiles = true;
        break;
    case FORMAT_SMILES:
        //molfile = ob_smiles_to_molfile (input);
        smiles = input;
        //freemolfile = true;
        break;
    default:
        if (input != NULL)
            pfree(input);
        elog (ERROR, "Input format is not V2000/V3000 molfile or InChI or SMILES: %s", input);
    }

    if(smiles == NULL || !strlen(smiles))
    {
        if (smiles != NULL && freesmiles)
            free (smiles);
        /*if (molfile != NULL && freemolfile)
            free (molfile);*/
        elog (ERROR, "Internal SMILES generation failed! Offender was :\n %s",input);
    }

    result = new_molecule (smiles, input);

    result->original_format = input_format;

    if (smiles != NULL && freesmiles)
        free (smiles);
    /*if (molfile != NULL && freemolfile)
        free (molfile);*/

    //if (input != NULL)
    //    pfree(input);

    return result;
}

/*
* Convert a molecule in text form (V2000, SMILES, InChI) into a molecule
*/
PG_FUNCTION_INFO_V1 (molecule_in);

Datum
molecule_in (PG_FUNCTION_ARGS)
{
    char *input = PG_GETARG_CSTRING (0);
    int size = strlen(input);

    //printf("molecule_in\n");

    PG_RETURN_MOLECULE_P (make_molecule(input,size));
}

PG_FUNCTION_INFO_V1 (molecule_in_text);

Datum
molecule_in_text (PG_FUNCTION_ARGS)
{
    text *x = PG_GETARG_TEXT_P (0);
    char *input = VARDATA(x);
    int size = VARSIZE(x)-VARHDRSZ;

    PG_RETURN_MOLECULE_P (make_molecule(input,size));
}

PG_FUNCTION_INFO_V1 (molecule_in_varchar);

Datum
molecule_in_varchar (PG_FUNCTION_ARGS)
{
    VarChar *x = PG_GETARG_VARCHAR_P (0);
    char *input = VARDATA(x);
    int size = VARSIZE(x)-VARHDRSZ;

    PG_RETURN_MOLECULE_P (make_molecule(input,size));
}

PG_FUNCTION_INFO_V1 (molecule_in_bytea);

Datum
molecule_in_bytea (PG_FUNCTION_ARGS)
{
    bytea *x = PG_GETARG_BYTEA_P (0);
    char *input = VARDATA(x);
    int size = VARSIZE(x)-VARHDRSZ;

    PG_RETURN_MOLECULE_P (make_molecule(input,size));
}

/*
* Output a molecule in cstring form as molfile
*/
PG_FUNCTION_INFO_V1 (molecule_out);

Datum
molecule_out (PG_FUNCTION_ARGS)
{
    MOLECULE *molecule = PG_GETARG_MOLECULE_P (0);
    DECOMPRESSED_DATA *original_data=NULL;

    char *result=NULL;

    if(molecule->original_format == FORMAT_SMILES)
    {
        result = (char *) palloc0 (molecule->sizesmi);

        //memset(result,0x0,molecule->sizesmi);

        strncpy (result, SMIPTR(molecule), molecule->sizesmi);
    }
    else
    {
        result = (char *) palloc0 (molecule->sizeo);

        original_data = decompress_data(CIPTR (molecule), molecule->compressed_sizeo, molecule->sizeo);

        //memset(result,0x0,molecule->sizeo);

        strncpy (result, original_data->decompressed_data, molecule->sizeo);
    }

    //for(i=0;i<arr[0];i++) printf("%d ",arr[i]);
    //printf("\n");

    PG_RETURN_CSTRING (result);
}

/*****************************************************************************
 * Binary Input/Output functions
 *
 * These are optional.
 *****************************************************************************/

PG_FUNCTION_INFO_V1 (molecule_recv);

Datum
molecule_recv (PG_FUNCTION_ARGS)
{
    StringInfo buf = (StringInfo) PG_GETARG_POINTER (0);
    int len = buf->len;
    const char *str = pq_getmsgbytes (buf, len);
    MOLECULE *result = (MOLECULE *) palloc0 (len);

    //SET_VARSIZE (result,(buf->len + VARHDRSZ));

    //memset(result,0x0,len);

    memcpy (result, str, len);

    PG_RETURN_POINTER (result);
}

PG_FUNCTION_INFO_V1 (molecule_send);

Datum
molecule_send (PG_FUNCTION_ARGS)
{
    MOLECULE *molecule = PG_GETARG_MOLECULE_P (0);

    /*StringInfoData buf;

    pq_begintypsend(&buf);

    pq_sendbytes(&buf,(const char*) molecule,sizeof(molecule));

    PG_RETURN_BYTEA_P(pq_endtypsend(&buf));*/

    PG_RETURN_BYTEA_P(molecule);
}

/*
* Strip smaller fragments from a molecule, leaving only the larggest one
*/
PG_FUNCTION_INFO_V1 (pgchem_strip_salts);

Datum
pgchem_strip_salts (PG_FUNCTION_ARGS)
{
    char *molfile = NULL;
    char *smiles = NULL;
    //unsigned int *efa_array = NULL;
    MOLECULE *retval;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
    bool neutralize_residue = PG_GETARG_BOOL (1);

    smiles = ob_strip_salts (SMIPTR (arg_molecule), neutralize_residue ? 1 : 0);

    molfile = ob_smiles_to_V2000 (smiles);

    if(molfile == NULL || !strlen(molfile))
    {
        elog (ERROR, "Molfile generation failed! Offender was :\n %s",smiles);
    }

    //efa_array = ob_efa_array(smiles);

    retval = new_molecule (smiles, molfile);

    free (molfile);
    free (smiles);
    //free (efa_array);

    PG_RETURN_MOLECULE_P (retval);
}

/*
* Add hydrogens to a molecule
*/
PG_FUNCTION_INFO_V1 (pgchem_add_hydrogens);

Datum
pgchem_add_hydrogens (PG_FUNCTION_ARGS)
{
    char *molfile = NULL;
    char *smiles = NULL;
    //unsigned int *efa_array = NULL;
    MOLECULE *retval;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
    bool polaronly = PG_GETARG_BOOL (1);
    bool correct4PH = PG_GETARG_BOOL (2);
    float8 PH = PG_GETARG_FLOAT8 (3);

    smiles =
        ob_add_hydrogens (SMIPTR (arg_molecule),
                          polaronly ? 1 : 0, correct4PH ? 1 : 0, PH);

    molfile = ob_smiles_to_V2000 (smiles);

    if(molfile == NULL || !strlen(molfile))
    {
        elog (ERROR, "Molfile generation failed! Offender was :\n %s",smiles);
    }

    //efa_array = ob_efa_array(smiles);

    retval = new_molecule (smiles, molfile);

    free (molfile);
    free (smiles);
    //free (efa_array);

    PG_RETURN_MOLECULE_P (retval);
}

/*
* Remove hydrogens from a molecule
*/
PG_FUNCTION_INFO_V1 (pgchem_remove_hydrogens);

Datum
pgchem_remove_hydrogens (PG_FUNCTION_ARGS)
{
    char *molfile = NULL;
    char *smiles = NULL;
    //unsigned int *efa_array = NULL;
    MOLECULE *retval;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);
    bool nonpolaronly = PG_GETARG_BOOL (1);

    smiles =
        ob_delete_hydrogens (SMIPTR (arg_molecule),
                             nonpolaronly ? 1 : 0);

    molfile = ob_smiles_to_V2000 (smiles);

    if(molfile == NULL || !strlen(molfile))
    {
        elog (ERROR, "Molfile generation failed! Offender was :\n %s",smiles);
    }
    //efa_array = ob_efa_array(smiles);

    retval = new_molecule (smiles, molfile);

    free (molfile);
    free (smiles);
    //free (efa_array);

    PG_RETURN_MOLECULE_P (retval);
}

PG_FUNCTION_INFO_V1 (pgchem_original_format);

Datum
pgchem_original_format (PG_FUNCTION_ARGS)
{
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

    switch (arg_molecule->original_format)
    {
    case FORMAT_V3000:
        PG_RETURN_CSTRING ("V3000");
        break;
    case FORMAT_V2000:
        PG_RETURN_CSTRING ("V2000");
        break;
    case FORMAT_INCHI:
        PG_RETURN_CSTRING ("INCHI");
        break;
    case FORMAT_SMILES:
        PG_RETURN_CSTRING ("SMILES");
        break;
    default:
        PG_RETURN_CSTRING ("UNKNOWN");
    }
}
