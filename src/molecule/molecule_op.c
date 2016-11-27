/************************************************************************
 * molecule_op.c molecule operator support functions
 *
 * Copyright (c) 2007,2013 by Ernst-G. Schmid
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
#include "executor/executor.h"
#include "fmgr.h"
#include "molecule.h"
#include "barsoi/barsoi.h"
#include "obwrapper.h"


PG_FUNCTION_INFO_V1 (molecule_maybe_contained_in);
PG_FUNCTION_INFO_V1 (molecule_maybe_contains);
PG_FUNCTION_INFO_V1 (molecule_contained_in);
PG_FUNCTION_INFO_V1 (molecule_contains);
PG_FUNCTION_INFO_V1 (molecule_equals);
PG_FUNCTION_INFO_V1 (molecule_similarity);
//PG_FUNCTION_INFO_V1 (molecule_similarity_gist);

Datum molecule_maybe_contained_in (PG_FUNCTION_ARGS);
Datum molecule_maybe_contains (PG_FUNCTION_ARGS);
Datum molecule_contained_in (PG_FUNCTION_ARGS);
Datum molecule_contains (PG_FUNCTION_ARGS);
Datum molecule_equals (PG_FUNCTION_ARGS);
Datum molecule_similarity (PG_FUNCTION_ARGS);
//Datum molecule_similarity_gist (PG_FUNCTION_ARGS);

/*
* Always returns TRUE for debugging purposes
*/
/*Datum
molecule_alwaystrue (PG_FUNCTION_ARGS)
{
    PG_RETURN_BOOL (true);
}*/

Datum
molecule_maybe_contained_in (PG_FUNCTION_ARGS)
{
    MOLECULE *query = PG_GETARG_MOLECULE_P (0);
    MOLECULE *predicate = PG_GETARG_MOLECULE_P (1);
    MOLFP allzero;
    uint64 *zero;

    int i= OB_FPSIZE/2;

    allzero.v = (predicate->fp.v & query->fp.v) ^ query->fp.v;

    zero = &allzero.qwords[0];

    while(i--)
    {
        if (*zero != 0x0ULL)
        {
            PG_RETURN_BOOL (false);
        }
        zero++;
    }

    //elog(INFO,"true");
    PG_RETURN_BOOL (true);
}

Datum
molecule_maybe_contains (PG_FUNCTION_ARGS)
{
    MOLECULE *predicate = PG_GETARG_MOLECULE_P (0);
    MOLECULE *query = PG_GETARG_MOLECULE_P (1);
    MOLFP allzero;
    uint64 *zero;

    int i= OB_FPSIZE/2;

    allzero.v = (predicate->fp.v & query->fp.v) ^ query->fp.v;

    zero = &allzero.qwords[0];

    while(i--)
    {
        if (*zero != 0x0ULL)
        {
            PG_RETURN_BOOL (false);
        }
        zero++;
    }

    //elog(INFO,"true");
    PG_RETURN_BOOL (true);
}


/*
* Check if a query molecule is contained in a predicate molecule by performing a graph isomorphism check.
* If the query molecule is disconnected (has fragments), the check is done with checkmol/barsoi, bacause the OpenBabel matcher does not support this.
* Otherwise the faster openBabel matcher is used.
*/
Datum
molecule_contained_in (PG_FUNCTION_ARGS)
{
    MOLECULE *query = PG_GETARG_MOLECULE_P (0);
    MOLECULE *predicate = PG_GETARG_MOLECULE_P (1);
    int match=0;//, i, j, offset = 1;
    char* querysmi = SMIPTR(query);
    //unsigned int *efa_q, *efa_p;


    if (query->disconnected == 1)
    {
        match = ob_SSS (querysmi, SMIPTR(predicate));
    }
    else
    {
        match = ob_SSS_SMARTS_native (querysmi, SMIPTR(predicate));
    }

    if(match==-1) elog (ERROR, "Invalid SMARTS pattern: %s",SMIPTR(query));
    if(match==-2) elog (ERROR, "Deserialize failed");
    if(match==-3) elog (ERROR, "Mol empty");

    if (match != 0)
        PG_RETURN_BOOL (true);
    //}
//printf("EFA kill type 3\n");
    PG_RETURN_BOOL (false);
}

/*
* Check if a query molecule is contained in a predicate molecule by performing a graph isomorphism check.
* If the query molecule is disconnected (has fragments), the check is done with checkmol/barsoi, bacause the OpenBabel matcher does not support this.
* Otherwise the faster openBabel matcher is used.
*/
Datum
molecule_contains (PG_FUNCTION_ARGS)
{
    MOLECULE *query = PG_GETARG_MOLECULE_P (1);
    MOLECULE *predicate = PG_GETARG_MOLECULE_P (0);
    int match=0;
    char* querysmi = SMIPTR(query);
    //,i,j, offset=1;
    //unsigned int *efa_q, *efa_p;

    if (query->disconnected == 1)
    {
        match = ob_SSS (querysmi, SMIPTR(predicate));
    }
    else
    {
        match = ob_SSS_SMARTS_native (querysmi, SMIPTR(predicate));
    }

    if(match==-1) elog (ERROR, "Invalid SMARTS pattern: %s",SMIPTR(query));
    if(match==-2) elog (ERROR, "Deserialize failed");
    if(match==-3) elog (ERROR, "Mol empty");

    if (match != 0)
        PG_RETURN_BOOL (true);
    //}

//printf("EFA kill type 3\n");
    PG_RETURN_BOOL (false);
}

/*
* Check if a query molecule is equal to a predicate molecule by comparing their molhashes.
*/
Datum
molecule_equals (PG_FUNCTION_ARGS)
{
    MOLECULE *query = PG_GETARG_MOLECULE_P (0);
    MOLECULE *predicate = PG_GETARG_MOLECULE_P (1);

    if (memcmp (query->inchikey, predicate->inchikey, INCHIKEYSZ) != 0)
        PG_RETURN_BOOL (false);

    //Further matching could be done here:
    //OBQuery *query = CompileSmilesQuery("c1ccccc1");
//OBIsomorphismMapper *mapper = OBIsomorphismMapper::GetInstance(query);
//OBIsomorphismMapper::Mappings maps = mapper->MapUnique(mol);


    PG_RETURN_BOOL (true);
}


/*
* Returns the Tanimoto similarity of two molecules.
*/
Datum
molecule_similarity (PG_FUNCTION_ARGS)
{
    MOLECULE *mol1 = PG_GETARG_MOLECULE_P (0);
    MOLECULE *mol2 = PG_GETARG_MOLECULE_P (1);

    PG_RETURN_FLOAT8 (ob_tanimoto
                      (mol1->fp.dwords,  mol2->fp.dwords, OB_FPSIZE2));
}
