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
#ifdef BUILD_WITH_INDIGO
#include "mingw/indigo.h"
#endif

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
#ifdef BUILD_WITH_INDIGO
    int i=IN_SUBFPSIZE/2;
#else
    int i= OB_FPSIZE/2;
#endif

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
#ifdef BUILD_WITH_INDIGO
    int i=IN_SUBFPSIZE/2;
#else
    int i= OB_FPSIZE/2;
#endif

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

#ifdef BUILD_WITH_INDIGO
Datum
molecule_contained_in (PG_FUNCTION_ARGS)
{
    MOLECULE *query = PG_GETARG_MOLECULE_P (0);
    MOLECULE *predicate = PG_GETARG_MOLECULE_P (1);
    int matchhandle;
    int matcherhandle, molhandle, queryhandle;
    bool match;



    /*elog (INFO, "SMILES: %s\n",SMIPTR(query));
    elog (INFO, "SMILES: %s\n",SMIPTR(predicate));
    elog (INFO, "MF: %s\n",MFPTR(query));
    elog (INFO, "MF: %s\n",MFPTR(predicate));*/

    molhandle = indigoUnserialize((unsigned char*) ANCPTR(predicate), predicate->sizeanc);
    //molhandle = indigoLoadMoleculeFromString(SMIPTR(predicate));

    //elog (INFO, "SMILES: %i\n",query->sizeanc);
    //elog (INFO, "SMILES: %i\n",predicate->sizeanc);
    //elog (INFO, "SMILES: %i\n",ANCPTR(predicate)[predicate->sizeanc]);

    //if(molhandle<0) elog (INFO, "1. Error was: %s",indigoGetLastError());

    queryhandle = indigoLoadQueryMoleculeFromString((const char*) SMIPTR(query));

    indigoAromatize(queryhandle);

    //if(queryhandle<0) elog (INFO, "2. Error was: %s",indigoGetLastError());

    matcherhandle = indigoSubstructureMatcher(molhandle,NULL);

    //if(matcherhandle<0) elog (INFO, "3. Error was: %s",indigoGetLastError());

    matchhandle = indigoMatch(matcherhandle,queryhandle);

    //if(matchhandle<0) elog (INFO, "4. Error was: %s",indigoGetLastError());

    match = (matchhandle != 0);

    indigoFree(molhandle);
    indigoFree(queryhandle);
    indigoFree(matcherhandle);
    indigoFree(matchhandle);

    PG_RETURN_BOOL (match);
}


Datum
molecule_contains (PG_FUNCTION_ARGS)
{
    MOLECULE *query = PG_GETARG_MOLECULE_P (1);
    MOLECULE *predicate = PG_GETARG_MOLECULE_P (0);
    int matchhandle;
    int matcherhandle, molhandle, queryhandle;
    bool match;



    /*elog (INFO, "SMILES: %s\n",SMIPTR(query));
    elog (INFO, "SMILES: %s\n",SMIPTR(predicate));
    elog (INFO, "MF: %s\n",MFPTR(query));
    elog (INFO, "MF: %s\n",MFPTR(predicate));*/

    molhandle = indigoUnserialize((unsigned char*) ANCPTR(predicate), predicate->sizeanc);
    //molhandle = indigoLoadMoleculeFromString(SMIPTR(predicate));

    //elog (INFO, "SMILES: %i\n",query->sizeanc);
    //elog (INFO, "SMILES: %i\n",predicate->sizeanc);
    //elog (INFO, "SMILES: %i\n",ANCPTR(predicate)[predicate->sizeanc]);

    //if(molhandle<0) elog (INFO, "1. Error was: %s",indigoGetLastError());

    queryhandle = indigoLoadQueryMoleculeFromString((const char*) SMIPTR(query));

    indigoAromatize(queryhandle);

    //if(queryhandle<0) elog (INFO, "2. Error was: %s",indigoGetLastError());

    matcherhandle = indigoSubstructureMatcher(molhandle,NULL);

    //if(matcherhandle<0) elog (INFO, "3. Error was: %s",indigoGetLastError());

    matchhandle = indigoMatch(matcherhandle,queryhandle);

    //if(matchhandle<0) elog (INFO, "4. Error was: %s",indigoGetLastError());

    match = (matchhandle != 0);

    indigoFree(molhandle);
    indigoFree(queryhandle);
    indigoFree(matcherhandle);
    indigoFree(matchhandle);

    PG_RETURN_BOOL (match);
}

/*
* Check if a query molecule is equal to a predicate molecule by comparing their molhashes.
*/
Datum
molecule_equals (PG_FUNCTION_ARGS)
{
    MOLECULE *query = PG_GETARG_MOLECULE_P (0);
    MOLECULE *predicate = PG_GETARG_MOLECULE_P (1);
    //int molhandle, queryhandle, matchhandle;
    //bool match;



    if (memcmp (query->inchikey, predicate->inchikey, INCHIKEYSZ) != 0)
        PG_RETURN_BOOL (false);


    /*molhandle = indigoUnserialize((unsigned char*) ANCPTR(predicate), predicate->sizeanc);

    //if(molhandle<0) elog (INFO, "1. Error was: %s",indigoGetLastError());

    queryhandle = indigoUnserialize((unsigned char*) ANCPTR(query), query->sizeanc);

    //if(queryhandle<0) elog (INFO, "2. Error was: %s",indigoGetLastError());

    matchhandle = indigoExactMatch(molhandle,queryhandle,"NONE");

    match = (matchhandle != 0);

    indigoFree(molhandle);
    indigoFree(queryhandle);
    indigoFree(matchhandle);*/

    PG_RETURN_BOOL (true);
}
#else
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



    //printf("Query: %s\n Pred: %s\n",query->data, predicate->data);

    /*  if (query->isbz == true)
        {
          if (predicate->nobz == true)
    	PG_RETURN_BOOL (false);
          else
    	PG_RETURN_BOOL (true);
        } */
    //if (strchr (querysmi, '@') != NULL)
    /*if(1==2)
    {
        if (query->disconnected == 1)
        {
            match = ob_SSS (querysmi, SMIPTR(predicate));
        }
        else
        {
            match = ob_SSS_SMARTS_native (querysmi, SMIPTR(predicate));
        }
    }
    else
    {*/
    if (query->disconnected == 1)
    {
        match = ob_SSS_bin (querysmi, ANCPTR(predicate));
    }
    else
    {
        match = ob_SSS_SMARTS_native_bin (querysmi, ANCPTR(predicate));
    }
    //}
    //if (query->disconnected == true)
    //{
    //  elog (ERROR, "Disconnected molecules as query input are not supported!");

    /*xm_set_ring_perception_algorithm (RPA_SAR);
    xm_set_strict_typing (FEATURE_OFF);
    mm_set_r_s_check (FEATURE_OFF);
    mm_set_e_z_check (FEATURE_OFF);
    mm_set_chg_check (FEATURE_ON);
    mm_set_iso_check (FEATURE_ON);
    mm_set_rad_check (FEATURE_ON);
    mm_set_exact_match (FEATURE_OFF);

    mm_set_mol (MFPTR(query));
    mm_set_current_mol_as_query ();

    printf("%s\n",SMIPTR(query));

    mm_set_mol (MFPTR(predicate));

    printf("%s\n",SMIPTR(predicate));

    if (mm_match () != 0) {
        printf("Match TRUE!\n");
    	PG_RETURN_BOOL (true);
    }
        }
      else
        {*/

    /*  efa_q = (unsigned int*) EFAPTR(query);
      efa_p = (unsigned int*) EFAPTR(predicate);

     if(efa_q[0] > efa_p[0]) {printf("EFA kill type 1"); PG_RETURN_BOOL (false);}
      else {
          for(i=1;i<efa_q[0];i+=2) {
              match = 0;
              for(j=offset;j<efa_p[0];j+=2) {
                  if(efa_q[i]==efa_p[j] && efa_q[++i] <= efa_p[++j]) {match = 1; offset=j; break;}
                  //printf("%d %d %d %d\n",efa_q[i],efa_p[j],efa_q[i+1],efa_p[j+1]);
              }
              if(match==0) {printf("EFA kill type 2");PG_RETURN_BOOL (false);}
          }
      }  */

    //if(strstr(querysmi,"@")!=NULL || strstr(querysmi,"/")!=NULL || strstr(querysmi,"\\")!=NULL)
    //match = ob_SSS (SMIPTR(query), MFPTR(predicate));
    //else
    //match = ob_SSS_SMARTS_native_bin (querysmi, ANCPTR(predicate));

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


//printf("Query: %s\n Pred: %s\n",query->data, predicate->data);

    /* if (query->isbz == true)
       {
         if (predicate->nobz == true)
    PG_RETURN_BOOL (false);
         else
    PG_RETURN_BOOL (true);
       } */

    //if (strchr (querysmi, '@') != NULL)
    /*if(1==2)
    {
        if (query->disconnected == 1)
        {
            match = ob_SSS (querysmi, SMIPTR(predicate));
        }
        else
        {
            match = ob_SSS_SMARTS_native (querysmi, SMIPTR(predicate));
        }
    }
    else
    {*/
    if (query->disconnected == 1)
    {
        match = ob_SSS_bin (querysmi, ANCPTR(predicate));
    }
    else
    {
        match = ob_SSS_SMARTS_native_bin (querysmi, ANCPTR(predicate));
    }
    //}
    //  elog (ERROR, "Disconnected molecules as query input are not supported!");

    /*xm_set_ring_perception_algorithm (RPA_SAR);
    xm_set_strict_typing (FEATURE_OFF);
    mm_set_r_s_check (FEATURE_OFF);
    mm_set_e_z_check (FEATURE_OFF);
    mm_set_chg_check (FEATURE_ON);
    mm_set_iso_check (FEATURE_ON);
    mm_set_rad_check (FEATURE_ON);
    mm_set_exact_match (FEATURE_OFF);

         mm_set_mol (MFPTR(query));
    mm_set_current_mol_as_query ();

    printf("%s\n",SMIPTR(query));

    mm_set_mol (MFPTR(predicate));

     printf("%s\n",SMIPTR(predicate));

    if (mm_match () != 0) {
    printf("Match TRUE!\n");
    	PG_RETURN_BOOL (true);}
        }
      else
        {*/

    /*efa_q = (unsigned int*) EFAPTR(query);
    efa_p = (unsigned int*) EFAPTR(predicate);

    if(efa_q[0] > efa_p[0]) {printf("EFA kill type 1"); PG_RETURN_BOOL (false);}
    else {
        for(i=1;i<efa_q[0];i+=2) {
            match = 1;
            for(j=offset;j<efa_p[0];j+=2) {
                if(efa_q[i] == efa_p[j] && efa_q[i+1] <= efa_p[j+1]) {match = 1; offset=j; break;}

               // printf("%d %d %d %d\n",efa_q[i],efa_p[j],efa_q[i+1],efa_p[j+1]);
            }
            if(match==0) {printf("EFA kill type 2");PG_RETURN_BOOL (false);}
        }
    } */


    // if(strstr(querysmi,"@")!=NULL || strstr(querysmi,"/")!=NULL || strstr(querysmi,"\\")!=NULL)
    // match = ob_SSS (SMIPTR(query), MFPTR(predicate));
    //else
    //match = ob_SSS_SMARTS_native_bin (querysmi, ANCPTR(predicate));

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


    /* if (query->isbz == true) {
         if (predicate->isbz == true)
       PG_RETURN_BOOL (true);
         else
       PG_RETURN_BOOL (false);
       }    */
    //if(query->popcount != predicate->popcount) PG_RETURN_BOOL (false);

    if (memcmp (query->inchikey, predicate->inchikey, INCHIKEYSZ) != 0)
        PG_RETURN_BOOL (false);

    //Further matching could be done here:
    //OBQuery *query = CompileSmilesQuery("c1ccccc1");
//OBIsomorphismMapper *mapper = OBIsomorphismMapper::GetInstance(query);
//OBIsomorphismMapper::Mappings maps = mapper->MapUnique(mol);


    PG_RETURN_BOOL (true);
}
#endif

/*
* Returns the Tanimoto similarity of two molecules.
*/
Datum
molecule_similarity (PG_FUNCTION_ARGS)
{
    MOLECULE *mol1 = PG_GETARG_MOLECULE_P (0);
    MOLECULE *mol2 = PG_GETARG_MOLECULE_P (1);

#ifdef BUILD_WITH_INDIGO
    PG_RETURN_FLOAT8 (ob_tanimoto
                      (mol1->fp.dwords+IN_SIMOFFSET,  mol2->fp.dwords+IN_SIMOFFSET, IN_SIMFPSIZE));
#else
    PG_RETURN_FLOAT8 (ob_tanimoto
                      (mol1->fp.dwords,  mol2->fp.dwords, OB_FPSIZE2));
#endif
}

/*
* Returns the Tanimoto similarity of two molecules with index support.
*/
/* Datum
molecule_similarity_gist (PG_FUNCTION_ARGS)
{
    MOLECULE *mol1 = PG_GETARG_MOLECULE_P (0);
    HeapTupleHeader  t = PG_GETARG_HEAPTUPLEHEADER(1);
    bool isnull;
    float8 similarity;
    char *oper;
    float8 threshold;
    MOLECULE *mol2 = (MOLECULE*) DatumGetPointer(GetAttributeByName(t, "_mq", &isnull));
    if(isnull) elog (ERROR, "Query molecule must not be NULL");
    oper = (char*) DatumGetPointer(GetAttributeByName(t, "_op", &isnull));
    if(isnull) elog (ERROR, "Query operator must not be NULL");
    threshold = DatumGetFloat8(GetAttributeByName(t, "_threshold", &isnull));
    if(isnull) elog (ERROR, "Query threshold must not be NULL");

#ifdef BUILD_WITH_INDIGO
    similarity = ob_tanimoto
                 (mol1->fp.dwords+IN_SIMOFFSET,  mol2->fp.dwords+IN_SIMOFFSET, IN_SIMFPSIZE);
#else
    similarity = ob_tanimoto
                 (mol1->fp.dwords,  mol2->fp.dwords, OB_FPSIZE2);
#endif

    if(oper[0]=='>' && oper[1]=='=') PG_RETURN_BOOL(similarity >= threshold);
    if(oper[0]=='<' && oper[1]=='=') PG_RETURN_BOOL(similarity <= threshold);
    if(oper[1]=='>') PG_RETURN_BOOL(similarity > threshold);
    if(oper[1]=='<') PG_RETURN_BOOL(similarity < threshold);
    if(oper[1]=='=') PG_RETURN_BOOL(similarity == threshold);

    PG_RETURN_BOOL(false);
} */
