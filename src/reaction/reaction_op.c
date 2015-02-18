/************************************************************************
 * reaction_op.c reaction operator support functions
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
#include "reaction.h"
//#include "barsoi/barsoi.h"
#include "obwrapper.h"
#include "fmgr.h"

#define CHECK_BIT(var,pos) ((var) & (1<<(pos)))
#define SET_BIT(var,pos) ((var) |= (1<<(pos)))
#define MAX_ELEM 32

PG_FUNCTION_INFO_V1 (reaction_alwaystrue);
PG_FUNCTION_INFO_V1 (reaction_contained_in);
PG_FUNCTION_INFO_V1 (reaction_contains);
PG_FUNCTION_INFO_V1 (reaction_equals);
//PG_FUNCTION_INFO_V1 (reaction_equals_exact);
//PG_FUNCTION_INFO_V1 (reaction_equals_products_exact);
PG_FUNCTION_INFO_V1 (reaction_similarity);
PG_FUNCTION_INFO_V1 (reaction_similarity_reactants);
PG_FUNCTION_INFO_V1 (reaction_similarity_products);

Datum reaction_alwaystrue (PG_FUNCTION_ARGS);
Datum reaction_contained_in (PG_FUNCTION_ARGS);
Datum reaction_contains (PG_FUNCTION_ARGS);
//Datum reaction_equals_exact (PG_FUNCTION_ARGS);
//Datum reaction_equals_products_exact (PG_FUNCTION_ARGS);
Datum reaction_equals (PG_FUNCTION_ARGS);
Datum reaction_similarity (PG_FUNCTION_ARGS);
Datum reaction_similarity_reactants (PG_FUNCTION_ARGS);
Datum reaction_similarity_products (PG_FUNCTION_ARGS);

inline static Datum rss_match(REACTION *query, REACTION *predicate)
{
    MOLECULE *tmpMol_q = NULL, *tmpMol_p = NULL;
    char *offset_q = MOLARRAYPTR(query);
    char *offset_p = MOLARRAYPTR(predicate);
    char *offset_q_products = NULL;
    int i,j,m=0;
    int match;
    int32 q_num_products = query->num_products;
    int32 q_num_reactants = query->num_reactants;
    int32 p_num_products = predicate->num_products;
    int32 p_num_reactants = predicate->num_reactants;
    //int q_r_matches[p_num_reactants];
    //int q_p_matches[p_num_products];
    //int p_r_matches[p_num_reactants];
    //int p_p_matches[p_num_products];
    uint32 q_r_matches=0x0;
    uint32 q_p_matches=0x0;
    uint32 p_r_matches=0x0;
    uint32 p_p_matches=0x0;

    if(q_num_products > MAX_ELEM || p_num_products > MAX_ELEM || q_num_reactants > MAX_ELEM || p_num_reactants > MAX_ELEM)
    {
        elog (WARNING, "Only partial reaction matching");
        if(q_num_products > MAX_ELEM) q_num_products = MAX_ELEM;
        if(p_num_products > MAX_ELEM) p_num_products = MAX_ELEM;
        if(q_num_reactants > MAX_ELEM) q_num_reactants = MAX_ELEM;
        if(p_num_reactants > MAX_ELEM) p_num_reactants = MAX_ELEM;
    }

    if(q_num_products > p_num_products || q_num_reactants > p_num_reactants) PG_RETURN_BOOL (false);

    if(q_num_products +  p_num_products + q_num_reactants + p_num_reactants == 0) PG_RETURN_BOOL (true);

    if(q_num_products + q_num_reactants == 0) PG_RETURN_BOOL (false);

    //memset(&q_r_matches,0x0,q_num_reactants*sizeof(int));
    //memset(&q_p_matches,0x0,q_num_products*sizeof(int));
    //memset(&p_r_matches,0x0,p_num_reactants*sizeof(int));
    //memset(&p_p_matches,0x0,p_num_products*sizeof(int));

    if (q_num_reactants != 0)
    {
        for (i=0; i<p_num_reactants; i++)
        {

            offset_q = MOLARRAYPTR(query);
            tmpMol_p = (MOLECULE*) offset_p;

            for(j=0; j<q_num_reactants; j++)
            {

                tmpMol_q = (MOLECULE*) offset_q;

                if (tmpMol_q->disconnected == true)
                {
                    match = ob_SSS (SMIPTR(tmpMol_q), SMIPTR(tmpMol_p));
                }
                else
                {
                    match = ob_SSS_SMARTS_native (SMIPTR(tmpMol_q), SMIPTR(tmpMol_p));
                }

                if(match<0) elog (ERROR, "Invalid SMARTS pattern: %s",SMIPTR(tmpMol_q));

                if ( match != 0)
                {
                    //p_r_matches[i]=1;
                    //q_r_matches[j]=1;
                    SET_BIT(p_r_matches,i);
                    SET_BIT(q_r_matches,j);
                }

                offset_q+=VARSIZE(tmpMol_q)*sizeof(char);
            }

            offset_p+=VARSIZE(tmpMol_p)*sizeof(char);
        }
    }
    else
    {
        for (i=0; i<p_num_reactants; i++)
        {
            tmpMol_p = (MOLECULE*) offset_p;
            offset_p+=VARSIZE(tmpMol_p)*sizeof(char);  //Must advance to products in any case
        }
    }

//j=0;

    for(i=0; i<p_num_reactants; i++)
    {
        //if(p_r_matches[i] == 1) j++;
        //m+=p_r_matches[i];
        if(CHECK_BIT(p_r_matches,i)) m++;
    }

    for(i=0; i<q_num_reactants; i++)
    {
        //if(p_r_matches[i] == 1) j++;
        //m+=q_r_matches[i];
        if(CHECK_BIT(q_r_matches,i)) m++;
    }

    if(m<(q_num_reactants*2)) PG_RETURN_BOOL (false);

    offset_q_products = offset_q;
    //offset_p = MOLARRAYPTR(predicate)+p_num_reactants*sizeof(MOLECULE*);

    if(q_num_products != 0)
    {
        for (i=0; i<p_num_products; i++)
        {

            //offset_q = MOLARRAYPTR(query);
            offset_q = offset_q_products;
            tmpMol_p = (MOLECULE*) offset_p;

            /*for(j=0;j<q_num_reactants;j++) {
                tmpMol_q = (MOLECULE*) offset_q;
                offset_q+=tmpMol_q->len*sizeof(char);
                } */

            for(j=0; j<q_num_products; j++)
            {

                tmpMol_q = (MOLECULE*) offset_q;

                if (tmpMol_q->disconnected == true)
                {
                    match = ob_SSS (SMIPTR(tmpMol_q), SMIPTR(tmpMol_p));
                }
                else
                {
                    match = ob_SSS_SMARTS_native (SMIPTR(tmpMol_q), SMIPTR(tmpMol_p));
                }

                if(match<0) elog (ERROR, "Invalid SMARTS pattern: %s",SMIPTR(tmpMol_q));

                if (match != 0)
                {
                    //p_p_matches[i]=1;
                    //q_p_matches[j]=1;
                    SET_BIT(p_p_matches,i);
                    SET_BIT(q_p_matches,j);
                }

                offset_q+=VARSIZE(tmpMol_q)*sizeof(char);
            }

            offset_p+=VARSIZE(tmpMol_p)*sizeof(char);
        }
    }

    m=0;

    for(i=0; i<p_num_products; i++)
    {
        //if(p_matches[i] == 1) j++;
        //m+=p_p_matches[i];
        if(CHECK_BIT(p_p_matches,i)) m++;
    }

    for(i=0; i<q_num_products; i++)
    {
        //if(p_r_matches[i] == 1) j++;
        //m+=q_p_matches[i];
        if(CHECK_BIT(q_p_matches,i)) m++;
    }

    if(m<(q_num_products*2)) PG_RETURN_BOOL (false);

    PG_RETURN_BOOL (true);
}

/*
* Always returns TRUE for debugging purposes
*/
Datum
reaction_alwaystrue (PG_FUNCTION_ARGS)
{
    PG_RETURN_BOOL (true);
}

/*
* Check if a query reaction is contained in a predicate reaction by performing a graph isomorphism check.
* If the query reaction is disconnected (has fragments), the check is done with checkmol/barsoi, bacause the OpenBabel matcher does not support this.
* Otherwise the faster openBabel matcher is used.
*/
Datum
reaction_contained_in (PG_FUNCTION_ARGS)
{
    REACTION *query = PG_GETARG_REACTION_P (0);
    REACTION *predicate = PG_GETARG_REACTION_P (1);


    PG_RETURN_BOOL (rss_match(query,predicate));
}

/*
* Check if a query reaction is contained in a predicate reaction by performing a graph isomorphism check.
* If the query reaction is disconnected (has fragments), the check is done with checkmol/barsoi, bacause the OpenBabel matcher does not support this.
* Otherwise the faster openBabel matcher is used.
*/
Datum
reaction_contains (PG_FUNCTION_ARGS)
{
    REACTION *query = PG_GETARG_REACTION_P (1);
    REACTION *predicate = PG_GETARG_REACTION_P (0);

    PG_RETURN_BOOL (rss_match(query,predicate));
}

/*
* Check if a query reaction is exactly equal to a predicate reaction by comparing their hashes and positions.
*/
//Datum
//reaction_equals_exact (PG_FUNCTION_ARGS)
//{
//  REACTION *query = PG_GETARG_REACTION_P (0);
//  REACTION *predicate = PG_GETARG_REACTION_P (1);
//  MOLECULE *tmpMol_q, *tmpMol_p;
//  char *offset_q = MOLARRAYPTR(query);
//  char *offset_p = MOLARRAYPTR(predicate);
//  int i;
//
//  if(query->num_products != predicate->num_products || query->num_reactants != predicate->num_reactants) PG_RETURN_BOOL (false);
//
//  for (i=0;i<query->num_products+query->num_reactants;i++) {
//      tmpMol_q = (MOLECULE*) offset_q;
//      tmpMol_p = (MOLECULE*) offset_p;
//
//      if (memcmp (tmpMol_q->inchikey, tmpMol_p->inchikey, INCHIKEYSZ) != 0) PG_RETURN_BOOL (false);
//
//      offset_q+=VARSIZE(tmpMol_q)*sizeof(char);
//      offset_p+=VARSIZE(tmpMol_p)*sizeof(char);
//  }
//
//  PG_RETURN_BOOL (true);
//}

/*
* Check if a query reaction is exactly equal on the product side to a predicate reaction by comparing all hashes and positions of the products.
*/
//Datum
//reaction_equals_products_exact (PG_FUNCTION_ARGS)
//{
//  REACTION *query = PG_GETARG_REACTION_P (0);
//  REACTION *predicate = PG_GETARG_REACTION_P (1);
//  MOLECULE *tmpMol_q = NULL, *tmpMol_p = NULL;
//  char *offset_q = MOLARRAYPTR(query);
//  char *offset_p = MOLARRAYPTR(predicate);
//  int i,j;
//  int r_matches[predicate->num_reactants];
//
//  if(query->num_products != predicate->num_products || query->num_reactants != predicate->num_reactants) PG_RETURN_BOOL (false);
//
//  memset(&r_matches,0x0,predicate->num_reactants*sizeof(int));
//
//  if (query->num_reactants != 0) {
//  for (i=0;i<predicate->num_reactants;i++) {
//
//      offset_q = MOLARRAYPTR(query);
//      tmpMol_p = (MOLECULE*) offset_p;
//
//      for(j=0;j<query->num_reactants;j++) {
//
//            tmpMol_q = (MOLECULE*) offset_q;
//
//            if (memcmp (tmpMol_q->inchikey, tmpMol_p->inchikey, INCHIKEYSZ) == 0) {
//                r_matches[i]=1;
//            }
//
//            offset_q+=VARSIZE(tmpMol_q)*sizeof(char);
//            }
//
//      offset_p+=VARSIZE(tmpMol_p)*sizeof(char);
//  }
//}  else {
//    for (i=0;i<predicate->num_reactants;i++) {
//        tmpMol_p = (MOLECULE*) offset_p;
//        offset_p+=VARSIZE(tmpMol_p)*sizeof(char);  //Must advance to products in any case
//}
//}
//
// j=0;
//
//  for(i=0;i<predicate->num_reactants;i++) {
//      if(r_matches[i] == 1) j++;
//  }
//
//  if(j!=query->num_reactants) PG_RETURN_BOOL (false);
//
//   /*for(j=0;j<query->num_reactants;j++) {
//          tmpMol_q = (MOLECULE*) offset_q;
//          offset_q+=tmpMol_q->len*sizeof(char);
//          }*/
//
//   if(query->num_products != 0) {
//  for (i=0;i<query->num_products;i++) {
//      tmpMol_q = (MOLECULE*) offset_q;
//      tmpMol_p = (MOLECULE*) offset_p;
//
//      if (memcmp (tmpMol_q->inchikey, tmpMol_p->inchikey, INCHIKEYSZ) != 0) PG_RETURN_BOOL (false);
//
//      offset_q+=VARSIZE(tmpMol_q)*sizeof(char);
//      offset_p+=VARSIZE(tmpMol_p)*sizeof(char);
//  }
//}
//
//  PG_RETURN_BOOL (true);
//}

/*
* Check if a query reaction is equal on the product side to a predicate reaction by comparing all hashes.
*/
Datum
reaction_equals (PG_FUNCTION_ARGS)
{
    REACTION *query = PG_GETARG_REACTION_P (0);
    REACTION *predicate = PG_GETARG_REACTION_P (1);
    MOLECULE *tmpMol_q = NULL, *tmpMol_p = NULL;
    char *offset_q = MOLARRAYPTR(query);
    char *offset_p = MOLARRAYPTR(predicate);
    char *offset_q_products = NULL;
    int i,j,m=0;//,query_len = predicate->num_products+predicate->num_reactants;
    int32 q_num_products = query->num_products;
    int32 q_num_reactants = query->num_reactants;
    int32 p_num_products = predicate->num_products;
    int32 p_num_reactants = predicate->num_reactants;
    //int q_r_matches[p_num_reactants];
    //int q_p_matches[p_num_products];
    //int p_r_matches[p_num_reactants];
    //int p_p_matches[p_num_products];
    uint32 q_r_matches=0x0;
    uint32 q_p_matches=0x0;
    uint32 p_r_matches=0x0;
    uint32 p_p_matches=0x0;

    if(q_num_products > MAX_ELEM || p_num_products > MAX_ELEM || q_num_reactants > MAX_ELEM || p_num_reactants > MAX_ELEM)
    {
        elog (WARNING, "Only partial reaction matching");
        if(q_num_products > MAX_ELEM) q_num_products = MAX_ELEM;
        if(p_num_products > MAX_ELEM) p_num_products = MAX_ELEM;
        if(q_num_reactants > MAX_ELEM) q_num_reactants = MAX_ELEM;
        if(p_num_reactants > MAX_ELEM) p_num_reactants = MAX_ELEM;
    }

    if(q_num_products != p_num_products || q_num_reactants != p_num_reactants) PG_RETURN_BOOL (false);

    //memset(&q_r_matches,0x0,q_num_reactants*sizeof(int));
    //memset(&q_p_matches,0x0,q_num_products*sizeof(int));
    //memset(&p_r_matches,0x0,p_num_reactants*sizeof(int));
    //memset(&p_p_matches,0x0,p_num_products*sizeof(int));

    if (q_num_reactants != 0)
    {
        for (i=0; i<p_num_reactants; i++)
        {

            offset_q = MOLARRAYPTR(query);
            tmpMol_p = (MOLECULE*) offset_p;

            for(j=0; j<q_num_reactants; j++)
            {

                tmpMol_q = (MOLECULE*) offset_q;

                if (memcmp (tmpMol_q->inchikey, tmpMol_p->inchikey, INCHIKEYSZ) == 0)
                {
                    //p_r_matches[i]=1;
                    //q_r_matches[j]=1;
                    SET_BIT(p_r_matches,i);
                    SET_BIT(q_r_matches,j);
                }

                offset_q+=VARSIZE(tmpMol_q)*sizeof(char);
            }

            offset_p+=VARSIZE(tmpMol_p)*sizeof(char);
        }
    }
    else
    {
        for (i=0; i<p_num_reactants; i++)
        {
            tmpMol_p = (MOLECULE*) offset_p;
            offset_p+=VARSIZE(tmpMol_p)*sizeof(char);  //Must advance to products in any case
        }
    }

    for(i=0; i<p_num_reactants; i++)
    {
        //if(p_r_matches[i] == 1) j++;
        //m+=p_r_matches[i];
        if(CHECK_BIT(p_r_matches,i)) m++;
    }

    for(i=0; i<q_num_reactants; i++)
    {
        //if(p_r_matches[i] == 1) j++;
//      m+=q_r_matches[i];
        if(CHECK_BIT(q_r_matches,i)) m++;
    }

    if(m!=(q_num_reactants+p_num_reactants)) PG_RETURN_BOOL (false);

    offset_q_products = offset_q;

    if(q_num_products != 0)
    {
        for (i=0; i<p_num_products; i++)
        {

            offset_q = offset_q_products;
            tmpMol_p = (MOLECULE*) offset_p;

            for(j=0; j<q_num_products; j++)
            {

                tmpMol_q = (MOLECULE*) offset_q;

                if ( memcmp (tmpMol_p->inchikey, tmpMol_p->inchikey, INCHIKEYSZ) == 0)
                {
                    //p_p_matches[i]=1;
                    //q_p_matches[j]=1;
                    SET_BIT(p_p_matches,i);
                    SET_BIT(q_p_matches,j);
                }

                offset_q+=VARSIZE(tmpMol_q)*sizeof(char);
            }

            offset_p+=VARSIZE(tmpMol_p)*sizeof(char);
        }
    }

    m=0;

    for(i=0; i<p_num_products; i++)
    {
        //if(p_matches[i] == 1) j++;
        //m+=p_p_matches[i];
        if(CHECK_BIT(p_p_matches,i)) m++;
    }

    for(i=0; i<q_num_products; i++)
    {
        //if(p_r_matches[i] == 1) j++;
        //m+=q_p_matches[i];
        if(CHECK_BIT(q_p_matches,i)) m++;
    }

    if(m!=(q_num_products+p_num_products)) PG_RETURN_BOOL (false);

    PG_RETURN_BOOL (true);
}

/*
* Returns the Tanimoto similarity of two reactions.
*/
Datum
reaction_similarity (PG_FUNCTION_ARGS)
{
    REACTION *rxn1 = PG_GETARG_REACTION_P (0);
    REACTION *rxn2 = PG_GETARG_REACTION_P (1);


    PG_RETURN_FLOAT8 (ob_tanimoto
                      ((const uint32 *) rxn1->fp,  (const uint32 *) rxn2->fp, 2*OB_FPSIZE2));
}

Datum
reaction_similarity_reactants (PG_FUNCTION_ARGS)
{
    REACTION *rxn1 = PG_GETARG_REACTION_P (0);
    REACTION *rxn2 = PG_GETARG_REACTION_P (1);


    PG_RETURN_FLOAT8 (ob_tanimoto
                      ((const uint32 *) rxn1->fp,  (const uint32 *) rxn2->fp, OB_FPSIZE2));
}

Datum
reaction_similarity_products (PG_FUNCTION_ARGS)
{
    REACTION *rxn1 = PG_GETARG_REACTION_P (0);
    REACTION *rxn2 = PG_GETARG_REACTION_P (1);


    PG_RETURN_FLOAT8 (ob_tanimoto
                      ((const uint32 *) rxn1->fp+OB_FPSIZE2,  (const uint32 *) rxn2->fp+OB_FPSIZE2, OB_FPSIZE2));
}
