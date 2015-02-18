/************************************************************************
 * reaction_gist.c reaction GiST support functions
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
#include "access/skey.h"
#include "access/gist.h"
#include "c.h"
#include "libpq/md5.h"
#include "obwrapper.h"
#include "reaction.h"

#define GETENTRY(vec,pos) ((RXNFP *) DatumGetPointer((vec)->vector[(pos)].key))

Datum rxnfp_in (PG_FUNCTION_ARGS);
Datum rxnfp_out (PG_FUNCTION_ARGS);

PG_FUNCTION_INFO_V1 (rxnfp_in);
PG_FUNCTION_INFO_V1 (rxnfp_out);

/*
rxnfp input function. Dummy because it's formally needed but should never be called
*/
Datum
rxnfp_in (PG_FUNCTION_ARGS)
{
    elog (ERROR, "Not implemented");

    PG_RETURN_DATUM (PG_GETARG_DATUM (0));
}

/*
rxnfp output function. Dummy because it's formally needed but should never be called
*/
Datum
rxnfp_out (PG_FUNCTION_ARGS)
{
    elog (ERROR, "Not implemented");

    PG_RETURN_DATUM (PG_GETARG_DATUM (0));
}

PG_FUNCTION_INFO_V1 (rxnfp_consistent);
PG_FUNCTION_INFO_V1 (rxnfp_same);
PG_FUNCTION_INFO_V1 (rxnfp_compress);
PG_FUNCTION_INFO_V1 (rxnfp_decompress);
PG_FUNCTION_INFO_V1 (rxnfp_penalty);
PG_FUNCTION_INFO_V1 (rxnfp_picksplit);
PG_FUNCTION_INFO_V1 (rxnfp_union);

Datum rxnfp_consistent (PG_FUNCTION_ARGS);
Datum rxnfp_same (PG_FUNCTION_ARGS);
Datum rxnfp_compress (PG_FUNCTION_ARGS);
Datum rxnfp_decompress (PG_FUNCTION_ARGS);
Datum rxnfp_penalty (PG_FUNCTION_ARGS);
Datum rxnfp_picksplit (PG_FUNCTION_ARGS);
Datum rxnfp_union (PG_FUNCTION_ARGS);

inline static RXNFP *
new_rxnfp ()
{
    RXNFP *fp = (RXNFP *) palloc0 (sizeof(RXNFP));
    //memset(fp,0x0,sizeof(RXNFP));
    return fp;
}

inline static void
union_internal (RXNFP * result, RXNFP * element)
{
    register int i=2*OB_FPSIZE2;
    register uint32 *r = result->fp;
    register uint32 *e = element->fp;

    while (i--)
    {
        *r = *r | *e;
        r++;
        e++;
    }
}

inline static float
soergel_distance (RXNFP * fp1, RXNFP * fp2)
{
    return 1.0f - ob_tanimoto((const uint32 *) fp1->fp, (const uint32 *) fp2->fp, 2*OB_FPSIZE2);
}

/*
* Check, if a query molecule is contained in an index entry, i.e. if molecule is substructure or equal to predicate.
* If entry is a LEAF page an the strategy is RTSame, a equality check is performed, a substructure check in all
* other valid cases.
*/
Datum
rxnfp_consistent (PG_FUNCTION_ARGS)
{
    GISTENTRY *entry = (GISTENTRY *) PG_GETARG_POINTER (0);
    RXNFP *predicate = (RXNFP *) entry->key;
    REACTION *query = PG_GETARG_REACTION_P (1);
    StrategyNumber strategy = (StrategyNumber) PG_GETARG_UINT16 (2);
    bool *recheck = (bool *) PG_GETARG_POINTER(4);
    register int i=2*OB_FPSIZE2;
    register uint32 *p = predicate->fp;
    register uint32 *q = query->fp;

    *recheck = true; //Always lossy

    if(strategy == RTSameStrategyNumber && GIST_LEAF(entry))
    {

        while(i--)
        {
            if (*p != *q)
            {
                PG_RETURN_BOOL (false);
            }
            p++;
            q++;
        }

    }
    else
    {

        while(i--)
        {
            if ((*p & *q) != *q)
            {
                PG_RETURN_BOOL (false);
            }
            p++;
            q++;
        }
    }

    PG_RETURN_BOOL (true);
}

/*
* Check, if a query index entry is equal to another index entry.
* I haven't seen this being called yet.
*/
Datum
rxnfp_same (PG_FUNCTION_ARGS)
{
    RXNFP *query = (RXNFP *) PG_GETARG_POINTER (0);
    RXNFP *entry = (RXNFP *) PG_GETARG_POINTER (1);
    bool  *result = (bool *) PG_GETARG_POINTER(2);
    register int i=2*OB_FPSIZE2;
    register uint32 *q = query->fp;
    register uint32 *e = entry->fp;

    while (i--)
    {
        if (*e != *q)
        {
             *result = false;
            PG_RETURN_BOOL (result);
        }
        e++;
        q++;
    }

    *result=true;
    PG_RETURN_BOOL (result);
}

/*
* Compress a molecule of type molecule into an index entry of type rxnfp.
* This is the conversion that makes index entries out of table entries.
*/
Datum
rxnfp_compress (PG_FUNCTION_ARGS)
{
    GISTENTRY *entry = (GISTENTRY *) PG_GETARG_POINTER (0);
    GISTENTRY *retval = entry;

    if (entry->leafkey)
    {
        REACTION *rxn =
            (REACTION *) DatumGetPointer (PG_DETOAST_DATUM (entry->key));
        RXNFP *fp = new_rxnfp ();

        memcpy (fp->fp, rxn->fp, 2*OB_FPSIZE2*sizeof(uint32));

        retval = (GISTENTRY *) palloc (sizeof (GISTENTRY));

        gistentryinit (*retval, PointerGetDatum (fp),
                       entry->rel, entry->page, entry->offset, FALSE);
    }
    PG_RETURN_POINTER (retval);
}

/*
* Decompress an index entry of type rxnfp into a molecule.
* This is impossible here as the molecule cannot be restored from it's fingerprint, but
* the function has to be there for formal reasons.
*/
Datum
rxnfp_decompress (PG_FUNCTION_ARGS)
{
    PG_RETURN_DATUM (PG_GETARG_DATUM (0));
}

/*
* Calculate a penalty for inserting an index entry into an existing index page.
* GiST uses this to determine, if a new entry is inserted into the left or right branch
* of the index tree. The page with the lower penalty is chosen.
* Here, the penalty is the Soergel distance between the existing fingerprint in the page
* and the one to be newly inserted.
*/
Datum
rxnfp_penalty (PG_FUNCTION_ARGS)
{
    GISTENTRY *indexentry = (GISTENTRY *) PG_GETARG_POINTER (0);
    GISTENTRY *newentry = (GISTENTRY *) PG_GETARG_POINTER (1);
    float *penalty = (float *) PG_GETARG_POINTER (2);

    *penalty =
        soergel_distance ((RXNFP *) DatumGetPointer (indexentry->key),
                          (RXNFP *) DatumGetPointer (newentry->key));

    PG_RETURN_POINTER (penalty);
}


// Simple Split-By-Half
/*Datum
rxnfp_picksplit (PG_FUNCTION_ARGS)
{
  GistEntryVector *entryvec = (GistEntryVector *) PG_GETARG_POINTER (0);
  GIST_SPLITVEC *v = (GIST_SPLITVEC *) PG_GETARG_POINTER (1);
  int32 len = entryvec->n;
  int32 lenhalf = len/2;
  OffsetNumber i;
  RXNFP *datum_l, *datum_r;

  v->spl_nright = v->spl_nleft = 0;

  //printf ("picksplit\n");

  v->spl_left = (OffsetNumber *) palloc (len * sizeof (OffsetNumber));
  v->spl_right = (OffsetNumber *) palloc (len * sizeof (OffsetNumber));

  datum_l = new_rxnfp();
  datum_r = new_rxnfp();

  //memset (datum_l, 0, RFPSIZE2);
  //memset (datum_r, 0, RFPSIZE2);

  //printf("1\n");

  for (i = FirstOffsetNumber; i < len; i = OffsetNumberNext (i))
    {
      if (i < lenhalf)
		 		  		 		   {
		 		  		 		     //printf("4\n");
		 		  		 		     union_internal (datum_l, GETENTRY (entryvec, i));

		 		  		 		     //printf("5\n");
		 		  		 		     v->spl_left[v->spl_nleft++] = i;
		 		  		 		   }
      else
		 		  		 		   {

		 		  		 		     union_internal (datum_r, GETENTRY (entryvec, i));

		 		  		 		     //printf("7\n");
		 		  		 		     v->spl_right[v->spl_nright++] = i;
		 		  		 		   }
    }

  //printf("8\n");
  v->spl_ldatum = (Datum) datum_l;
  //printf("9\n");
  v->spl_rdatum = (Datum) datum_r;

  //printf("10\n");
  PG_RETURN_POINTER (v);

}*/

/*
* Decide which index entries are moved to which page, when a page split becomes neccessary, i.e.
* the current page would overflow with the new entry to be inserted.
* This function uses a Guttmann like Quadratic-Split algorithm, based on Soergel distance.
* The idea behind this is to minimize the chance of having to search the left and right
* tree branch, by first splitting along the two least similar entries, which become the initial (seed)
* entries of the new pages and then assigning the remaining ones by higher similarity to either the left
* or right initial entry.
*/
Datum
rxnfp_picksplit (PG_FUNCTION_ARGS)
{
    GistEntryVector *entryvec = (GistEntryVector *) PG_GETARG_POINTER (0);
    GIST_SPLITVEC *v = (GIST_SPLITVEC *) PG_GETARG_POINTER (1);
    int32 len = entryvec->n;
    RXNFP *entry, *entry_l, *entry_r;
    OffsetNumber i, j;
    RXNFP *datum_l, *datum_r;
    int seed_l = 0, seed_r = 0;
    float new_dist, dist = -1.0f, delta_r, delta_l;

    v->spl_nright = v->spl_nleft = 0;

    v->spl_left = (OffsetNumber *) palloc (len * sizeof (OffsetNumber));
    v->spl_right = (OffsetNumber *) palloc (len * sizeof (OffsetNumber));

    datum_l = new_rxnfp ();
    datum_r = new_rxnfp ();

// Initially, find the two most distant elements

    for (i = FirstOffsetNumber; i < len; i = OffsetNumberNext (i))
    {
        entry = GETENTRY (entryvec, i);

        for (j = OffsetNumberNext (i); j < len; j = OffsetNumberNext (j))
        {
            new_dist = soergel_distance (entry, GETENTRY (entryvec, j));
            //new_dist = new_dist * abs(entry->hvycount - (GETENTRY (entryvec, j))->hvycount);

            if (new_dist > dist)
            {
                dist = new_dist;
                seed_l = i;
                seed_r = j;
            }
        }
    }

    if (seed_l == 0 || seed_r == 0)
    {
        seed_l = 1;
        seed_r = 2;
    }

// Assign the seeds as left and right first elements to the split vector.
// Initialize the union pages

    entry_l = GETENTRY (entryvec, seed_l);
    entry_r = GETENTRY (entryvec, seed_r);

    v->spl_left[v->spl_nleft++] = seed_l;
    memcpy (datum_l, entry_l, sizeof(RXNFP));

    v->spl_right[v->spl_nright++] = seed_r;
    memcpy (datum_r, entry_r, sizeof(RXNFP));

    // Assign other elements to seed element with smaller distance

    for (i = FirstOffsetNumber; i < len; i = OffsetNumberNext (i))
    {
        if (i != seed_l && i != seed_r)
        {
            delta_l = soergel_distance (entry_l, GETENTRY (entryvec, i));
            delta_r = soergel_distance (entry_r, GETENTRY (entryvec, i));

            if (delta_l < delta_r)
            {
                union_internal (datum_l, GETENTRY (entryvec, i));
                v->spl_left[v->spl_nleft++] = i;
            }
            else
            {
                union_internal (datum_r, GETENTRY (entryvec, i));
                v->spl_right[v->spl_nright++] = i;
            }
        }
    }

    v->spl_ldatum = (Datum) datum_l;
    v->spl_rdatum = (Datum) datum_r;

    PG_RETURN_POINTER (v);
}

/*
* Given a set of entries, this function generates a new predicate that is true for all the entries.
* Here, this is achieved by simply OR'ing the bits of all index entries together.
*/
Datum
rxnfp_union (PG_FUNCTION_ARGS)
{
    GistEntryVector *entryvec = (GistEntryVector *) PG_GETARG_POINTER (0);
    int32 i;
    int32 len = entryvec->n;
    int *size = (int *) PG_GETARG_POINTER (1);
    RXNFP *result = new_rxnfp ();

    for (i = 0; i < len; i++)
    {
        union_internal (result, GETENTRY (entryvec, i));
    }

    *size = 2*OB_FPSIZE2*sizeof(uint32);

    PG_RETURN_POINTER (result);
}


