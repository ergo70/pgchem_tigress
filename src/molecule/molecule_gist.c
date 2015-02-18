/************************************************************************
 * molecule_gist.c molecule GiST support functions
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
#include "molecule.h"
#ifdef BUILD_WITH_INDIGO
#include "mingw/indigo.h"
#endif
#include <assert.h>

#define GETENTRY_MOLFP(vec,pos) ((MOLFP *) DatumGetPointer((vec)->vector[(pos)].key))

Datum molfp_in (PG_FUNCTION_ARGS);
Datum molfp_out (PG_FUNCTION_ARGS);

PG_FUNCTION_INFO_V1 (molfp_in);
PG_FUNCTION_INFO_V1 (molfp_out);

/*
molfp input function. Dummy because it's formally needed but should never be called
*/
Datum
molfp_in (PG_FUNCTION_ARGS)
{
    elog (ERROR, "Not implemented");

    PG_RETURN_DATUM (PG_GETARG_DATUM (0));
}

/*
molfp output function. Dummy because it's formally needed but should never be called
*/
Datum
molfp_out (PG_FUNCTION_ARGS)
{
    elog (ERROR, "Not implemented");

    PG_RETURN_DATUM (PG_GETARG_DATUM (0));
}

PG_FUNCTION_INFO_V1 (molfp_consistent);
PG_FUNCTION_INFO_V1 (molfp_same);
PG_FUNCTION_INFO_V1 (molfp_compress);
PG_FUNCTION_INFO_V1 (molfp_decompress);
PG_FUNCTION_INFO_V1 (molfp_penalty);
PG_FUNCTION_INFO_V1 (molfp_picksplit);
PG_FUNCTION_INFO_V1 (molfp_union);

Datum molfp_consistent (PG_FUNCTION_ARGS);
Datum molfp_same (PG_FUNCTION_ARGS);
Datum molfp_compress (PG_FUNCTION_ARGS);
Datum molfp_decompress (PG_FUNCTION_ARGS);
Datum molfp_penalty (PG_FUNCTION_ARGS);
Datum molfp_picksplit (PG_FUNCTION_ARGS);
Datum molfp_union (PG_FUNCTION_ARGS);

inline static MOLFP *
new_molfp ()
{
    MOLFP *fp = (MOLFP *) palloc0 (sizeof(MOLFP));
    //memset(fp,0x0,sizeof(MOLFP));
    return fp;
}

inline static void
union_internal (MOLFP * result, const MOLFP * element)
{
    /*#ifdef BUILD_WITH_INDIGO
        int i=IN_FPSIZE;
    #else
        int i=OB_FPSIZE;
    #endif
        const uint32 *e =   element->dwords;
        uint32 *r =  result->dwords;

        while (i--)
        {
            *r = (*r | *e);
            r++;
            e++;
        }*/

    result->v = result->v | element->v;
}

inline static float
soergel_distance (const MOLFP * fp1, const MOLFP * fp2)
{
#ifdef BUILD_WITH_INDIGO
    return 1.0f - ob_tanimoto
           (fp1->dwords+IN_SIMOFFSET,  fp2->dwords+IN_SIMOFFSET, IN_SIMFPSIZE);
#else
    return 1.0f - ob_tanimoto
           (fp1->dwords,  fp2->dwords, OB_FPSIZE2);
#endif
}

/*
* Check, if a query molecule is contained in an index entry, i.e. if molecule is substructure or equal to predicate.
* If entry is a LEAF page an the strategy is RTSame, a equality check is performed, a substructure check in all
* other valid cases.
*/

/*Datum
molfp_consistent (PG_FUNCTION_ARGS)
{
    GISTENTRY *entry = (GISTENTRY *) PG_GETARG_POINTER (0);
    MOLECULE *query = PG_GETARG_MOLECULE_P (1);
    StrategyNumber strategy = (StrategyNumber) PG_GETARG_UINT16 (2);
    // Oid subtype = PG_GETARG_OID(3);
    bool *recheck = (bool *) PG_GETARG_POINTER(4);
    MOLFP *predicate = (MOLFP *) DatumGetPointer(entry->key);
#ifdef BUILD_WITH_INDIGO
    int i = IN_SUBFPSIZE;
#else
int i = OB_FPSIZE;
#endif
    uint32 *p = predicate->dwords;
    uint32 *q = query->fp.dwords;

    *recheck = true; //Always lossy

    if(strategy == RTSameStrategyNumber && GIST_LEAF(entry))
    {

        while(i--)
        {
            if (*p != *q)
            {
                //elog(INFO,"p=%i q=%i",*p,*q);
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
                //elog(INFO,"p&q= %i p=%i q=%i",(*p & *q),*p,*q);
                PG_RETURN_BOOL (false);
            }
            p++;
            q++;
        }
    }

    //elog(INFO,"true");
    PG_RETURN_BOOL (true);
}*/

Datum
molfp_consistent (PG_FUNCTION_ARGS)
{
    GISTENTRY *entry = (GISTENTRY *) PG_GETARG_POINTER (0);
    MOLECULE *query = PG_GETARG_MOLECULE_P (1);
    StrategyNumber strategy = (StrategyNumber) PG_GETARG_UINT16 (2);
    /* Oid subtype = PG_GETARG_OID(3); */
    bool *recheck = (bool *) PG_GETARG_POINTER(4);
    MOLFP *predicate = (MOLFP *) DatumGetPointer(entry->key);
    MOLFP allzero;
    uint64 *zero;
#ifdef BUILD_WITH_INDIGO
    int i=IN_SUBFPSIZE/2;
#else
    int i= OB_FPSIZE/2;
#endif

    *recheck = true; //Always recheck except

    if(strategy > RTContainedByStrategyNumber) *recheck = false;

    if(strategy == RTSameStrategyNumber && GIST_LEAF(entry))
    {
        allzero.v = predicate->v ^ query->fp.v;

    }
    else
    {
        allzero.v = (predicate->v & query->fp.v) ^ query->fp.v;
    }

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
* Check, if a query index entry is equal to another index entry.
* I haven't seen this being called yet.
*/

Datum
molfp_same (PG_FUNCTION_ARGS)
{
    MOLFP *query = (MOLFP *) PG_GETARG_POINTER (0);
    MOLFP *entry = (MOLFP *) PG_GETARG_POINTER (1);
    bool  *result = (bool *) PG_GETARG_POINTER(2);
    MOLFP allzero;
    uint64 *zero;
#ifdef BUILD_WITH_INDIGO
    int i = IN_SUBFPSIZE/2;
#else
    int i = OB_FPSIZE/2;
#endif

    allzero.v = entry->v ^ query->v;

    zero = &allzero.qwords[0];

    while (i--)
    {
        if (*zero != 0x0ULL)
        {
            *result = false;
            PG_RETURN_BOOL (result);
        }
        zero++;
    }

    *result = true;
    PG_RETURN_BOOL (result);
}

/*Datum
molfp_same (PG_FUNCTION_ARGS)
{
    MOLFP *query = (MOLFP *) PG_GETARG_POINTER (0);
    MOLFP *entry = (MOLFP *) PG_GETARG_POINTER (1);
    bool  *result = (bool *) PG_GETARG_POINTER(2);
#ifdef BUILD_WITH_INDIGO
    int i = IN_SUBFPSIZE;
#else
    int i = OB_FPSIZE;
#endif
    uint32 *e = entry->dwords;
    uint32 *q = query->dwords;

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

    *result = true;
    PG_RETURN_BOOL (result);
}*/

/*
* Compress a molecule of type molecule into an index entry of type molfp.
* This is the conversion that makes index entries out of table entries.
*/
Datum
molfp_compress (PG_FUNCTION_ARGS)
{
    GISTENTRY *entry = (GISTENTRY *) PG_GETARG_POINTER (0);
    GISTENTRY *retval;

    if (entry->leafkey)
    {
        MOLECULE *mol =
            (MOLECULE *) DatumGetPointer (PG_DETOAST_DATUM (entry->key));

        MOLFP *fp = new_molfp();

        memcpy (fp->bytes, mol->fp.bytes, sizeof(MOLFP));

        retval = (GISTENTRY *) palloc (sizeof (GISTENTRY));

        gistentryinit (*retval, PointerGetDatum (fp),
                       entry->rel, entry->page, entry->offset, FALSE);
    }
    else
    {
        retval = entry;
    }

    PG_RETURN_POINTER (retval);
}
/*
* Decompress an index entry of type molfp into a molecule.
* This is impossible here as the molecule cannot be restored from it's fingerprint, but
* the function has to be there for formal reasons.
*/
Datum
molfp_decompress (PG_FUNCTION_ARGS)
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
molfp_penalty (PG_FUNCTION_ARGS)
{
    GISTENTRY *indexentry = (GISTENTRY *) PG_GETARG_POINTER (0);
    GISTENTRY *newentry = (GISTENTRY *) PG_GETARG_POINTER (1);
    float *penalty = (float *) PG_GETARG_POINTER (2);

    *penalty =
        soergel_distance ((MOLFP *) DatumGetPointer (indexentry->key),
                          (MOLFP *) DatumGetPointer (newentry->key));

    PG_RETURN_POINTER (penalty);
}


// Simple Split-By-Half
/*Datum
molfp_picksplit (PG_FUNCTION_ARGS)
{
  GistEntryVector *entryvec = (GistEntryVector *) PG_GETARG_POINTER (0);
  GIST_SPLITVEC *v = (GIST_SPLITVEC *) PG_GETARG_POINTER (1);
  int32 len = entryvec->n;
  int32 lenhalf = len/2;
  OffsetNumber i;
  MOLFP *datum_l, *datum_r;

  v->spl_nright = v->spl_nleft = 0;

  //printf ("picksplit\n");

  v->spl_left = (OffsetNumber *) palloc (len * sizeof (OffsetNumber));
  v->spl_right = (OffsetNumber *) palloc (len * sizeof (OffsetNumber));

  datum_l = new_molfp();
  datum_r = new_molfp();

  //memset (datum_l, 0, FPSIZE2);
  //memset (datum_r, 0, FPSIZE2);

  //printf("1\n");

  for (i = FirstOffsetNumber; i < len; i = OffsetNumberNext (i))
    {
      if (i < lenhalf)
		 		  		 		   {
		 		  		 		     //printf("4\n");
		 		  		 		     union_internal (datum_l, GETENTRY_MOLFP (entryvec, i));

		 		  		 		     //printf("5\n");
		 		  		 		     v->spl_left[v->spl_nleft++] = i;
		 		  		 		   }
      else
		 		  		 		   {

		 		  		 		     union_internal (datum_r, GETENTRY_MOLFP (entryvec, i));

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
molfp_picksplit (PG_FUNCTION_ARGS)
{
    GistEntryVector *entryvec = (GistEntryVector *) PG_GETARG_POINTER (0);
    GIST_SPLITVEC *v = (GIST_SPLITVEC *) PG_GETARG_POINTER (1);
    int32 len = entryvec->n;
    MOLFP *entry, *entry_l, *entry_r;
    OffsetNumber i, j;
    MOLFP *datum_l, *datum_r;
    int seed_l = 0, seed_r = 0;
    float new_dist, dist = -1.0f, delta_r, delta_l;

    v->spl_nright = v->spl_nleft = 0;

    v->spl_left = (OffsetNumber *) palloc (len * sizeof (OffsetNumber));
    v->spl_right = (OffsetNumber *) palloc (len * sizeof (OffsetNumber));

    datum_l = new_molfp ();
    datum_r = new_molfp ();

// Initially, find the two most distant elements

    for (i = FirstOffsetNumber; i < len; i = OffsetNumberNext (i))
    {
        entry = GETENTRY_MOLFP (entryvec, i);

        for (j = OffsetNumberNext (i); j < len; j = OffsetNumberNext (j))
        {
            new_dist = soergel_distance (entry, GETENTRY_MOLFP (entryvec, j));
            //new_dist = new_dist * abs(entry->hvycount - (GETENTRY_MOLFP (entryvec, j))->hvycount);

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

    entry_l = GETENTRY_MOLFP (entryvec, seed_l);
    entry_r = GETENTRY_MOLFP (entryvec, seed_r);

    v->spl_left[v->spl_nleft++] = seed_l;
    memcpy (datum_l, entry_l, sizeof(MOLFP));

    v->spl_right[v->spl_nright++] = seed_r;
    memcpy (datum_r, entry_r, sizeof(MOLFP));

    // Assign other elements to seed element with smaller distance

    for (i = FirstOffsetNumber; i < len; i = OffsetNumberNext (i))
    {
        if (i != seed_l && i != seed_r)
        {
            delta_l = soergel_distance (entry_l, GETENTRY_MOLFP (entryvec, i));
            delta_r = soergel_distance (entry_r, GETENTRY_MOLFP (entryvec, i));

            if (delta_l < delta_r)
            {
                union_internal (datum_l, GETENTRY_MOLFP (entryvec, i));
                v->spl_left[v->spl_nleft++] = i;
            }
            else
            {
                union_internal (datum_r, GETENTRY_MOLFP (entryvec, i));
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
molfp_union (PG_FUNCTION_ARGS)
{
    GistEntryVector *entryvec = (GistEntryVector *) PG_GETARG_POINTER (0);
    int *size = (int *) PG_GETARG_POINTER (1);
    int range = entryvec->n, i;

    MOLFP *result = new_molfp ();
    memcpy(result, GETENTRY_MOLFP (entryvec, 0),sizeof(MOLFP));

    *size = sizeof(MOLFP);

    if(range==1)
    {
        PG_RETURN_POINTER (result);
    }

    for (i = 1; i < range; i++)
    {
        union_internal (result, GETENTRY_MOLFP (entryvec, i));
    }

    PG_RETURN_POINTER (result);
}


