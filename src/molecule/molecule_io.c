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
#ifdef BUILD_WITH_INDIGO
#include "mingw/indigo.h"
#endif
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

///*static uint8 *
//createPaddedCopyWithLength(uint8 *b, uint32 *l)
//{
//	uint8	   *ret;
//	uint32		q;
//	uint32		len,
//				newLen448;
//	uint32		len_high,
//				len_low;		// 64-bit value split into 32-bit sections
//
//	len = ((b == NULL) ? 0 : *l);
//	newLen448 = len + 64 - (len % 64) - 8;
//	if (newLen448 <= len)
//		newLen448 += 64;
//
//	*l = newLen448 + 8;
//	if ((ret = (uint8 *) malloc(sizeof(uint8) * *l)) == NULL)
//		return NULL;
//
//	if (b != NULL)
//		memcpy(ret, b, sizeof(uint8) * len);
//
//	// pad
//	ret[len] = 0x80;
//	for (q = len + 1; q < newLen448; q++)
//		ret[q] = 0x00;
//
//	// append length as a 64 bit bitcount
//	len_low = len;
//	// split into two 32-bit values
//	// we only look at the bottom 32-bits
//	len_high = len >> 29;
//	len_low <<= 3;
//	q = newLen448;
//	ret[q++] = (len_low & 0xff);
//	len_low >>= 8;
//	ret[q++] = (len_low & 0xff);
//	len_low >>= 8;
//	ret[q++] = (len_low & 0xff);
//	len_low >>= 8;
//	ret[q++] = (len_low & 0xff);
//	ret[q++] = (len_high & 0xff);
//	len_high >>= 8;
//	ret[q++] = (len_high & 0xff);
//	len_high >>= 8;
//	ret[q++] = (len_high & 0xff);
//	len_high >>= 8;
//	ret[q] = (len_high & 0xff);
//
//	return ret;
//}
//
//#define F(x, y, z) (((x) & (y)) | (~(x) & (z)))
//#define G(x, y, z) (((x) & (z)) | ((y) & ~(z)))
//#define H(x, y, z) ((x) ^ (y) ^ (z))
//#define I(x, y, z) ((y) ^ ((x) | ~(z)))
//#define ROT_LEFT(x, n) (((x) << (n)) | ((x) >> (32 - (n))))
//
//static void
//doTheRounds(uint32 X[16], uint32 state[4])
//{
//	uint32		a,
//				b,
//				c,
//				d;
//
//	a = state[0];
//	b = state[1];
//	c = state[2];
//	d = state[3];
//
//	// round 1
//	a = b + ROT_LEFT((a + F(b, c, d) + X[0] + 0xd76aa478), 7);	/* 1 */
//	d = a + ROT_LEFT((d + F(a, b, c) + X[1] + 0xe8c7b756), 12); /* 2 */
//	c = d + ROT_LEFT((c + F(d, a, b) + X[2] + 0x242070db), 17); /* 3 */
//	b = c + ROT_LEFT((b + F(c, d, a) + X[3] + 0xc1bdceee), 22); /* 4 */
//	a = b + ROT_LEFT((a + F(b, c, d) + X[4] + 0xf57c0faf), 7);	/* 5 */
//	d = a + ROT_LEFT((d + F(a, b, c) + X[5] + 0x4787c62a), 12); /* 6 */
//	c = d + ROT_LEFT((c + F(d, a, b) + X[6] + 0xa8304613), 17); /* 7 */
//	b = c + ROT_LEFT((b + F(c, d, a) + X[7] + 0xfd469501), 22); /* 8 */
//	a = b + ROT_LEFT((a + F(b, c, d) + X[8] + 0x698098d8), 7);	/* 9 */
//	d = a + ROT_LEFT((d + F(a, b, c) + X[9] + 0x8b44f7af), 12); /* 10 */
//	c = d + ROT_LEFT((c + F(d, a, b) + X[10] + 0xffff5bb1), 17);		/* 11 */
//	b = c + ROT_LEFT((b + F(c, d, a) + X[11] + 0x895cd7be), 22);		/* 12 */
//	a = b + ROT_LEFT((a + F(b, c, d) + X[12] + 0x6b901122), 7); /* 13 */
//	d = a + ROT_LEFT((d + F(a, b, c) + X[13] + 0xfd987193), 12);		/* 14 */
//	c = d + ROT_LEFT((c + F(d, a, b) + X[14] + 0xa679438e), 17);		/* 15 */
//	b = c + ROT_LEFT((b + F(c, d, a) + X[15] + 0x49b40821), 22);		/* 16 */
//
//	/* round 2 */
//	a = b + ROT_LEFT((a + G(b, c, d) + X[1] + 0xf61e2562), 5);	/* 17 */
//	d = a + ROT_LEFT((d + G(a, b, c) + X[6] + 0xc040b340), 9);	/* 18 */
//	c = d + ROT_LEFT((c + G(d, a, b) + X[11] + 0x265e5a51), 14);		/* 19 */
//	b = c + ROT_LEFT((b + G(c, d, a) + X[0] + 0xe9b6c7aa), 20); /* 20 */
//	a = b + ROT_LEFT((a + G(b, c, d) + X[5] + 0xd62f105d), 5);	/* 21 */
//	d = a + ROT_LEFT((d + G(a, b, c) + X[10] + 0x02441453), 9); /* 22 */
//	c = d + ROT_LEFT((c + G(d, a, b) + X[15] + 0xd8a1e681), 14);		/* 23 */
//	b = c + ROT_LEFT((b + G(c, d, a) + X[4] + 0xe7d3fbc8), 20); /* 24 */
//	a = b + ROT_LEFT((a + G(b, c, d) + X[9] + 0x21e1cde6), 5);	/* 25 */
//	d = a + ROT_LEFT((d + G(a, b, c) + X[14] + 0xc33707d6), 9); /* 26 */
//	c = d + ROT_LEFT((c + G(d, a, b) + X[3] + 0xf4d50d87), 14); /* 27 */
//	b = c + ROT_LEFT((b + G(c, d, a) + X[8] + 0x455a14ed), 20); /* 28 */
//	a = b + ROT_LEFT((a + G(b, c, d) + X[13] + 0xa9e3e905), 5); /* 29 */
//	d = a + ROT_LEFT((d + G(a, b, c) + X[2] + 0xfcefa3f8), 9);	/* 30 */
//	c = d + ROT_LEFT((c + G(d, a, b) + X[7] + 0x676f02d9), 14); /* 31 */
//	b = c + ROT_LEFT((b + G(c, d, a) + X[12] + 0x8d2a4c8a), 20);		/* 32 */
//
//	/* round 3 */
//	a = b + ROT_LEFT((a + H(b, c, d) + X[5] + 0xfffa3942), 4);	/* 33 */
//	d = a + ROT_LEFT((d + H(a, b, c) + X[8] + 0x8771f681), 11); /* 34 */
//	c = d + ROT_LEFT((c + H(d, a, b) + X[11] + 0x6d9d6122), 16);		/* 35 */
//	b = c + ROT_LEFT((b + H(c, d, a) + X[14] + 0xfde5380c), 23);		/* 36 */
//	a = b + ROT_LEFT((a + H(b, c, d) + X[1] + 0xa4beea44), 4);	/* 37 */
//	d = a + ROT_LEFT((d + H(a, b, c) + X[4] + 0x4bdecfa9), 11); /* 38 */
//	c = d + ROT_LEFT((c + H(d, a, b) + X[7] + 0xf6bb4b60), 16); /* 39 */
//	b = c + ROT_LEFT((b + H(c, d, a) + X[10] + 0xbebfbc70), 23);		/* 40 */
//	a = b + ROT_LEFT((a + H(b, c, d) + X[13] + 0x289b7ec6), 4); /* 41 */
//	d = a + ROT_LEFT((d + H(a, b, c) + X[0] + 0xeaa127fa), 11); /* 42 */
//	c = d + ROT_LEFT((c + H(d, a, b) + X[3] + 0xd4ef3085), 16); /* 43 */
//	b = c + ROT_LEFT((b + H(c, d, a) + X[6] + 0x04881d05), 23); /* 44 */
//	a = b + ROT_LEFT((a + H(b, c, d) + X[9] + 0xd9d4d039), 4);	/* 45 */
//	d = a + ROT_LEFT((d + H(a, b, c) + X[12] + 0xe6db99e5), 11);		/* 46 */
//	c = d + ROT_LEFT((c + H(d, a, b) + X[15] + 0x1fa27cf8), 16);		/* 47 */
//	b = c + ROT_LEFT((b + H(c, d, a) + X[2] + 0xc4ac5665), 23); /* 48 */
//
//	/* round 4 */
//	a = b + ROT_LEFT((a + I(b, c, d) + X[0] + 0xf4292244), 6);	/* 49 */
//	d = a + ROT_LEFT((d + I(a, b, c) + X[7] + 0x432aff97), 10); /* 50 */
//	c = d + ROT_LEFT((c + I(d, a, b) + X[14] + 0xab9423a7), 15);		/* 51 */
//	b = c + ROT_LEFT((b + I(c, d, a) + X[5] + 0xfc93a039), 21); /* 52 */
//	a = b + ROT_LEFT((a + I(b, c, d) + X[12] + 0x655b59c3), 6); /* 53 */
//	d = a + ROT_LEFT((d + I(a, b, c) + X[3] + 0x8f0ccc92), 10); /* 54 */
//	c = d + ROT_LEFT((c + I(d, a, b) + X[10] + 0xffeff47d), 15);		/* 55 */
//	b = c + ROT_LEFT((b + I(c, d, a) + X[1] + 0x85845dd1), 21); /* 56 */
//	a = b + ROT_LEFT((a + I(b, c, d) + X[8] + 0x6fa87e4f), 6);	/* 57 */
//	d = a + ROT_LEFT((d + I(a, b, c) + X[15] + 0xfe2ce6e0), 10);		/* 58 */
//	c = d + ROT_LEFT((c + I(d, a, b) + X[6] + 0xa3014314), 15); /* 59 */
//	b = c + ROT_LEFT((b + I(c, d, a) + X[13] + 0x4e0811a1), 21);		/* 60 */
//	a = b + ROT_LEFT((a + I(b, c, d) + X[4] + 0xf7537e82), 6);	/* 61 */
//	d = a + ROT_LEFT((d + I(a, b, c) + X[11] + 0xbd3af235), 10);		/* 62 */
//	c = d + ROT_LEFT((c + I(d, a, b) + X[2] + 0x2ad7d2bb), 15); /* 63 */
//	b = c + ROT_LEFT((b + I(c, d, a) + X[9] + 0xeb86d391), 21); /* 64 */
//
//	state[0] += a;
//	state[1] += b;
//	state[2] += c;
//	state[3] += d;
//}
//
//static int
//calculateDigestFromBuffer(uint8 *b, uint32 len, uint8 sum[16])
//{
//	register uint32 i,
//				j,
//				k,
//				newI;
//	uint32		l;
//	uint8	   *input;
//	register uint32 *wbp;
//	uint32		workBuff[16],
//				state[4];
//
//	l = len;
//
//	state[0] = 0x67452301;
//	state[1] = 0xEFCDAB89;
//	state[2] = 0x98BADCFE;
//	state[3] = 0x10325476;
//
//	if ((input = createPaddedCopyWithLength(b, &l)) == NULL)
//		return 0;
//
//	for (i = 0;;)
//	{
//		if ((newI = i + 16 * 4) > l)
//			break;
//		k = i + 3;
//		for (j = 0; j < 16; j++)
//		{
//			wbp = (workBuff + j);
//			*wbp = input[k--];
//			*wbp <<= 8;
//			*wbp |= input[k--];
//			*wbp <<= 8;
//			*wbp |= input[k--];
//			*wbp <<= 8;
//			*wbp |= input[k];
//			k += 7;
//		}
//		doTheRounds(workBuff, state);
//		i = newI;
//	}
//	free(input);
//
//	j = 0;
//	for (i = 0; i < 4; i++)
//	{
//		k = state[i];
//		sum[j++] = (k & 0xff);
//		k >>= 8;
//		sum[j++] = (k & 0xff);
//		k >>= 8;
//		sum[j++] = (k & 0xff);
//		k >>= 8;
//		sum[j++] = (k & 0xff);
//	}
//	return 1;
//}*/

/*inline static char *
get_bzhash ()
{
  if (inited == false)
    {
      char *inchi;
      char *molfile;

      molfile = ob_smiles_to_mol (BZSMI);

      inchi = ob_mol_to_inchi (molfile);

      //pg_md5_hash (inchi, strlen (inchi) + 1, bzhash);
      calculateDigestFromBuffer(inchi, strlen (inchi) + 1, bzhash);

      bzhash[16]='\0';

      free (molfile);
      free (inchi);

      inited = true;
    }
  return bzhash;
}*/

/*inline static void merge_fps(unsigned int *fp2, unsigned int *fp3)
{
        fp2[0] |= fp3[0];
        fp2[8] |= fp3[1];
        fp2[16] |= fp3[2];
        fp2[24] |= fp3[3];
}  */

#ifdef BUILD_WITH_INDIGO
MOLECULE *
new_molecule (char *smiles, char *original_data)
{
    unsigned int sizeo;
    unsigned int sizesmi;
    size_t totalsize;
    MOLECULE *result;
    char *inchikey = NULL;
    unsigned char *ancillarydata;
    int sizeanc;
    int molhandle, fphandle, fpsize;
    char *fpbuf, *offset_fp;
    COMPRESSED_DATA *compressed_data;
    char checkresult[255];

    molhandle = indigoLoadMoleculeFromString(smiles);

    if(molhandle < 0)
    {
        indigoFree(molhandle);
        elog (ERROR, "1. Molecule generation failed! Offender was :\n %s\nError was: %s",original_data, indigoGetLastError());
    }

    strncpy(checkresult,indigoCheckBadValence(molhandle),255*sizeof(char));

    if(strlen(checkresult) != 0)
    {

        indigoFree(molhandle);
        elog (ERROR, "2. Molecule generation failed! Offender was :\n %s\nError was: %s",original_data, checkresult);
    }

    strncpy(checkresult,indigoCheckAmbiguousH(molhandle),255*sizeof(char));

    if (strlen(checkresult) != 0)
    {
        indigoAromatize(molhandle);
        indigoDearomatize(molhandle);

        strncpy(checkresult,indigoCheckAmbiguousH(molhandle),255*sizeof(char));

        if(strlen(checkresult) != 0)
        {
            indigoFree(molhandle);
            elog (ERROR, "3. Molecule generation failed! Offender was :\n %s\nError was: %s",original_data, checkresult);
        }
    }

    //indigoSetOptionInt("layout-max-iterations",1);

    //indigoLayout(molhandle);

    //newmolfile = strdup(indigoMolfile(molhandle));

    sizeo = strlen (original_data)+sizeof(char);

    sizesmi = strlen (smiles)+sizeof(char);

    compressed_data = compress_data(original_data, sizeo);

    indigoAromatize(molhandle);

    indigoSerialize(molhandle,&ancillarydata, &sizeanc);

    totalsize = CALCDATASZ (compressed_data->compressed_size, sizesmi, sizeanc);

    result = (MOLECULE *) palloc0 (totalsize);

    //memset (result, 0x0, totalsize);

    if (strchr (smiles, '.') != NULL)
        result->disconnected = true;

    result->sizeo = sizeo;
    result->compressed_sizeo=compressed_data->compressed_size;
    result->sizesmi = sizesmi;
    result->sizeanc = sizeanc;

    strncpy (SMIPTR(result), smiles, sizesmi);

    //strncpy (MFPTR(result), molfile, sizemf);
    memcpy(CIPTR(result),compressed_data->compressed_data, compressed_data->compressed_size);

    memcpy(ANCPTR(result), ancillarydata, sizeanc);

    fphandle = indigoFingerprint(molhandle,"full");

    indigoToBuffer (fphandle, &fpbuf, &fpsize);

    assert(fpsize*sizeof(char) <= sizeof(MOLFP));

    offset_fp = result->fp.bytes;

    offset_fp+IN_FPSTART;

    memcpy(offset_fp,fpbuf,fpsize*sizeof(char));

    //result->fp.popcount = ob_popcount(result->fp.dwords,IN_FPSIZE);

    indigoFree(fphandle);

    indigoFree(molhandle);

    inchikey = ob_smiles_to_inchikey (smiles);

    if(inchikey == NULL)
    {
        goto inchikey_fail;
    }
    else if (strlen(inchikey) != INCHIKEYSZ)
    {
        elog (ERROR, "%s INCHIKEYSIZE %d != %d\n",inchikey, INCHIKEYSZ, strlen(inchikey));
        free(inchikey);
inchikey_fail:
        elog (ERROR, "4. Molecule generation failed! Offender was :\n %s",original_data);
    }

    memcpy(result->inchikey, inchikey, INCHIKEYSZ);

    //free(newmolfile);
    free(inchikey);

    SET_VARSIZE (result,totalsize);

    return result;
}
#else
MOLECULE *
new_molecule (char *smiles, char *original_data, char *ancillarydata)
{
    unsigned int sizeo;
    unsigned int sizesmi;
    size_t totalsize;
    MOLECULE *result;
    char *inchikey = NULL;
    //char *ancillarydata = NULL;
    //uint32 *offset;
    char *aidata = NULL;
    uint32 ancsize = 0;
    COMPRESSED_DATA *compressed_data;

    if(ancillarydata==NULL)
    {
        ancillarydata = serializeOBMol(smiles, FORMAT_SMILES);

        if(ancillarydata == NULL)
        {
            elog (ERROR, "Molecule generation failed! Offender was :\n %s",original_data);
        }
    }

    ancsize = *(unsigned int*) ancillarydata;

    sizeo = strlen (original_data)+1;

    sizesmi = strlen (smiles)+1;

    compressed_data = compress_data(original_data, sizeo);

    totalsize = CALCDATASZ (compressed_data->compressed_size, sizesmi, ancsize);

    result = (MOLECULE *) palloc0 (totalsize);

    //memset (result, 0x0, totalsize);

    if (strchr (smiles, '.') != NULL)
        result->disconnected = true;

    result->sizeanc = ancsize;
    result->sizeo = sizeo;
    result->compressed_sizeo=compressed_data->compressed_size;
    result->sizesmi = sizesmi;

    strncpy (SMIPTR(result), smiles, sizesmi);

    //strncpy (MFPTR(result), molfile, sizemf);
    memcpy(CIPTR(result),compressed_data->compressed_data, compressed_data->compressed_size);

    aidata = (char*) &((unsigned int*)ancillarydata)[1];

    memcpy(ANCPTR(result), aidata, ancsize);

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
        elog (ERROR, "%s INCHIKEYSIZE %d != %d\n",inchikey, INCHIKEYSZ, strlen(inchikey));
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

    if(ancillarydata != NULL)
    {
        //printf("%d %d %d %d\n",ancsize,offset[0],offset[1],offset[2]);
        free(ancillarydata);
    }

    //memset(fp3,0x0,FPSIZE3*sizeof(unsigned int));
    //ob_fp2 (molfile, result->fp);

    //offset = result->fp+OFFSET;

    //ob_fp3 (molfile, offset);

    //result->popcount = ob_popcount((uint8 *)result->fp,FPSIZE*sizeof(uint32));

    //printf("popcount: %i\n",result->popcount);

    //merge_fps(result->fp2, fp3);

    /*if (strncmp (result->molhash, get_bzhash (), MOLHASHSZ) == 0)
      {
        result->nobz = false;
        result->isbz = true;
      }
    else if (ob_SSS_SMARTS_native (BZSMI, smiles) != 0)
      {
        result->nobz = false;
      }*/

    SET_VARSIZE (result,totalsize);

    return result;
}
#endif

/* static MOLECULE *make_molecule(char *raw_input, int size) {
  MOLECULE *result;
  char *input = NULL;
  char *molfile = NULL;
  char *smiles = NULL;
  char *endptr;
  bool freemolfile = false;
  bool freesmiles = false;
  unsigned int new_len;
  //unsigned int *efa_array = NULL;

  if(strstr (raw_input, "M  END") != NULL) {
      input = palloc (size+sizeof(char));
      memcpy (input, raw_input, size);
      endptr = strstr (input, "M  END") + strlen("M  END")*sizeof(char);
      *endptr = 0x0;
      new_len = strlen(input);
      pfree (input);
      input = palloc (new_len + 1);
      strncpy(input,raw_input,new_len);
      input[new_len] = 0x0;
      //strncat(input,raw_input,new_len);
  } else {
      //input = raw_input;
      input = palloc (size+1);
      memcpy (input, raw_input, size);
      input[size]=0x0;
  }

  if (strstr (input, "V2000") == NULL || strstr (input, "M  END") == NULL)
    {				//TODO: Bad Hack. If more fmts are needed, we need a dedicated converter function
      if (strstr (input, "V3000") != NULL && strstr (input, "M  END") != NULL) //V3000?
	{
	  molfile = ob_V3000_to_mol (input);

      if(molfile == NULL || !strlen(molfile) || strstr(molfile,"V3000")==NULL) {
          if(molfile!=NULL) free (molfile);
      elog (ERROR, "Molfile generation failed! Offender was :\n %s",input);
      }

	  smiles = ob_mol_to_smiles (input,0);

	  if(smiles == NULL || !strlen(smiles)) {
         elog (ERROR, "SMILES generation failed! Offender was :\n %s",input);
      }  else if (!strlen(smiles)) {
          elog (WARNING, "SMILES generation failed! Trying fallback...");
          free (smiles);
           smiles = ob_mol_to_canonical_smiles (input,0);
            if(smiles == NULL) {
                elog (ERROR, "SMILES generation finally failed! Offender was :\n %s",input);
            } else if (!strlen(smiles)) {
                free(smiles);
                elog (ERROR, "SMILES generation finally failed! Offender was :\n %s",input);
           }
           elog (WARNING, "Fallback OK");
      }

	  freemolfile = true;
	  freesmiles = true;
	}
      else if (strstr (input, "InChI=") != NULL) //InChI?
	{
	  molfile = ob_inchi_to_mol (input);

	  if(molfile == NULL || !strlen(molfile) || strstr(molfile,"V2000")==NULL) {
	      if(molfile!=NULL) free (molfile);
      elog (ERROR, "Molfile generation failed! Offender was :\n %s",input);
      }

	  smiles = ob_mol_to_smiles (molfile,0);

	  if(smiles == NULL || !strlen(smiles)) {
         elog (ERROR, "SMILES generation failed! Offender was :\n %s",input);
      }  else if (!strlen(smiles)) {
          elog (WARNING, "SMILES generation failed! Trying fallback...");
          free (smiles);
           smiles = ob_mol_to_canonical_smiles (input,0);
            if(smiles == NULL) {
                elog (ERROR, "SMILES generation finally failed! Offender was :\n %s",input);
            } else if (!strlen(smiles)) {
                free(smiles);
                elog (ERROR, "SMILES generation finally failed! Offender was :\n %s",input);
           }
            elog (WARNING, "Fallback OK");
      }

	  freemolfile = true;
	  freesmiles = true;
	}
      else //SMILES?
	{
	  molfile = ob_smiles_to_mol (input);

      if(molfile == NULL || !strlen(molfile)) {
          if(molfile!=NULL) free (molfile);
      elog (ERROR, "Molfile generation failed! Offender was :\n %s",input);
      }

	  smiles = input;

	  freemolfile = true;
	}
    }
  else //V2000
    {
      smiles = ob_mol_to_smiles (input,0);

      if(smiles == NULL || !strlen(smiles)) {
         elog (ERROR, "SMILES generation failed! Offender was :\n %s",input);
      }  else if (!strlen(smiles)) {
          elog (WARNING, "SMILES generation failed! Trying fallback...");
          free (smiles);
           smiles = ob_mol_to_canonical_smiles (input,0);
            if(smiles == NULL) {
                elog (ERROR, "SMILES generation finally failed! Offender was :\n %s",input);
            } else if (!strlen(smiles)) {
                free(smiles);
                elog (ERROR, "SMILES generation finally failed! Offender was :\n %s",input);
           }
            elog (WARNING, "Fallback OK");
      }

      molfile = input;
      freesmiles = true;
    }

  if (molfile == NULL || smiles == NULL)
    {
      if (smiles != NULL && freesmiles)
	free (smiles);
      if (molfile != NULL && freemolfile)
	free (molfile);
      elog (ERROR,
	    "Input is not a V2000/V3000 molfile or InChI or SMILES: %s",
	    input);
    }

  //efa_array = ob_efa_array(smiles);

  result = new_molecule (smiles, molfile);

  //if (efa_array != NULL) free(efa_array);

  if (smiles != NULL && freesmiles)
    free (smiles);
  if (molfile != NULL && freemolfile)
    free (molfile);

  if (input != NULL)
    pfree(input);

  return result;
}    */

static MOLECULE *make_molecule(const char *raw_input, int size)
{
    MOLECULE *result = NULL;
    char *input = NULL;
#ifndef BUILD_WITH_INDIGO
    char *ancillarydata = NULL;
#endif
    char *smiles = NULL;
    char *endptr = NULL;
    //bool freemolfile = false;
    bool freesmiles = false;
    unsigned int new_len = 0;
    unsigned int input_format = FORMAT_SMILES;

    if(strstr (raw_input, "M  END") != NULL)
    {
        input = palloc (size+sizeof(char));
        memcpy (input, raw_input, size);
        endptr = strstr (input, "M  END") + strlen("M  END")*sizeof(char);
        *endptr = 0x0;
        new_len = strlen(input);
        pfree (input);
        input = palloc (new_len + 1);
        strncpy(input,raw_input,new_len);
        input[new_len] = 0x0;
    }
    else
    {
        input = palloc (size+1);
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

#ifdef BUILD_WITH_INDIGO
    result = new_molecule (smiles, input);
#else
    ancillarydata = serializeOBMol(input, input_format);

    result = new_molecule (smiles, input, ancillarydata);
#endif

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
* Convert an old pgchem::tigress bytea molecule to a new one of type molecule
*/

/*PG_FUNCTION_INFO_V1 (pgchem_molecule_to_new_molecule);

Datum
pgchem_molecule_to_new_molecule (PG_FUNCTION_ARGS)
{
  bytea *old = PG_GETARG_BYTEA_P (0);
  char *molfile;
  char *endptr;
  MOLECULE *result;
  size_t totalsize;
  uint32 sizemf;
  uint32 sizesmi;
  char *smiles = NULL;
  char *inchikey = NULL;
  //uint32 *offset;
  char *ancillarydata = NULL;
  char *aidata = NULL;
  uint32 ancsize = 0;

  molfile = palloc (VARSIZE (old) - VARHDRSZ + 1);

  memcpy (molfile, VARDATA (old), VARSIZE (old) - VARHDRSZ);

  if (strstr (molfile, "V2000") == NULL || strstr (molfile, "M  END") == NULL)
    elog (ERROR, "Input is not a V2000 molfile: %s", molfile);

  endptr = strstr (molfile, "M  END") + 6;

  *endptr = 0x0;

  sizemf = strlen (molfile) + 1;

  smiles = ob_mol_to_smiles (molfile,0);

 if(smiles == NULL || !strlen(smiles)) {
         elog (ERROR, "SMILES generation failed! Offender was :\n %s",molfile);
      }  else if (!strlen(smiles)) {
          elog (WARNING, "SMILES generation failed! Trying fallback...");
          free (smiles);
           smiles = ob_mol_to_canonical_smiles (molfile,0);
            if(smiles == NULL) {
                elog (ERROR, "SMILES generation finally failed! Offender was :\n %s",molfile);
            } else if (!strlen(smiles)) {
                free(smiles);
                free(ancillarydata);
                elog (ERROR, "SMILES generation finally failed! Offender was :\n %s",molfile);
           }
           elog (WARNING, "Fallback OK");
      }

  sizesmi = strlen (smiles) + 1;

  ancillarydata = ob_lyophilize_molecule(smiles);

  ancsize = *(unsigned int*) ancillarydata;

  if(ancillarydata == NULL)
      elog (ERROR, "Molecule generation failed! Offender was :\n %s",molfile);

  totalsize = CALCDATASZ (sizemf, sizesmi, ancsize);

  result = (MOLECULE *) palloc (totalsize);

  memset (result, 0x0, totalsize);

  result->sizemf = sizemf;
  result->sizesmi = sizesmi;
  //result->nobz = true;

  strncpy (SMIPTR(result), smiles, sizesmi);

  strncpy (MFPTR(result), molfile, sizemf);

  aidata = (char*) &((unsigned int*)ancillarydata)[1];

  memcpy(ANCPTR(result), aidata, ancsize);

  //free(efa_array);

  inchikey = ob_molfile_to_inchikey (molfile);

   if(inchikey == NULL) {
              goto inchikey_fail;
        }  else if (strlen(inchikey) != INCHIKEYSZ) {
              free(inchikey);
              inchikey_fail: elog (ERROR, "Molecule generation failed! Offender was :\n %s",molfile);
        }

  //pg_md5_hash (inchi, strlen (inchi) + 1, result->molhash);

  //calculateDigestFromBuffer(inchi, strlen (inchi), result->molhash);

  memcpy(result->inchikey, inchikey, INCHIKEYSZ);

  //calculateDigestFromBuffer(smiles, sizesmi, result->molhash);

  //result->molhash[16]='\0';

  free (inchikey);

  if (strchr (smiles, '.') != NULL)
    result->disconnected = true;

  //memset(fp3,0x0,FPSIZE3*sizeof(unsigned int));
  //ob_fp2 (molfile, result->fp);

  //offset = result->fp+OFFSET;

  //ob_fp3 (molfile, offset);

  ob_fp_bin(aidata, result->fp);

  if(ancillarydata != NULL) {
   //printf("%d %d %d %d\n",ancsize,offset[0],offset[1],offset[2]);
   free(ancillarydata);
  }

  //result->popcount = ob_popcount((uint8 *)result->fp,FPSIZE*sizeof(uint32));

  //merge_fps(result->fp2, fp3);

 /* if (strncmp (result->molhash, get_bzhash (), MOLHASHSZ) == 0)
    {
      result->nobz = false;
      result->isbz = true;
    }
  else if (ob_SSS_SMARTS_native (BZSMI, smiles) != 0)
    {
      result->nobz = false;
    }

  pfree (molfile);

  free (smiles);

  SET_VARSIZE (result,totalsize);

  PG_RETURN_MOLECULE_P (result);
}*/

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

/*PG_FUNCTION_INFO_V1 (pgchem_V3000_to_molecule);

Datum
pgchem_V3000_to_molecule (PG_FUNCTION_ARGS)
{
  MOLECULE *retval;
  text *arg_V3000;
  char *tmpV3000;
  char *molfile;
  char *smiles;
  int V3000_string_len;

  arg_V3000 = PG_GETARG_TEXT_P (0);

  V3000_string_len = VARSIZE (arg_V3000) - VARHDRSZ;

  tmpV3000 = (char *) palloc (V3000_string_len + 1);
  tmpV3000[0] = '\0';

  strncat (tmpV3000, VARDATA (arg_V3000), V3000_string_len);
  molfile = ob_V3000_to_mol (tmpV3000);
  smiles = ob_mol_to_canonical_smiles (molfile,1);

  retval = new_molecule(smiles,molfile);

  free(molfile);
  free(smiles);

  PG_RETURN_MOLECULE_P (retval);
}

PG_FUNCTION_INFO_V1 (pgchem_inchi_to_molecule);

Datum
pgchem_inchi_to_molecule (PG_FUNCTION_ARGS)
{
  MOLECULE *retval;
  text *arg_inchi;
  char *tmpinchi;
  char *molfile;
  char *smiles;
  int inchi_string_len;

  arg_inchi = PG_GETARG_TEXT_P (0);

  inchi_string_len = VARSIZE (arg_inchi) - VARHDRSZ;

  tmpinchi = (char *) palloc (inchi_string_len + 1);
  tmpinchi[0] = '\0';

  strncat (tmpinchi, VARDATA (arg_inchi), inchi_string_len);
  molfile = ob_inchi_to_mol (tmpinchi);
  smiles = ob_mol_to_canonical_smiles (molfile,1);

  retval = new_molecule(smiles,molfile);

  free(molfile);
  free(smiles);

  PG_RETURN_MOLECULE_P (retval);
} */

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

#ifdef BUILD_WITH_INDIGO
    retval = new_molecule (smiles, molfile);
#else
    retval = new_molecule (smiles, molfile, NULL);
#endif

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

#ifdef BUILD_WITH_INDIGO
    retval = new_molecule (smiles, molfile);
#else
    retval = new_molecule (smiles, molfile, NULL);
#endif

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

#ifdef BUILD_WITH_INDIGO
    retval = new_molecule (smiles, molfile);
#else
    retval = new_molecule (smiles, molfile, NULL);
#endif

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
