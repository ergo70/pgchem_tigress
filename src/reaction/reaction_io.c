/************************************************************************
 * reaction_io.c reaction input/output support functions
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
#include <time.h>
#include "postgres.h"
#include "libpq/md5.h"
#include "fmgr.h"
#include "libpq/pqformat.h"	/* needed for send/recv functions */
#include "reaction.h"
#include "obwrapper.h"

#ifdef BUILD_WITH_INDIGO
    #error No reaction support when compiling with Indigo!
#endif

Datum reaction_in (PG_FUNCTION_ARGS);
Datum reaction_in_text (PG_FUNCTION_ARGS);
Datum reaction_in_varchar (PG_FUNCTION_ARGS);
Datum reaction_in_bytea (PG_FUNCTION_ARGS);
Datum reaction_out (PG_FUNCTION_ARGS);
Datum reaction_recv (PG_FUNCTION_ARGS);
Datum reaction_send (PG_FUNCTION_ARGS);
Datum pgchem_reaction_mol_strip_rxninfo (PG_FUNCTION_ARGS);

static REACTION *
new_reaction (MOLECULE *mols[], int numreactants, int numproducts, int mode)
{
    size_t totalsize;
    REACTION *result;
    char *offset;
    int datasize=0;
    int i,j;

    for(i=0; i<(numreactants+numproducts); i++)
    {
        datasize+=VARSIZE(mols[i])*sizeof(char);
    }

    totalsize = (RCALCDATASZ (datasize))*sizeof(char);

    result = (REACTION *) palloc0 (totalsize);

    //memset (result, 0x0, totalsize);

    result->datasize = datasize;
    result->num_reactants = numreactants;
    result->num_products = numproducts;
    result->mode = mode;

    offset = MOLARRAYPTR(result);

    for(i=0; i<numreactants; i++)
    {

        memcpy((MOLECULE*) offset, mols[i],VARSIZE(mols[i])*sizeof(char));

        for(j=0; j<OB_FPSIZE2; j++)
        {
            result->fp[j] |= mols[i]->fp.dwords[j];
        }

        offset += (VARSIZE(mols[i]))*sizeof(char);
    }

    for(i=numreactants; i<(numproducts+numreactants); i++)
    {

        memcpy((MOLECULE*)offset, mols[i],VARSIZE(mols[i])*sizeof(char));

        for(j=0; j<OB_FPSIZE2; j++)
        {
            result->fp[j+OB_FPSIZE2] |= mols[i]->fp.dwords[j];
        }

        offset += (VARSIZE(mols[i]))*sizeof(char);
    }

    SET_VARSIZE (result,totalsize);

    return result;
}

static REACTION *make_reaction(const char *raw_input, const int size)
{

    int i;
    int r,p;
    int mode=0;
    int maxbuflen;
    int mcount=0;
    char *tmpbuf;
    char *rxnfile;
    char *workptrA;
    char *workptrB;
    char *smiles = NULL;
//unsigned int *efa_array = NULL;
    REACTION *result;
    MOLECULE **mols;

    if (strchr(raw_input,'\n')==NULL) elog (ERROR, "Reaction generation failed! No line separators found. Offender was :\n %s",raw_input);

//printf("-1");

    if (strstr(raw_input,"$RXN")==NULL || strstr(raw_input,"$MOL")==NULL || strstr(raw_input,"M  END")==NULL) elog (ERROR, "Reaction generation failed! Invalid Reactionfile. Offender was :\n %s",raw_input);

    if (strstr(raw_input,"\r\n")!=NULL) mode=1;

    maxbuflen = (size+(mode==1 ? 3 : 2))*sizeof(char);

    rxnfile=(char*)palloc0(maxbuflen);

    //memset(rxnfile,0x0,maxbuflen);

    memcpy(rxnfile,raw_input,size*sizeof(char));

    workptrA = rxnfile;

    for(i=0; i<4; i++)
    {
        workptrA=strchr(workptrA,'\n')+sizeof(char);
    }

    sscanf(workptrA,"%3d%3d",&r,&p);

    if(r<0 || p<0) elog (ERROR, "Negative count in counts line of rxnfile! Offender was :\n %s",rxnfile);

    mcount = r+p;

    i = 0;

//printf("%s\n",workptrA);

    while((workptrA=strstr(workptrA,"$MOL")))
    {
        workptrA++;
        i++;
    }

    if(i != mcount) elog (ERROR, "Count mismatch in rxnfile! Offender was :\n %s\nDeclared: %d\nFound: %d\n",rxnfile,mcount,i);

//printf("%d %d\n",size,maxbuflen);

    if(mode==1) rxnfile[maxbuflen-3] = '\r';

    rxnfile[maxbuflen-2] = '\n';
    rxnfile[maxbuflen-1] = '\0';

    mols = palloc(mcount*sizeof(MOLECULE*));

    for(i=0; i<mcount; i++)
    {
        mols[i]=NULL;
    }

    workptrA = rxnfile;
    workptrB = workptrA;

    tmpbuf=(char*) palloc(maxbuflen);

    for(i=0; i<r; i++)
    {

        memset(tmpbuf,0x0,maxbuflen);

        if(mode==1)
        {
            workptrA=strstr(workptrA,"$MOL\r\n")+strlen("$MOL\r\n")*sizeof(char);
            workptrB = strstr(workptrA,"M  END\r\n")+strlen("M  END\r\n")*sizeof(char);
        }
        else
        {
            workptrA=strstr(workptrA,"$MOL\n")+strlen("$MOL\n")*sizeof(char);
            workptrB = strstr(workptrA,"M  END\n")+strlen("M  END\n")*sizeof(char);
        }

//printf("9");

        memcpy(tmpbuf,workptrA,(workptrB-workptrA)*sizeof(char));

//printf("10");

        workptrA=workptrB;

        smiles = ob_molfile_to_canonical_smiles (tmpbuf,0);

        if(smiles == NULL)
        {
            goto smiles_fail;
        }
        else if (!strlen(smiles))
        {
            free (smiles);
smiles_fail:
            elog (WARNING, "SMILES generation failed! Trying fallback...");

        }

        if(ob_is_nostruct(tmpbuf) != 0)
        {

            for(i=0; i<mcount; i++)
            {
                if(mols[i] != NULL) pfree(mols[i]);
            }

            free(smiles);
            pfree(mols);
            pfree(rxnfile);
            pfree(tmpbuf);

            elog (ERROR, "NoStructures are not allowed in reactions! Offender was :\n %s",raw_input);
        }

//efa_array = ob_efa_array(smiles);

        mols[i]=new_molecule(smiles,tmpbuf,NULL);

//free(efa_array);

        if (smiles != NULL) free(smiles);
    }

    for(i=r; i<(r+p); i++)
    {

        memset(tmpbuf,0x0,maxbuflen);

        if(mode==1)
        {
            workptrA=strstr(workptrA,"$MOL\r\n")+strlen("$MOL\r\n")*sizeof(char);
            workptrB = strstr(workptrA,"M  END\r\n")+strlen("M  END\r\n")*sizeof(char);
        }
        else
        {
            workptrA=strstr(workptrA,"$MOL\n")+strlen("$MOL\n")*sizeof(char);
            workptrB = strstr(workptrA,"M  END\n")+strlen("M  END\n")*sizeof(char);
        }

//printf("12");

        memcpy(tmpbuf,workptrA,((int)workptrB-(int)workptrA)*sizeof(char));

//printf("13");

        workptrA=workptrB;

        smiles = ob_molfile_to_canonical_smiles (tmpbuf,0);

        if(smiles == NULL)
        {
            goto smiles_fail2;
        }
        else if (!strlen(smiles))
        {
            free (smiles);
smiles_fail2:
            elog (WARNING, "SMILES generation failed! Trying fallback...");
        }

        if(ob_is_nostruct(tmpbuf) != 0)
        {

            for(i=0; i<mcount; i++)
            {
                if(mols[i] != NULL) pfree(mols[i]);
            }

            free(smiles);
            pfree(mols);
            pfree(rxnfile);
            pfree(tmpbuf);

            elog (ERROR, "NoStructures are not allowed in reactions! Offender was :\n %s",raw_input);
        }

//efa_array = ob_efa_array(smiles);

        mols[i]=new_molecule(smiles,tmpbuf,NULL);

//free(efa_array);

        if (smiles != NULL) free(smiles);
    }

//printf("15");

    pfree(tmpbuf);

//printf("16");

    pfree(rxnfile);

//printf("17");

    /*for(i=0;i<mcount;i++) {
        ////printf("%d\n",mols[i]->len);
        ////printf("%s\n",SMIPTR(mols[i]));
        ////printf("%s\n",MFPTR(mols[i]));
        rxnsize+=mols[i]->len;
    }*/

//printf("%d\n",mcount);

    result = new_reaction(mols,r,p,mode);

    for(i=0; i<mcount; i++)
    {
        if(mols[i] != NULL) pfree(mols[i]);
    }

    pfree(mols);

    return result;
}

/*
* Convert a molecule in text form (V2000, SMILES, InChI) into a molecule
*/
PG_FUNCTION_INFO_V1 (reaction_in);

Datum
reaction_in (PG_FUNCTION_ARGS)
{
    char *input = PG_GETARG_CSTRING (0);
    int size = strlen(input);

    //printf("reaction_in %d\n",size);

    PG_RETURN_REACTION_P (make_reaction(input,size));
}

PG_FUNCTION_INFO_V1 (reaction_in_text);

Datum
reaction_in_text (PG_FUNCTION_ARGS)
{
    char *input = VARDATA(PG_GETARG_TEXT_P (0));
    int size = VARSIZE(PG_GETARG_TEXT_P (0))-VARHDRSZ;

    PG_RETURN_REACTION_P (make_reaction(input,size));
}

PG_FUNCTION_INFO_V1 (reaction_in_varchar);

Datum
reaction_in_varchar (PG_FUNCTION_ARGS)
{
    VarChar *x = PG_GETARG_VARCHAR_P (0);
    char *input = VARDATA(x);
    int size = VARSIZE(x)-VARHDRSZ;

    PG_RETURN_REACTION_P (make_reaction(input,size));
}

PG_FUNCTION_INFO_V1 (reaction_in_bytea);

Datum
reaction_in_bytea (PG_FUNCTION_ARGS)
{
    char *input = VARDATA(PG_GETARG_BYTEA_P (0));
    int size = VARSIZE(PG_GETARG_BYTEA_P (0))-VARHDRSZ;

    PG_RETURN_REACTION_P (make_reaction(input,size));
}

/*
* Output a molecule in cstring form as molfile
*/
PG_FUNCTION_INFO_V1 (reaction_out);

Datum
reaction_out (PG_FUNCTION_ARGS)
{
    REACTION *reaction = PG_GETARG_REACTION_P (0);
    MOLECULE *tmpMol;
    char* offset = MOLARRAYPTR(reaction);
    char *result;
    char *molfile;
    int i, slack;
    char nowstr[13];
    time_t nowbin;
    const struct tm *nowstruct;

    if (time(&nowbin) == (time_t) - 1) elog(WARNING,"Could not get time of day from time()");

    nowstruct = localtime(&nowbin);

    if (strftime(nowstr, 13, "%m%d%Y%H%M", nowstruct) == (size_t) 0) elog(WARNING,"Could not get string from strftime()");

    slack = (37 + (reaction->num_products+reaction->num_reactants)*8)*sizeof(char);

    result = (char *) palloc0 (slack*reaction->datasize*sizeof(char));

    //memset(result,0x0,slack*reaction->datasize*sizeof(char));

    /*  for(i=0;i<reaction->num_reactants;i++) {
      tmpMol = (MOLECULE*) offset;
      if(strstr(SMIPTR(tmpMol),"\r\n") != NULL) {
      strncat(result,SMIPTR(tmpMol),tmpMol->sizesmi-3);
    } else if(strstr(SMIPTR(tmpMol),"\n") != NULL) {
      strncat(result,SMIPTR(tmpMol),tmpMol->sizesmi-2);
    }
      if(i<reaction->num_reactants-1) strncat(result,".",sizeof(char));
      offset+=tmpMol->len*sizeof(char);
    }

    strncat(result,">>",2*sizeof(char));

    for(i=0;i<(reaction->num_products);i++) {
      tmpMol = (MOLECULE*) offset;
      if(strstr(SMIPTR(tmpMol),"\r\n") != NULL) {
      strncat(result,SMIPTR(tmpMol),tmpMol->sizesmi-3);
    } else if(strstr(SMIPTR(tmpMol),"\n") != NULL) {
      strncat(result,SMIPTR(tmpMol),tmpMol->sizesmi-2);
    }
      if(i<reaction->num_products-1) strncat(result,".",sizeof(char));
      offset+=tmpMol->len*sizeof(char);
    }  */

    if(reaction->mode == 1)
        sprintf(result,"$RXN\r\n\r\n      pgchem   %s\r\n\r\n%3d%3d\r\n",nowstr,reaction->num_reactants, reaction->num_products);
    else
        sprintf(result,"$RXN\n\n      pgchem   %s\n\n%3d%3d\n",nowstr,reaction->num_reactants, reaction->num_products);

    for(i=0; i<reaction->num_reactants+reaction->num_products; i++)
    {
        tmpMol = (MOLECULE*) offset;

        if (reaction->mode == 1)
            strncat(result,"$MOL\r\n",strlen("$MOL\r\n")*sizeof(char));
        else
            strncat(result,"$MOL\n",strlen("$MOL\n")*sizeof(char));

        //printf("%s\n",SMIPTR(tmpMol));

        molfile=ob_smiles_to_V2000(SMIPTR(tmpMol));

        strncat(result,molfile,strlen(molfile)+1);

        free(molfile);

        offset+=VARSIZE(tmpMol)*sizeof(char);
    }

    //memset(result,0x0,molecule->sizemf);

    ////printf(molecule->data+molecule->sizesmi);

    //strncpy (result, MFPTR(molecule), molecule->sizemf);

    PG_RETURN_CSTRING (result);
}

/*****************************************************************************
 * Binary Input/Output functions
 *
 * These are optional.
 *****************************************************************************/

PG_FUNCTION_INFO_V1 (reaction_recv);

Datum
reaction_recv (PG_FUNCTION_ARGS)
{
    StringInfo buf = (StringInfo) PG_GETARG_POINTER (0);
    int len = buf->len;
    const char *str = pq_getmsgbytes (buf, len);
    REACTION *result = (REACTION *) palloc0 (len);

    //SET_VARSIZE (result,(buf->len + VARHDRSZ));

    //memset(result,0x0,len);

    memcpy (result, str, len);

    PG_RETURN_POINTER (result);
}

PG_FUNCTION_INFO_V1 (reaction_send);

Datum
reaction_send (PG_FUNCTION_ARGS)
{
    REACTION *reaction = PG_GETARG_REACTION_P (0);

    /*StringInfoData buf;

    pq_begintypsend(&buf);

    pq_sendbytes(&buf,(const char*) molecule,sizeof(molecule));

    PG_RETURN_BYTEA_P(pq_endtypsend(&buf));*/

    PG_RETURN_BYTEA_P(reaction);
}

PG_FUNCTION_INFO_V1 (pgchem_reaction_mol_strip_rxninfo);

Datum
pgchem_reaction_mol_strip_rxninfo (PG_FUNCTION_ARGS)
{
    char *molfile = NULL;
    char *smiles = NULL;
    //unsigned int *efa_array = NULL;
    MOLECULE *retval;
    MOLECULE *arg_molecule = PG_GETARG_MOLECULE_P (0);

    smiles = SMIPTR(arg_molecule);

    molfile =
        ob_smiles_to_V2000 (smiles);

    //efa_array = ob_efa_array(smiles);

    retval = new_molecule (smiles, molfile,NULL);

    //free(efa_array);

    free (molfile);

    PG_RETURN_MOLECULE_P (retval);
}

