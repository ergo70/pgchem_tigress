#include "postgres.h"
#include "molecule/molecule.h"
#include "molecule/molecule_limits.h"

/*#ifndef SET_VARSIZE // < 8.3 compatibility
#define SET_VARSIZE(v,l) (VARATT_SIZEP(v) = (l))
#endif*/

typedef struct
{
    uint32 fp[2*OB_FPSIZE2];
} RXNFP;

typedef struct
{
    char vl_len[4];
    uint32 datasize;
    uint32 num_reactants;
    uint32 num_products;
    uint32 mode;
    uint32 fp[2*OB_FPSIZE2];
    char data[1];
} REACTION;

//#define RXNHDRSZ		             (4*sizeof(int4))
#define MOLARRAYPTR(x)		 	         ((char*)(x->data))
//#define MFPTR(x)		 		     ((char*)(x->data+x->sizesmi))
//#define RCALCDATASZ(datasize)  (((2*FPSIZE2)*sizeof(uint32)) + RXNHDRSZ + datasize)
#define RCALCDATASZ(datasize)  ((sizeof(REACTION)) + datasize - sizeof(char))

#define PG_GETARG_RXNFP_P(n)         (RXNFP *) DatumGetPointer(PG_GETARG_DATUM(n))
#define PG_RETURN_RXNFP_P(n)         PG_RETURN_POINTER(n)

#define DatumGetReactionP(n)         ((REACTION *) PG_DETOAST_DATUM(n))
#define PG_GETARG_REACTION_P(n)      DatumGetReactionP(PG_GETARG_DATUM(n))
#define PG_RETURN_REACTION_P(n)      PG_RETURN_POINTER(n)

