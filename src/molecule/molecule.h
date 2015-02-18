#include "postgres.h"
#include "molecule_limits.h"

/*#ifndef SET_VARSIZE // < 8.3 compatibility
#define SET_VARSIZE(v,l) (VARATT_SIZEP(v) = (l))
#endif*/

#ifdef BUILD_WITH_INDIGO
typedef int vFPsi __attribute__ ((vector_size (512)));
#else
typedef int vFPsi __attribute__ ((vector_size (256)));
#endif

typedef union
{
    vFPsi v;
#ifdef BUILD_WITH_INDIGO
    uint32 dwords[IN_FPSIZE];
    uint64 qwords[IN_FPSIZE/2];
    char bytes[IN_FPSIZE*sizeof(uint32)];
#else
    uint32 dwords[OB_FPSIZE];
    uint64 qwords[OB_FPSIZE/2];
    char bytes[OB_FPSIZE*sizeof(uint32)];
#endif
} MOLFP;

typedef struct
{
    char vl_len[4];
    uint32 sizeo;
    uint32 compressed_sizeo;
    uint32 sizesmi;
    uint32 sizeanc;
    uint32 disconnected;
    uint32 original_format;
    MOLFP fp;
    char inchikey[INCHIKEYSZ];
    char data[1];
} MOLECULE;

#define SMIPTR(x)		 	         ((char*)(x->data))
#define CIPTR(x)		 		     ((char*)(x->data+x->sizesmi))
#define ANCPTR(x)		 		     ((char*)(x->data+x->sizesmi+x->compressed_sizeo))
#define CALCDATASZ(compressed_sizeo, sizesmi, sizeanc)  (sizeof(MOLECULE) + compressed_sizeo + sizesmi + sizeanc)

#define PG_GETARG_MOLFP_P(n)         (MOLFP *) DatumGetPointer(PG_GETARG_DATUM(n))
#define PG_RETURN_MOLFP_P(n)         PG_RETURN_POINTER(n)

#define DatumGetMoleculeP(n)         ((MOLECULE *) PG_DETOAST_DATUM(n))
#define PG_GETARG_MOLECULE_P(n)      DatumGetMoleculeP(PG_GETARG_DATUM(n))
#define PG_RETURN_MOLECULE_P(n)      PG_RETURN_POINTER(n)

//#define FORMAT_OTHER 666

#ifdef BUILD_WITH_INDIGO
MOLECULE *new_molecule (char *smiles, char *original_data);
#else
MOLECULE *new_molecule (char *smiles, char *original_data, char *ancillarydata);
#endif
