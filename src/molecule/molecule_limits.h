#define INCHIKEYSZ 27 // bytes

#ifdef BUILD_WITH_INDIGO
#define IN_SUBFPSIZE 54 //dwords
#define IN_SIMOFFSET IN_SUBFPSIZE //dwords
#define IN_SIMFPSIZE 16 //dwords
#define IN_FPSIZE 120 //dwords
#define IN_FPSTART 13 //bytes
#endif

#define OB_FPSIZE2 32  //dwords
#define OB_FPSIZE3 16   //dwords
#define OB_FPSIZE_MACCS 8 //dwords
//#define OB_FPSIZE (OB_FPSIZE2+OB_FPSIZE3) //dwords
#define OB_FPSIZE 64 //dwords
#define OB_OFFSET OB_FPSIZE2 //dwords

#define FORMAT_V2000 1
#define FORMAT_V3000 2
#define FORMAT_INCHI 3
#define FORMAT_SMILES 4
