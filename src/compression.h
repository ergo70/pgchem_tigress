#include "postgres.h"
#include <zlib.h>

typedef struct
{
  uint32 compressed_size;
  char *compressed_data;
} COMPRESSED_DATA;

typedef struct
{
  uint32 decompressed_size;
  char *decompressed_data;
} DECOMPRESSED_DATA;

COMPRESSED_DATA *compress_data (const char *data, uint32 decompressed_size);
DECOMPRESSED_DATA *decompress_data (const char *data, uint32 compressed_size, uint32 decompressed_size);
