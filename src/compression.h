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

inline COMPRESSED_DATA *compress_data (const char *data, uint32 decompressed_size);
inline DECOMPRESSED_DATA *decompress_data (const char *data, uint32 compressed_size, uint32 decompressed_size);
