/************************************************************************
 * compression.c compression support functions
 *
 * Copyright (c) 2007,2012 by Ernst-G. Schmid
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
#include "compression.h"

inline COMPRESSED_DATA *compress_data (const char *data, uint32 decompressed_size)
{
    Byte *compressed_data=NULL;
    uLong compressed_size=0UL;
    int err=0;
    COMPRESSED_DATA *retval=NULL;

    //elog (INFO, "Size :\n %u",decompressed_size);

    compressed_size = compressBound(decompressed_size);

    //elog (INFO, "Size :\n %lu",decompressed_size);

    compressed_data = (Byte*) palloc0(compressed_size*sizeof(Byte));

    err = compress2(compressed_data, &compressed_size, (const Bytef*) data, (uLong) decompressed_size, Z_BEST_COMPRESSION);

    if(err!=Z_OK)
    {
        elog (ERROR, "Compression failed: %i", err);
    }

    retval = (COMPRESSED_DATA*) palloc0(sizeof(COMPRESSED_DATA));

    //elog (INFO, "Size :\n %lu",compressed_size);
    retval->compressed_size=(uint32) compressed_size;
    retval->compressed_data=(char*) compressed_data;

    //elog (INFO, "Size :\n %u",retval->compressed_size);

    return retval;
}


inline DECOMPRESSED_DATA *decompress_data (const char *data, uint32 compressed_size, uint32 decompressed_size)
{
    Byte *decompressed_data=NULL;
    uLong size_after_decompression=decompressed_size;
    int err=0;
    DECOMPRESSED_DATA *retval=NULL;

    /*data_to_decompress = malloc(compressed_size*sizeof(Bytef));*/
    decompressed_data = palloc0(decompressed_size*sizeof(Byte));

    /*memset(data_to_decompress,0x0, compressed_size);
    memcpy(data_to_decompress,data,compressed_size);*/

    err = uncompress(decompressed_data, &size_after_decompression, (const Bytef*) data, compressed_size);

    //free(data_to_decompress);

    if(err!=Z_OK)
    {
        elog (ERROR, "Decompression failed: %i", err);
    }

    if(size_after_decompression != decompressed_size)
    {
        elog (WARNING, "Decompressed size != expected size: %u != %u",(uint32) size_after_decompression, decompressed_size);
    }

    retval = (DECOMPRESSED_DATA*) palloc(sizeof(DECOMPRESSED_DATA));

    retval->decompressed_size=(uint32) size_after_decompression;
    retval->decompressed_data=(char*) decompressed_data;

    return retval;
}
