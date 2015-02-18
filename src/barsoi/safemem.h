/************************************************************************
 * safemem.h checking memory allocation functions
 *  
 * Copyright (c) 2005 by Ernst-G. Schmid
 *   
 * This file is part of the xchem::tigress project.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the lesser GNU General Public License as published by
 * the Free Software Foundation version 2.1 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * lesser GNU General Public License for more details.
 ************************************************************************/

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>

void *safe_calloc (size_t nmemb, size_t size);
void *safe_malloc (size_t size);
