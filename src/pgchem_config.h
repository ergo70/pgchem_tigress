/************************************************************************
 * pgchem_config.h compile-time configuration settings
 *
 * Copyright (c) 2016 by Ernst-G. Schmid
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

#ifdef WIN32
#define PGCHEM_VERSION "pgchem_tigress_3.3 GiST WIN32"
#else
#define PGCHEM_VERSION "pgchem_tigress_3.3 GiST POSIX"
#endif
