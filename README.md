## Build notice

If you build against PostgreSQL 10.x, make sure to adjust the includes for md5.h in molecule_io.c, molecule_gist.c, reaction_io.c, and reaction_gist.c accordingly. md5.h has moved from libpg/md5.h to common/md5.h!
