@echo off
set SUPERUSER=postgres
set PGPASSWORD=postgres
set DATABASE=mol

psql -U %SUPERUSER%  %DATABASE% -f user.sql.chembank
psql -U %SUPERUSER%  %DATABASE% -f schema.sql.chembank
psql -U %SUPERUSER%  %DATABASE% -f chembank_exampletable.sql
psql -U %SUPERUSER%  %DATABASE% -f tables.sql.chembank
psql -U %SUPERUSER%  %DATABASE% -f triggers.sql.chembank
