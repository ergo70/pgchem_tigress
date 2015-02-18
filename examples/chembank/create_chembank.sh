#!/bin/sh
SUPERUSER=postgres
PGPASSWORD=postgres
DATABASE=mol
export PGPASSWORD
psql -U ${SUPERUSER}  ${DATABASE} -f user.sql.chembank
psql -U ${SUPERUSER}  ${DATABASE} -f schema.sql.chembank
psql -U ${SUPERUSER}  ${DATABASE} -f chembank_exampletable.sql
psql -U ${SUPERUSER}  ${DATABASE} -f tables.sql.chembank
psql -U ${SUPERUSER}  ${DATABASE} -f triggers.sql.chembank
