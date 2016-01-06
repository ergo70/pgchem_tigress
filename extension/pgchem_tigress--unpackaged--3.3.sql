/**********************************************************************
 * pgchem_setup.sql unified type, operator and function setup file for pgchem GiST
 *
 * Copyright (c) 2004,2016 by Ernst-G. Schmid
 *
 * This file is part of the pgchem::tigress project.
 * For more information, see
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 2 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 ************************************************************************/


DROP TYPE IF EXISTS molecule CASCADE;
DROP TYPE IF EXISTS molfp CASCADE;

CREATE TYPE molecule;
CREATE TYPE molfp;

CREATE OR REPLACE FUNCTION molfp_in(cstring)
    RETURNS molfp
    AS '$libdir/libpgchem'
    LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION molfp_out(molfp)
    RETURNS cstring
   AS '$libdir/libpgchem'
    LANGUAGE C IMMUTABLE STRICT;

CREATE TYPE molfp (
   input = molfp_in,
   output = molfp_out,
   alignment = int4,
   /* OB: internallength = 256 */
   /* IN: internallength = 512 */
   internallength = 256,
   storage = PLAIN
   );

CREATE OR REPLACE FUNCTION molecule_in(cstring)
    RETURNS molecule
    AS '$libdir/libpgchem'
    LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION molecule_out(molecule)
    RETURNS cstring
   AS '$libdir/libpgchem'
    LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION molecule_recv(internal)
   RETURNS molecule
  AS '$libdir/libpgchem'
   LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION molecule_send(molecule)
   RETURNS bytea
  AS '$libdir/libpgchem'
   LANGUAGE C IMMUTABLE STRICT;

CREATE TYPE molecule (
   input = molecule_in,
   output = molecule_out,
   alignment = int4,
   --receive = molecule_recv,
   --send = molecule_send,
   internallength = VARIABLE,
   storage = MAIN
   );

CREATE OR REPLACE FUNCTION molfp_compress(internal)
RETURNS internal
AS '$libdir/libpgchem'
LANGUAGE C;

CREATE OR REPLACE FUNCTION molfp_decompress(internal)
RETURNS internal
AS '$libdir/libpgchem'
LANGUAGE C;

CREATE OR REPLACE FUNCTION molfp_penalty(internal,internal,internal)
RETURNS internal
AS '$libdir/libpgchem'
LANGUAGE C WITH (isstrict);

CREATE OR REPLACE FUNCTION molfp_picksplit(internal,internal)
RETURNS internal
AS '$libdir/libpgchem'
LANGUAGE C WITH (isstrict);

CREATE OR REPLACE FUNCTION molfp_union(internal, internal)
RETURNS internal
AS '$libdir/libpgchem'
LANGUAGE C;

CREATE OR REPLACE FUNCTION molfp_same(internal, internal, internal)
RETURNS internal
AS '$libdir/libpgchem'
LANGUAGE C;

CREATE OR REPLACE FUNCTION molfp_consistent(internal,internal,int4)
RETURNS bool
AS '$libdir/libpgchem'
LANGUAGE C;

CREATE OR REPLACE FUNCTION molecule_contained_in(molecule,molecule)
RETURNS bool
AS '$libdir/libpgchem'
LANGUAGE C with (isstrict);

CREATE OR REPLACE FUNCTION molecule_contains(molecule,molecule)
RETURNS bool
AS '$libdir/libpgchem'
LANGUAGE C with (isstrict);

CREATE OR REPLACE FUNCTION molecule_maybe_contained_in(molecule,molecule)
RETURNS bool
AS '$libdir/libpgchem'
LANGUAGE C with (isstrict);

/*CREATE OR REPLACE FUNCTION molecule_maybe_contains(molecule,molecule)
RETURNS bool
AS '$libdir/libpgchem'
LANGUAGE C with (isstrict);*/

CREATE OR REPLACE FUNCTION molecule_equals(molecule,molecule)
RETURNS bool
AS '$libdir/libpgchem'
LANGUAGE C with (isstrict);

CREATE OR REPLACE FUNCTION molecule_similarity(molecule,molecule)
RETURNS double precision
AS '$libdir/libpgchem'
LANGUAGE C with (isstrict);

CREATE OPERATOR < (
		 LEFTARG = molecule,
		 RIGHTARG = molecule,
		 PROCEDURE = molecule_maybe_contained_in,
		 COMMUTATOR = <,
		 RESTRICT = contsel,
		 JOIN = contjoinsel
);

/*CREATE OPERATOR > (
		 LEFTARG = molecule,
		 RIGHTARG = molecule,
		 PROCEDURE = molecule_maybe_contains,
		 COMMUTATOR = <,
		 RESTRICT = contsel,
		 JOIN = contjoinsel
);*/

CREATE OPERATOR <= (
		 LEFTARG = molecule,
		 RIGHTARG = molecule,
		 PROCEDURE = molecule_contained_in,
                 COMMUTATOR = >=,
		 RESTRICT = contsel,
		 JOIN = contjoinsel
);

CREATE OPERATOR >= (
		 LEFTARG = molecule,
		 RIGHTARG = molecule,
		 PROCEDURE = molecule_contains,
                 COMMUTATOR = <=,
		 RESTRICT = contsel,
		 JOIN = contjoinsel
);

CREATE OPERATOR = (
		 LEFTARG = molecule,
		 RIGHTARG = molecule,
		 PROCEDURE = molecule_equals,
		 COMMUTATOR = =,
		 RESTRICT = eqsel,
		 JOIN = eqjoinsel
);

CREATE OPERATOR @ (
    leftarg = molecule,
    rightarg = molecule,
    procedure = molecule_similarity
);

CREATE OPERATOR CLASS gist_molecule_ops
DEFAULT FOR TYPE molecule USING gist
AS
	/*OPERATOR	3	@	 ,*/
        OPERATOR        6       =        ,
        OPERATOR        7       >=       ,
        OPERATOR        8       <=       ,
        OPERATOR	9	<	,
        /*OPERATOR	10	>	,*/
        FUNCTION        1       molfp_consistent (internal, internal, int4),
        FUNCTION        2       molfp_union (internal, internal),
        FUNCTION        3       molfp_compress (internal),
        FUNCTION        4       molfp_decompress (internal),
        FUNCTION        5       molfp_penalty (internal, internal, internal),
        FUNCTION        6       molfp_picksplit (internal, internal),
        FUNCTION        7       molfp_same (internal, internal, internal),
        STORAGE		molfp;

CREATE OR REPLACE FUNCTION molecule_in(text)
  RETURNS molecule AS
'$libdir/libpgchem', 'molecule_in_text'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION molecule_in(character varying)
  RETURNS molecule AS
'$libdir/libpgchem', 'molecule_in_varchar'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION molecule_in(bytea)
  RETURNS molecule AS
'$libdir/libpgchem', 'molecule_in_bytea'
  LANGUAGE C IMMUTABLE STRICT;

CREATE CAST (text AS molecule)
  WITH FUNCTION molecule_in(text);

CREATE CAST (character varying AS molecule)
  WITH FUNCTION molecule_in(character varying);

CREATE CAST (bytea AS molecule)
  WITH FUNCTION molecule_in(bytea);

CREATE OR REPLACE FUNCTION add_hydrogens(molecule, boolean, boolean, double precision default 7.4)
  RETURNS molecule AS
'$libdir/libpgchem', 'pgchem_add_hydrogens'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION exactmass(molecule)
  RETURNS double precision AS
'$libdir/libpgchem', 'pgchem_exactmass'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION fgroup_codes(molecule)
  RETURNS text AS
'$libdir/libpgchem', 'pgchem_fgroup_codes_a'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION is_2d(molecule)
  RETURNS boolean AS
'$libdir/libpgchem', 'pgchem_2D'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION is_3d(molecule)
  RETURNS boolean AS
'$libdir/libpgchem', 'pgchem_3D'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION is_chiral(molecule)
  RETURNS boolean AS
'$libdir/libpgchem', 'pgchem_is_chiral'
  LANGUAGE C IMMUTABLE STRICT;

/*CREATE OR REPLACE FUNCTION is_nostruct(molecule)
  RETURNS boolean AS
'$libdir/libpgchem', 'pgchem_is_nostruct'
  LANGUAGE C IMMUTABLE STRICT;*/

CREATE OR REPLACE FUNCTION inchi(molecule)
  RETURNS text AS
'$libdir/libpgchem', 'pgchem_molecule_to_inchi'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION canonical_smiles(molecule,boolean)
  RETURNS text AS
'$libdir/libpgchem', 'pgchem_molecule_to_canonical_smiles'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION smiles(molecule, boolean)
  RETURNS text AS
'$libdir/libpgchem', 'pgchem_molecule_to_smiles'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION v3000(molecule)
  RETURNS text AS
'$libdir/libpgchem', 'pgchem_molecule_to_V3000'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION v2000(molecule)
  RETURNS text AS
'$libdir/libpgchem', 'pgchem_molecule_to_V2000'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION molformula(molecule)
  RETURNS text AS
'$libdir/libpgchem', 'pgchem_hillformula'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION molweight(molecule)
  RETURNS double precision AS
'$libdir/libpgchem', 'pgchem_molweight'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION molkeys_long(molecule, boolean, boolean, boolean)
  RETURNS text AS
'$libdir/libpgchem', 'pgchem_ms_fingerprint_long_a'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION molkeys_long(molecule)
  RETURNS text AS
$BODY$
BEGIN
RETURN molkeys_long($1,false,false,false);
END;
$BODY$
  LANGUAGE 'plpgsql' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION molkeys_short(molecule)
  RETURNS text AS
'$libdir/libpgchem', 'pgchem_ms_fingerprint_short_a'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION number_of_atoms(molecule)
  RETURNS integer AS
'$libdir/libpgchem', 'pgchem_num_atoms'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION number_of_bonds(molecule)
  RETURNS integer AS
'$libdir/libpgchem', 'pgchem_num_bonds'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION number_of_heavy_atoms(molecule)
  RETURNS integer AS
'$libdir/libpgchem', 'pgchem_num_heavy_atoms'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION number_of_rotatable_bonds(molecule)
  RETURNS integer AS
'$libdir/libpgchem', 'pgchem_num_rotatable_bonds'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION pgchem_barsoi_version()
  RETURNS cstring AS
'$libdir/libpgchem', 'pgchem_barsoi_version'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION pgchem_version()
  RETURNS cstring AS
'$libdir/libpgchem', 'pgchem_version'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION remove_hydrogens(molecule, boolean)
  RETURNS molecule AS
'$libdir/libpgchem', 'pgchem_remove_hydrogens'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION strip_salts(molecule, boolean)
  RETURNS molecule AS
'$libdir/libpgchem', 'pgchem_strip_salts'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION total_charge(molecule)
  RETURNS integer AS
'$libdir/libpgchem', 'pgchem_total_charge'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION original_format(molecule)

RETURNS cstring AS

'$libdir/libpgchem', 'pgchem_original_format'

 LANGUAGE C IMMUTABLE STRICT
;

CREATE OR REPLACE FUNCTION validate_cas_no(character varying)
  RETURNS boolean AS
$BODY$
DECLARE checksum_from_cas_no varchar;
DECLARE cas_no_left varchar;
DECLARE cas_no_right varchar;
DECLARE cas_no_full varchar;
DECLARE tmpsum int;
DECLARE position_multiplier int;
DECLARE caslen int;
BEGIN
caslen:=length($1);

IF caslen<5 OR caslen>12 THEN RETURN FALSE;
END IF;

checksum_from_cas_no:=split_part($1,'-',3)::int;
cas_no_left:=split_part($1,'-',1);
cas_no_right:=split_part($1,'-',2);
cas_no_full:=cas_no_left || cas_no_right;

if(length(cas_no_left)>7 OR length(cas_no_right)>2 OR length(checksum_from_cas_no)!=1) THEN
return false;
END IF;

caslen:=length(cas_no_full);
tmpsum:=0;
position_multiplier:=1;

 FOR i IN REVERSE caslen..1 LOOP
  tmpsum:=tmpsum+substr(cas_no_full,i,1)::int*position_multiplier;
  position_multiplier:=position_multiplier+1;
 END LOOP;
 RETURN tmpsum % 10 = checksum_from_cas_no::int;
END;
$BODY$
  LANGUAGE 'plpgsql' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION logP(molecule)
  RETURNS double precision AS
'$libdir/libpgchem', 'pgchem_logP'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION MR(molecule)
  RETURNS double precision AS
'$libdir/libpgchem', 'pgchem_MR'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION TPSA(molecule)
  RETURNS double precision AS
'$libdir/libpgchem', 'pgchem_TPSA'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION SMARTSmatch(text, molecule)
  RETURNS boolean AS
'$libdir/libpgchem', 'pgchem_smartsfilter'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION SMARTSmatch_count(text, molecule)
  RETURNS integer AS
'$libdir/libpgchem', 'pgchem_smartsfilter_count'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION inchikey(molecule)
  RETURNS text AS
'$libdir/libpgchem', 'pgchem_molecule_to_inchikey'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION fpstring(molecule)
  RETURNS bit varying AS
'$libdir/libpgchem', 'pgchem_fp_out'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION number_of_h_acceptors(molecule)
  RETURNS integer AS
'$libdir/libpgchem', 'pgchem_num_H_acceptors'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION number_of_h_donors(molecule)
  RETURNS integer AS
'$libdir/libpgchem', 'pgchem_num_H_donors'
  LANGUAGE C IMMUTABLE STRICT;

-- Function: isotope_pattern(molecule, integer, double precision)

-- DROP FUNCTION isotope_pattern(molecule, integer, double precision);

CREATE OR REPLACE FUNCTION isotope_pattern(
    IN molecule,
    IN integer DEFAULT 0,
    IN double precision DEFAULT 100.0)
  RETURNS TABLE("m/z" double precision, intensity double precision, intensity_normalized double precision, mdtbp integer) AS
'$libdir/libpgchem', 'pgchem_isotope_pattern'
  LANGUAGE c IMMUTABLE STRICT
  COST 1
  ROWS 1000;

CREATE OR REPLACE FUNCTION lipinsky(mol molecule)
  RETURNS text AS
$BODY$
DECLARE parameters text;
BEGIN
parameters := '';

IF number_of_H_donors(mol) > 5 THEN
parameters := parameters || 'A';
END IF;
IF molweight(mol) > 500 THEN
parameters := parameters || 'B';
END IF;
IF logP(mol) > 5.0 THEN
parameters := parameters || 'C';
END IF;
IF number_of_H_acceptors(mol) > 10 THEN
parameters := parameters || 'D';
END IF;

RETURN parameters;
END;
$BODY$
  LANGUAGE 'plpgsql' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION mutabor(molecule)
  RETURNS molecule AS
'$libdir/libpgchem', 'pgchem_mutate_fp'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION delebor(molecule)
  RETURNS molecule AS
'$libdir/libpgchem', 'pgchem_blank_fp'
  LANGUAGE C IMMUTABLE STRICT;

-- Function: fp2string(molecule)

-- DROP FUNCTION fp2string(molecule);

CREATE OR REPLACE FUNCTION fp2string(struct molecule)
  RETURNS bit varying AS
$BODY$
DECLARE fp bit varying;
BEGIN
fp := fpstring(struct);

IF (length(fp) != 1536) THEN
RAISE EXCEPTION 'OPERATION NOT SUPPORTED';
END IF;

fp := substring(fp for 1024);
RETURN fp;
END;
$BODY$
  LANGUAGE plpgsql IMMUTABLE STRICT
  COST 100;

CREATE OR REPLACE FUNCTION fp3string(struct molecule)
  RETURNS bit varying AS
$BODY$
DECLARE fp bit varying;
BEGIN
fp := fpstring(struct);

IF (length(fp) != 1536) THEN
RAISE EXCEPTION 'OPERATION NOT SUPPORTED';
END IF;

fp := substring(fp from 1025);
RETURN fp;
END;
$BODY$
  LANGUAGE plpgsql IMMUTABLE STRICT
  COST 100;

/*CREATE OR REPLACE FUNCTION nse(tablename text, columnname text)
  RETURNS SETOF int AS
$BODY$
DECLARE fpAllOn bit varying;
DECLARE fpAllOff bit varying;
DECLARE row record;
DECLARE i int;
BEGIN
fpAllOn := B'11111111111111111111111111111111';
fpAllOff := B'00000000000000000000000000000000';

FOR row IN EXECUTE 'SELECT fp3string(' || columnname || ') as fp from ' || tablename LOOP
fpAllOn = fpAllOn & row.fp;
fpAllOff = fpAllOff | row.fp;
END LOOP;

FOR i IN 1..32 LOOP

IF (substring(fpAllOn from i for 1)='1') THEN RETURN NEXT i; END IF;
IF (substring(fpAllOff from i for 1)='0') THEN RETURN NEXT i; END IF;

END LOOP;

RETURN;
END;
$BODY$
  LANGUAGE 'plpgsql' IMMUTABLE STRICT;*/

CREATE OR REPLACE FUNCTION nbits_set(bit varying)
  RETURNS integer AS
'$libdir/libpgchem', 'pgchem_nbits_set'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION disconnected(molecule)
  RETURNS boolean AS
'$libdir/libpgchem', 'pgchem_disconnected'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION check_fingerprint_optimization(t_schema text, t_name text, s_column text)
  RETURNS double precision AS
$BODY$
DECLARE i float8;
DECLARE tablesize integer;
DECLARE t_q_name text;
DECLARE t_b_name text;
BEGIN

t_q_name := t_schema || '.q_c_' || t_name;
t_b_name := t_schema || '.' || t_name;

EXECUTE 'DROP TABLE IF EXISTS '||t_q_name;

EXECUTE 'CREATE TABLE ' || t_q_name || ' AS SELECT ' || s_column || ' FROM ' || t_b_name || ' WHERE fp2string('||s_column||') IN (SELECT fp2string('||s_column||') FROM ' || t_b_name || ' GROUP BY fp2string('||s_column||') HAVING (COUNT(fp2string('||s_column||'))>1))';

EXECUTE 'UPDATE '||t_q_name||' SET structure=mutabor('||s_column||')';

EXECUTE 'SELECT COUNT(1) FROM (SELECT count(fpstring('||s_column||')) as count FROM '||t_q_name||' GROUP BY fpstring('||s_column||')) as t WHERE t.count = 1' into i;

EXECUTE 'SELECT COUNT(1) FROM '||t_q_name into tablesize;

EXECUTE 'DROP TABLE IF EXISTS '||t_q_name;

RETURN i/tablesize;
END;
$BODY$
  LANGUAGE 'plpgsql' VOLATILE;

CREATE OR REPLACE FUNCTION create_random_sample_q_table(t_schema text, t_name text, s_column text)
  RETURNS void AS
$BODY$
DECLARE n numeric;
DECLARE i integer;
DECLARE largeN integer;
DECLARE ts_size integer;
DECLARE t numeric;
DECLARE p numeric;
DECLARE q numeric;
DECLARE d numeric;
DECLARE numerator numeric;
DECLARE denominator numeric;
DECLARE t_q_name text;
DECLARE t_s_q_name text;
BEGIN

t:=1.96;
p:=0.5;
q:=0.5;
d:=0.05;

t_q_name := t_schema || '.q_' || t_name;
t_s_q_name := t_schema || '.s_q_' || t_name;

EXECUTE 'DROP TABLE IF EXISTS '||t_s_q_name;

EXECUTE 'CREATE TABLE '||t_s_q_name||' ('||s_column||' molecule) WITH (OIDS=FALSE)';

EXECUTE 'SELECT count(1) from '||t_q_name into largeN;

numerator:=t^2*(p*q);

denominator:=d^2;

n:=round((numerator/denominator),0);

--raise info 'n=%',n;

IF (n/largeN >= 0.05) THEN
n := round(((numerator/denominator)/((((numerator/denominator)-1.0)/largeN)+1.0))::numeric,0);
END IF;

--raise info 'n=%',n;

FOR i IN 1..n LOOP

EXECUTE 'INSERT INTO '||t_s_q_name||' (SELECT structure from '||t_q_name||' WHERE fp2string('||s_column||')=(SELECT fp2string('||s_column||') FROM '||t_q_name||' ORDER BY RANDOM() LIMIT 1))';

END LOOP;

EXECUTE 'SELECT count(1) from '||t_s_q_name into ts_size;

--raise info 'ts_size=%',ts_size;

IF (ts_size>=largeN) THEN
BEGIN
RAISE INFO 'Sampling has no effect. Retaining original table.';
EXECUTE 'DROP TABLE IF EXISTS '||t_s_q_name;
END;
ELSE
BEGIN
EXECUTE 'DROP TABLE IF EXISTS '||t_q_name;
EXECUTE 'ALTER TABLE '||t_s_q_name||' RENAME TO q_' || t_name;
END;
END IF;

RETURN;
END;
$BODY$
  LANGUAGE 'plpgsql' VOLATILE;

CREATE OR REPLACE FUNCTION optimize_fingerprint(t_schema text, t_name text, s_column text, algorithm text, use_sampling boolean, basedict text, p_limit integer)
  RETURNS double precision AS
$BODY$
BEGIN

IF (algorithm='LP') THEN
BEGIN
PERFORM optimize_fingerprint_step_one(t_schema, t_name, s_column, use_sampling , basedict ,p_limit);
RETURN optimize_fingerprint_step_two(t_schema, t_name, s_column);
END;
ELSIF (algorithm='GA') THEN
RAISE EXCEPTION 'Genetic optimization only externally supported';
ELSE RAISE EXCEPTION 'Algorithm % not supported',algorithm;
END IF;

RAISE INFO 'Optimization on % complete. Achieved rate %',t_schema||'.'||t_name, rate;

RETURN 0.0;
END;
$BODY$
  LANGUAGE 'plpgsql' VOLATILE;

CREATE OR REPLACE FUNCTION optimize_fingerprint_step_one(t_schema text, t_name text, s_column text, use_sampling boolean, basedict text, p_limit integer)
  RETURNS void AS
$BODY$
DECLARE expansionarray smallint[9];
DECLARE curr_row record;
DECLARE curr_row_i record;
DECLARE tablesize integer;
DECLARE hitcount integer;
DECLARE matchcount integer;
DECLARE expcount integer;
DECLARE j integer;
DECLARE e float8;
DECLARE x float8;
DECLARE t_dict_tmp_name text;
DECLARE t_dict_name text;
DECLARE t_sel_name text;
DECLARE t_q_name text;
DECLARE t_b_name text;
DECLARE t_h_name text;
DECLARE arraytext text;
BEGIN

t_dict_tmp_name := t_schema || '.' || t_name || '_dictionary_tmp';
t_dict_name := t_schema || '.' || t_name || '_dictionary';
t_sel_name := t_schema || '.' || t_name || '_selectivity';
t_q_name := t_schema || '.q_' || t_name;
t_b_name := t_schema || '.' || t_name;
t_h_name := t_schema || '.' || t_name || '_helper';

EXECUTE 'DROP TABLE IF EXISTS ' || t_dict_tmp_name;

EXECUTE 'DROP TABLE IF EXISTS ' || t_dict_name;

EXECUTE 'DROP TABLE IF EXISTS ' || t_sel_name;

EXECUTE 'CREATE TABLE ' || t_dict_tmp_name || ' (id serial NOT NULL, "SMARTS" text NOT NULL) WITH (OIDS=FALSE)';

EXECUTE 'COPY ' || t_dict_tmp_name || ' ("SMARTS") FROM ' ||  quote_literal(basedict);

EXECUTE 'CREATE TABLE ' || t_dict_name || ' (id serial NOT NULL, "SMARTS" text NOT NULL) WITH (OIDS=FALSE)';

EXECUTE 'INSERT INTO ' || t_dict_name ||  ' ("SMARTS") (SELECT distinct("SMARTS") FROM ' || t_dict_tmp_name || ')';

EXECUTE 'DROP TABLE IF EXISTS ' || t_dict_tmp_name;

EXECUTE 'CREATE TABLE ' || t_sel_name || '(pattern_id integer NOT NULL,expansion smallint NOT NULL DEFAULT 0,coverage double precision NOT NULL DEFAULT 0.0,weight double precision,exparray integer[],_e double precision,_x double precision,CONSTRAINT pattern_id_unique UNIQUE (pattern_id)) WITH (OIDS=FALSE)';

EXECUTE 'DROP TABLE IF EXISTS ' || t_q_name;

EXECUTE 'CREATE TABLE ' || t_q_name || ' AS SELECT ' || s_column || ' FROM ' || t_b_name || ' WHERE fp2string('||s_column||') IN (SELECT fp2string('||s_column||') FROM ' || t_b_name || ' GROUP BY fp2string('||s_column||') HAVING (COUNT(fp2string('||s_column||'))>1))';

IF (use_sampling) THEN
	PERFORM create_random_sample_q_table(t_schema, t_name, s_column);
END IF;

EXECUTE 'DELETE FROM '|| t_q_name ||' WHERE inchikey(' || s_column || ') IN ((SELECT inchikey(' || s_column || ') FROM '|| t_q_name ||' GROUP BY inchikey(' || s_column || ') HAVING (COUNT(inchikey(' || s_column || '))>1)))';

EXECUTE 'SELECT count(1) FROM '|| t_q_name into tablesize;

e := tablesize / 9;

FOR curr_row IN EXECUTE 'SELECT id, "SMARTS" FROM ' || t_dict_name LOOP

x := 0.0;

expansionarray := '{0,0,0,0,0,0,0,0,0}';

hitcount := 0;

expcount := 0;

FOR curr_row_i IN EXECUTE 'SELECT smartsmatch_count('||quote_literal(curr_row."SMARTS")||', ' || s_column || ') as mc FROM '|| t_q_name ||' WHERE smartsmatch_count('||quote_literal(curr_row."SMARTS")||', ' || s_column || ')>0' LOOP

	hitcount:=hitcount+1;

	matchcount := curr_row_i.mc;

	IF (matchcount > 8) THEN matchcount := 8; END IF;

	expansionarray[matchcount+1] := expansionarray[matchcount+1] + 1;

END LOOP;

IF (hitcount > 0) THEN

expansionarray[1] = tablesize - hitcount;

FOR i IN 1..9 LOOP

	IF (expansionarray[i] > 0) THEN
		BEGIN
		expcount := expcount + 1;
		x := x+(power(expansionarray[i] - e,2)/e);
		END;
	END IF;

END LOOP;

arraytext :='{';

FOR i IN 1..9 LOOP
	arraytext := arraytext || expansionarray[i] || ',';
END LOOP;

arraytext := substring(arraytext FROM 1 for length(arraytext)-1) || '}';

EXECUTE 'INSERT INTO '||t_sel_name||' (pattern_id, expansion, coverage,exparray,_e,_x) VALUES ('||curr_row.id||','||expcount||','||hitcount/tablesize::float8||','||quote_literal(arraytext)||','||e||','||x||')';

END IF;

END LOOP;

EXECUTE 'UPDATE '||t_sel_name||' SET weight = _x';

EXECUTE 'CREATE TABLE '||t_h_name||' AS SELECT "SMARTS", expansion from '||t_sel_name||','|| t_dict_name||'  where id=pattern_id and pattern_id IN (SELECT pattern_id FROM '||t_sel_name||' ORDER BY weight ASC LIMIT '||p_limit||')' ;

EXECUTE 'COPY (SELECT ''#Comments after SMARTS'' UNION SELECT * from (SELECT "SMARTS" FROM '||t_h_name||') t) TO '||quote_literal('c:/tigress/obdata/dictionary.txt');

EXECUTE 'DROP TABLE IF EXISTS '||t_h_name;

EXECUTE 'DROP TABLE IF EXISTS '||t_dict_name;

EXECUTE 'DROP TABLE IF EXISTS '||t_sel_name;

END;
$BODY$
  LANGUAGE 'plpgsql' VOLATILE;

CREATE OR REPLACE FUNCTION optimize_fingerprint_step_two(t_schema text, t_name text, s_column text)
  RETURNS double precision AS
$BODY$
DECLARE i float8;
DECLARE tablesize integer;
DECLARE t_q_name text;
BEGIN

t_q_name := t_schema || '.q_' || t_name;

EXECUTE 'UPDATE '||t_q_name||' SET structure=mutabor('||s_column||')';

EXECUTE 'SELECT COUNT(1) FROM (SELECT count(fpstring('||s_column||')) as count FROM '||t_q_name||' GROUP BY fpstring('||s_column||')) as t WHERE t.count = 1' into i;

EXECUTE 'SELECT COUNT(1) FROM '||t_q_name into tablesize;

EXECUTE 'DROP TABLE IF EXISTS '||t_q_name;

RETURN i/tablesize;
END;
$BODY$
  LANGUAGE 'plpgsql' VOLATILE;

CREATE OR REPLACE FUNCTION baldi_tanimoto_hi(struct molecule, similarity double precision)
  RETURNS integer AS
$BODY$
BEGIN
RETURN ceil(nbits_set(fp2string(struct))/similarity)::integer;
END;
$BODY$
  LANGUAGE 'plpgsql' VOLATILE
  COST 100;

CREATE OR REPLACE FUNCTION baldi_tanimoto_lo(struct molecule, similarity double precision)
  RETURNS integer AS
$BODY$
BEGIN
RETURN floor(nbits_set(fp2string(struct))*similarity)::integer;
END;
$BODY$
  LANGUAGE 'plpgsql' VOLATILE
  COST 100;

CREATE OR REPLACE FUNCTION fpmaccsstring(molecule)
  RETURNS bit varying AS
'$libdir/libpgchem', 'pgchem_fp_MACCS'
  LANGUAGE C IMMUTABLE STRICT
  COST 1;

CREATE OR REPLACE FUNCTION tversky(molecule, molecule, double precision, double precision)
  RETURNS double precision AS
'$libdir/libpgchem', 'pgchem_tversky'
  LANGUAGE C VOLATILE STRICT
  COST 1;

CREATE OR REPLACE FUNCTION dice(prototype molecule, variant molecule)
  RETURNS double precision AS
$BODY$
BEGIN
RETURN tversky(prototype, variant, 0.5, 0.5);
END;
$BODY$
  LANGUAGE 'plpgsql' VOLATILE
  COST 100;

CREATE OR REPLACE FUNCTION spectrophore(molecule)
  RETURNS text AS
'$libdir/libpgchem', 'pgchem_spectrophore'
  LANGUAGE C IMMUTABLE STRICT
  COST 1;

CREATE OR REPLACE FUNCTION spectrophore_array(mol molecule)
  RETURNS double precision[] AS
$BODY$
BEGIN
RETURN string_to_array(spectrophore(mol),';')::float8[];
END;
$BODY$
  LANGUAGE plpgsql IMMUTABLE STRICT
  COST 100;

/*CREATE OR REPLACE FUNCTION reversestring(character varying)
  RETURNS character varying AS
'$libdir/libpgchem', 'reversestring'
  LANGUAGE C IMMUTABLE STRICT
  COST 1;

CREATE OR REPLACE FUNCTION fp3fps(struct molecule)
  RETURNS text AS
$BODY$
DECLARE fp varchar;
DECLARE fps text;
DECLARE len integer;
DECLARE i integer;
BEGIN
fp := substring(fpstring(struct) from 1025)::varchar;
fps:='';

FOR i IN 1..512 BY 8 LOOP

fps := fps || to_hex((reversestring(substring(fp from i for 8))::bit(8))::integer);

END LOOP;

RETURN fps;
END;
$BODY$
  LANGUAGE plpgsql IMMUTABLE STRICT
  COST 100;

CREATE OR REPLACE FUNCTION fp2fps(struct molecule)
  RETURNS text AS
$BODY$
DECLARE fp varchar;
DECLARE fps text;
DECLARE len integer;
DECLARE i integer;
BEGIN
len:=1024;
fp := substring(fpstring(struct) for len)::varchar;
fps:='';

FOR i IN 1..len BY 8 LOOP

fps := fps || to_hex((reversestring(substring(fp from i for 8))::bit(8))::integer);

END LOOP;

RETURN fps;
END;
$BODY$
  LANGUAGE plpgsql IMMUTABLE STRICT
  COST 100;*/

