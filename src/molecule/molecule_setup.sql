DROP TYPE molecule CASCADE;
DROP TYPE molfp CASCADE;

CREATE TYPE molecule;
CREATE TYPE molfp;

CREATE FUNCTION molfp_in(cstring)
    RETURNS molfp
    AS 'libpgchem'
    LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION molfp_out(molfp)
    RETURNS cstring
   AS 'libpgchem'
    LANGUAGE C IMMUTABLE STRICT;

CREATE TYPE molfp (
   input = molfp_in,
   output = molfp_out,
   internallength = 192,
   storage = PLAIN
   );

CREATE FUNCTION molecule_in(cstring)
    RETURNS molecule
    AS 'libpgchem'
    LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION molecule_out(molecule)
    RETURNS cstring
   AS 'libpgchem'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION molecule_recv(internal)
   RETURNS molecule
  AS 'libpgchem'
   LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION molecule_send(molecule)
   RETURNS bytea
  AS 'libpgchem'
   LANGUAGE C IMMUTABLE STRICT;

CREATE TYPE molecule (
   input = molecule_in,
   output = molecule_out,
   --receive = molecule_recv,
   --send = molecule_send,
   internallength = VARIABLE,
   storage = EXTENDED
   );

CREATE FUNCTION molfp_compress(internal)
RETURNS internal
AS 'libpgchem'
LANGUAGE 'C';

CREATE FUNCTION molfp_decompress(internal)
RETURNS internal
AS 'libpgchem'
LANGUAGE 'C';

CREATE FUNCTION molfp_penalty(internal,internal,internal)
RETURNS internal
AS 'libpgchem'
LANGUAGE 'C' WITH (isstrict);

CREATE FUNCTION molfp_picksplit(internal,internal)
RETURNS internal
AS 'libpgchem'
LANGUAGE 'C' WITH (isstrict);

CREATE FUNCTION molfp_union(internal, internal)
RETURNS internal
AS 'libpgchem'
LANGUAGE 'C';

CREATE FUNCTION molfp_same(internal, internal, internal)
RETURNS internal
AS 'libpgchem'
LANGUAGE 'C';

CREATE FUNCTION molfp_consistent(internal,internal,int4)
RETURNS bool
AS 'libpgchem'
LANGUAGE 'C';

CREATE FUNCTION molecule_contained_in(molecule,molecule)
RETURNS bool
AS 'libpgchem'
LANGUAGE 'C' with (isstrict);

CREATE FUNCTION molecule_contains(molecule,molecule)
RETURNS bool
AS 'libpgchem'
LANGUAGE 'C' with (isstrict);

CREATE FUNCTION molecule_equals(molecule,molecule)
RETURNS bool
AS 'libpgchem'
LANGUAGE 'C' with (isstrict);

CREATE FUNCTION molecule_similarity(molecule,molecule)
RETURNS double precision
AS 'libpgchem'
LANGUAGE 'C' with (isstrict);

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
        OPERATOR        6       =        RECHECK,
        OPERATOR        7       >=       RECHECK,
        OPERATOR        8       <=       RECHECK,
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
'libpgchem', 'molecule_in_text'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION molecule_in(character varying)
  RETURNS molecule AS
'libpgchem', 'molecule_in_varchar'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION molecule_in(bytea)
  RETURNS molecule AS
'libpgchem', 'molecule_in_bytea'
  LANGUAGE 'c' IMMUTABLE STRICT;

CREATE CAST (text AS molecule)
  WITH FUNCTION molecule_in(text);

CREATE CAST (character varying AS molecule)
  WITH FUNCTION molecule_in(character varying);

CREATE CAST (bytea AS molecule)
  WITH FUNCTION molecule_in(bytea);



