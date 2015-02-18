---------------------------------------------------------------------------
--
-- molfp.sql-
--    This file shows how to create a new user-defined type and how to
--    use this new type.
-- 
--
-- Portions Copyright (c) 1996-2006, PostgreSQL Global Development Group
-- Portions Copyright (c) 1994, Regents of the University of California
--
-- $PostgreSQL: pgsql/src/tutorial/molfp.source,v 1.20 2006/03/05 15:59:11 momjian Exp $
--
---------------------------------------------------------------------------

-----------------------------
-- Creating a new type:
--	We are going to create a new type called 'molfp' which represents
--	molfp numbers.
--	A user-defined type must have an input and an output function, and
--	optionally can have binary input and output functions.  All of these
--	are usually user-defined C functions.
-----------------------------

-- Assume the user defined functions are in _OBJWD_/molfp$DLSUFFIX
-- (we do not want to assume this is in the dynamic loader search path).
-- Look at $PWD/molfp.c for the source.  Note that we declare all of
-- them as STRICT, so we do not need to cope with NULL inputs in the
-- C code.  We also mark them IMMUTABLE, since they always return the
-- same outputs given the same inputs.

-- the input function 'molfp_in' takes a null-terminated string (the 
-- textual representation of the type) and turns it into the internal
-- (in memory) representation. You will get a message telling you 'molfp'
-- does not exist yet but that's okay.

CREATE FUNCTION molfp_in(cstring)
   RETURNS molfp
   AS '_OBJWD_/molfp'
   LANGUAGE C IMMUTABLE STRICT;

-- the output function 'molfp_out' takes the internal representation and
-- converts it into the textual representation.

CREATE FUNCTION molfp_out(molfp)
   RETURNS cstring
   AS '_OBJWD_/molfp'
   LANGUAGE C IMMUTABLE STRICT;

-- the binary input function 'molfp_recv' takes a StringInfo buffer
-- and turns its contents into the internal representation.

--CREATE FUNCTION molfp_recv(internal)
--   RETURNS molfp
--   AS '_OBJWD_/molfp'
--   LANGUAGE C IMMUTABLE STRICT;

-- the binary output function 'molfp_send' takes the internal representation
-- and converts it into a (hopefully) platform-independent bytea string.

--CREATE FUNCTION molfp_send(molfp)
--   RETURNS bytea
--   AS '_OBJWD_/molfp'
--   LANGUAGE C IMMUTABLE STRICT;


-- now, we can create the type. The internallength specifies the size of the
-- memory block required to hold the type (we need two 8-byte doubles).

CREATE TYPE molfp (
   internallength = 128, 
   input = molfp_in,
   output = molfp_out,
   alignment = integer
);


-----------------------------
-- Using the new type:
--	user-defined types can be used like ordinary built-in types.
-----------------------------

-- eg. we can use it in a table

CREATE TABLE test_molfp (
	a	molfp,
	b	molfp
);

-- data for user-defined types are just strings in the proper textual
-- representation. 

INSERT INTO test_molfp VALUES ('(1.0, 2.5)', '(4.2, 3.55 )');
INSERT INTO test_molfp VALUES ('(33.0, 51.4)', '(100.42, 93.55)');

SELECT * FROM test_molfp;

-----------------------------
-- Creating an operator for the new type:
--	Let's define an add operator for molfp types. Since POSTGRES
--	supports function overloading, we'll use + as the add operator.
--	(Operator names can be reused with different numbers and types of 
--	arguments.)
-----------------------------

-- first, define a function molfp_add (also in molfp.c)
--CREATE FUNCTION molfp_add(molfp, molfp)
--   RETURNS molfp
--   AS '_OBJWD_/molfp'
--   LANGUAGE C IMMUTABLE STRICT;

-- we can now define the operator. We show a binary operator here but you
-- can also define unary operators by omitting either of leftarg or rightarg.
--CREATE OPERATOR + ( 
--   leftarg = molfp,
--   rightarg = molfp,
--   procedure = molfp_add,
--   commutator = +
--);


--SELECT (a + b) AS c FROM test_molfp;

-- Occasionally, you may find it useful to cast the string to the desired
-- type explicitly. :: denotes a type cast.

--SELECT  a + '(1.0,1.0)'::molfp AS aa,
--        b + '(1.0,1.0)'::molfp AS bb
--   FROM test_molfp;


-----------------------------
-- Creating aggregate functions
--	you can also define aggregate functions. The syntax is somewhat
--	cryptic but the idea is to express the aggregate in terms of state
--	transition functions.
-----------------------------

CREATE AGGREGATE molfp_sum (
   sfunc = molfp_add,
   basetype = molfp,
   stype = molfp,
   initcond = '(0,0)'
);

SELECT molfp_sum(a) FROM test_molfp;


-----------------------------
-- Interfacing New Types with Indexes:
--	We cannot define a secondary index (eg. a B-tree) over the new type
--	yet. We need to create all the required operators and support
--      functions, then we can make the operator class.
-----------------------------

-- first, define the required operators
--CREATE FUNCTION molfp_abs_lt(molfp, molfp) RETURNS bool
--   AS '_OBJWD_/molfp' LANGUAGE C IMMUTABLE STRICT;
CREATE FUNCTION molfp_abs_le(molfp, molfp) RETURNS bool
   AS '_OBJWD_/molfp' LANGUAGE C IMMUTABLE STRICT;
CREATE FUNCTION molfp_abs_eq(molfp, molfp) RETURNS bool
   AS '_OBJWD_/molfp' LANGUAGE C IMMUTABLE STRICT;
CREATE FUNCTION molfp_abs_ge(molfp, molfp) RETURNS bool
   AS '_OBJWD_/molfp' LANGUAGE C IMMUTABLE STRICT;
--CREATE FUNCTION molfp_abs_gt(molfp, molfp) RETURNS bool
--   AS '_OBJWD_/molfp' LANGUAGE C IMMUTABLE STRICT;

--CREATE OPERATOR < (
--   leftarg = molfp, rightarg = molfp, procedure = molfp_abs_lt,
--   commutator = > , negator = >= ,
--   restrict = scalarltsel, join = scalarltjoinsel
--);
CREATE OPERATOR <= (
   leftarg = molfp, rightarg = molfp, procedure = molfp_abs_le,
   commutator = >= ,
   restrict = scalarltsel, join = scalarltjoinsel
);
CREATE OPERATOR = (
   leftarg = molfp, rightarg = molfp, procedure = molfp_abs_eq,
   commutator = = ,
   -- leave out negator since we didn't create <> operator
   -- negator = <> ,
   restrict = eqsel, join = eqjoinsel
);
CREATE OPERATOR >= (
   leftarg = molfp, rightarg = molfp, procedure = molfp_abs_ge,
   commutator = <= ,
   restrict = scalargtsel, join = scalargtjoinsel
);
--CREATE OPERATOR > (
--   leftarg = molfp, rightarg = molfp, procedure = molfp_abs_gt,
--   commutator = < , negator = <= ,
--   restrict = scalargtsel, join = scalargtjoinsel
--);

-- create the support function too
CREATE FUNCTION molfp_abs_cmp(molfp, molfp) RETURNS int4
   AS '_OBJWD_/molfp' LANGUAGE C IMMUTABLE STRICT;

-- now we can make the operator class
CREATE OPERATOR CLASS molfp_abs_ops
    DEFAULT FOR TYPE molfp USING btree AS
--        OPERATOR        1       < ,
        OPERATOR        2       <= ,
        OPERATOR        3       = ,
        OPERATOR        4       >= ,
--        OPERATOR        5       > ,
        FUNCTION        1       molfp_abs_cmp(molfp, molfp);


-- now, we can define a btree index on molfp types. First, let's populate
-- the table. Note that postgres needs many more tuples to start using the
-- btree index during selects.
INSERT INTO test_molfp VALUES ('(56.0,-22.5)', '(-43.2,-0.07)');
INSERT INTO test_molfp VALUES ('(-91.9,33.6)', '(8.6,3.0)');

CREATE INDEX test_cplx_ind ON test_molfp
   USING btree(a molfp_abs_ops);

SELECT * from test_molfp where a = '(56.0,-22.5)';
SELECT * from test_molfp where a < '(56.0,-22.5)';
SELECT * from test_molfp where a > '(56.0,-22.5)';


-- clean up the example
DROP TABLE test_molfp;
DROP TYPE molfp CASCADE;
