CREATE TABLE chembank.bioactives
(
  iiid serial NOT NULL,
  chembankid int4 NOT NULL,
  name text NOT NULL,
  compoundname text,
  molecule bytea NOT NULL,
  molformula text NOT NULL,
  molweight float4 NOT NULL,
  casno text,
  nscno int4,
  CONSTRAINT pk_bioactives PRIMARY KEY (iiid)
) WITH OIDS;
GRANT ALL ON TABLE chembank.bioactives TO chembank WITH GRANT OPTION;



