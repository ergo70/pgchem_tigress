-- Table: reactions

-- DROP TABLE reactions;

CREATE TABLE reactions
(
  rxn bytea NOT NULL,
  rxn_id int4 NOT NULL DEFAULT nextval('chembank.rxntest_rxn_id_seq'::text),
  CONSTRAINT rxntest_pkey PRIMARY KEY (rxn_id)
) 
WITHOUT OIDS;
ALTER TABLE reactions OWNER TO chembank;






