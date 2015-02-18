CREATE OR REPLACE FUNCTION chembank_find_matching_molecules(sarg bytea, "mode" bpchar, row_limit integer)
  RETURNS SETOF integer AS
$BODY$ 

-- Arguments:
-- sarg - the structure to search with NOT NULL
-- mode - the search mode 'S' -> Substructure, 'E' -> Exact NOT NULL
-- row_limit - limit result to n rows, 0 means no limit

DECLARE curr_row record;
DECLARE mol_fp text;
DECLARE mol_fp2 text;
DECLARE tmp_sql text;
--DECLARE tmp_smarts text; //If you want to use SMARTS instead
DECLARE tmp_sarg bytea;
BEGIN

IF (mode NOT IN ('E','S')) THEN
	RAISE EXCEPTION 'Invalid search mode: %', mode;
END IF;

tmp_sql='';

	tmp_sql:='SELECT c.iiid FROM bioactives c, molkeys ms, molfingerprints mt WHERE c.iiid=ms.iiid AND c.iiid=mt.iiid AND ' || tmp_sql;
	tmp_sarg:=sarg;

mol_fp:=molkeys_long(tmp_sarg);
--mol_fp2='';

IF (mode='S' AND mol_fp='n_atoms:6;n_bonds:6;n_rings:1;n_C2:6;n_C:6;n_b2:1;n_bar:6;n_r6:1;n_rar:1;n_rbz:1;n_psg14:6;') THEN
tmp_sql:=tmp_sql || 'n_rbz!=0';
ELSE

mol_fp:=replace(mol_fp,';',' AND ');

IF (mode='S') THEN
        --tmp_smarts := molecule_to_smiles(tmp_sarg_acd2d_moltable); //If you want to use SMARTS instead
	mol_fp:=replace(mol_fp,':','>=');
        mol_fp2:=encode(fingerprint2(tmp_sarg),'hex');
        mol_fp2:=' substruct_screen_fingerprint(decode('''||mol_fp2||''',''hex''),mt.fp2bitmap)=true AND ';
	tmp_sql:=tmp_sql || mol_fp || mol_fp2 || 'match_substruct(decode('''||encode(tmp_sarg,'hex')||''',''hex''),c.molecule)=TRUE';
        --tmp_sql:=tmp_sql || mol_fp || mol_fp2 || 'match_substruct_smarts('''||tmp_smarts||''',c.molecule)=TRUE'; //If you want to use SMARTS instead
ELSE
	mol_fp:=replace(mol_fp,':','=');
	tmp_sql:=tmp_sql || mol_fp || 'match_exact(decode('''||encode(tmp_sarg,'hex')||''',''hex''),c.molecule)=TRUE';
END IF;


END IF;

IF (row_limit > 0) THEN
	tmp_sql:=tmp_sql || ' LIMIT ' || row_limit;
END IF;

FOR curr_row IN EXECUTE tmp_sql LOOP
	RETURN NEXT curr_row.iiid;
END LOOP;

RETURN;
END;$BODY$
  LANGUAGE 'plpgsql' STABLE STRICT;

