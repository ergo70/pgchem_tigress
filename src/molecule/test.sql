insert into vc.moltest select iiid, molecule_to_new_molecule(molecule) from vc.moltable

truncate table vc.moltest

create index mol_idx on vc.moltest using gist(mol);

-- cluster mol_idx ON vc.moltest;

select count(*) from vc.moltest where  (select mol from vc.moltest where iiid=33) = mol

select count(*) from vc.moltest where 'c2ncc1ccccc1n2' <= mol and smartsfilter('[#6]c2nc([#7])c1ccccc1n2',mol)=true

select count(*) from vc.moltest where 'c2ccc1ccccc1c2' <= mol limit 50

select count(*) from vc.moltest where 'c3ccc2cc1ccccc1cc2c3' <= mol limit 50

select count(*) from vc.moltest where '[Cl-]' <= mol --LIMIT 500

select count(*) from vc.moltest where '[Cl]' <= mol and smartsfilter('[Cl-]',mol)=true

select * from vc.moltest where 'CCN' <= mol --limit 250

select count(*) from vc.moltest where smartsfilter('CCN',mol) = true

select count(*) from vc.moltest where hasbz(mol)=true

select count(*) from chemcollect.moltest where  smartsfilter('a',mol)=true

select mol @ (select mol from vc.moltest where iiid=17000)  from vc.moltest where iiid= 18000

select count(*) from vc.moltest where 'CCN' @ mol > .7

select exactmass(mol) from vc.moltest where iiid= 2

select molformula(mol) from vc.moltest where iiid= 2

select molhash(mol) from vc.moltest where iiid= 2

select molecule_to_smiles(mol) from vc.moltest where iiid=17000

"7bc09f01793e7b3c2b4cb08ffb6f8715"

