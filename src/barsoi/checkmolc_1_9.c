/*
   Modified to run w/o p2c libraries
   2007 Ernst-Georg Schmid
   Modified to run as shared library
   2007 Ernst-Georg Schmid
   pgchem@tuschehund.de
   
   compile with gcc without optimizations (except -march= -mcpu=) or it gets SLOWER!
   
   THIS VERSION HAS extended_molstat SWITCHED ON BY DEFAULT!
   
   This file is part of the xchem::tigress project.

   --- DUAL LICENSING ---
   
   If checkmol.c aka. pgchem::barsoi is compiled as a library, it is released under the terms of the
   lesser GNU General Public License. For the executable, the original GPL
   licensing applies!
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the lesser GNU General Public License as published by
   the Free Software Foundation version 2.1 of the License.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   lesser GNU General Public License for more details.
   
   --- DUAL LICENSING ---

   ORIGINAL HEADER:
   
checkmol/matchmol
Norbert Haider, University of Vienna, 2003-2007
norbert.haider@univie.ac.at

This software is published under the terms of the GNU General Public
License (GPL, see below). For a detailed description of this license,
see http://www.gnu.org/copyleft/gpl.html

If invoked as "checkmol", this program reads 2D and 3D molecular structure
files in various formats (Tripos Alchemy "mol", Tripos SYBYL "mol2", MDL "mol")
and describes the molecule either by its functional groups or by a set of
descriptors useful for database pre-screening (number of rings, sp2-hybridized
carbons, aromatic bonds, etc.; see definition of molstat_rec).

If invoked as "matchmol", the program reads two individual 2D or 3D molecular
structure files in various formats (see above) and checks if the first molecule
(the "needle") is a substructure of the second one (the "haystack").
"Haystack" can also be a MDL SD-file (containing multiple MOL files);
if invoked with "-" as file argument, both "needle" and "haystack" are
read as an SD-file from standard input, assuming the first entry in
the SDF to be the "needle"; output: entry number + ":F" (false) or ":T" (true)


Compile with fpc (Free Pascal, see http://www.freepascal.org), using
the -Sd or -S2 option (Delphi mode; IMPORTANT!)

example for compilation (with optimization) and installation:

fpc checkmol.pas -S2 -Xs -OG -Or -O2u -Op1

as "root", do the following:

cp checkmol /usr/local/bin
cd /usr/local/bin
ln checkmol matchmol


Version history

v0.1   extends "chkmol" utility (N. Haider, 2002), adds matching functionality;

v0.1a  minor bugfixes:
       stop ring search when max_rings is reached (function is_ringpath);
       fixed upper/lowercase of element symbol in MDL molfile output;
       added debug output for checkmol (-D option)

v0.2   added functionality to switch ring search from SAR (set of all rings) to
       SSR (set of small rings; see below) if number of rings exceeds max_rings
       (1024), e.g. with Buckminster-fullerenes (thanks to E.-G. Schmid for a hint);
       added -r command-line option to force ring search into SSR mode;
       as SSR, we regard the set of all rings with max. 10 ring members, and in
       which no ring is completely contained in another one

v0.2a  fixed a bug in function is_ringpath which could cause SAR ring search to
       overlook a few rings (e.g., the six 8-membered rings in cubane)

v0.2b  fixed unequal treatment of needle and haystack structures with respect
       to aromaticity check (trusting input file vs. full check)

v0.2c  modified the changes of v0.2b so that they affect only "matchmol" mode;
       added quick_match function

v0.2d  refined function bondtypes_OK so that it does not accept C=O fragments
       (and C=S, C=Se, C=Te) as equivalent to aromatic C-O (C-S, C-Se, C-Te);
       fixed function find_ndl_ref_atom to use only heavy atoms as "needle"
       reference atoms for atom-to-atom match; update function quick_match to
       handle queries with only one heavy atom (thanks to A. Barozza for a hint)

v0.2e  modified procedure get_molstat: adds 1 formal double bond for each aromatic
       ring to the "fingerprint" in order to allow an isolated double bond in the
       needle to be matched as a fragment of an aromatic ring
       ATTENTION: "fingerprints" (molecular statistics) generated with a previous
       version of checkmol(-x and -X option) must be updated for structure/
       substructure database searching with matchmaol.

v0.2f  modified procedures chk_ccx, chk_xccx, write_fg_text, write_fg_text_de,
       write_fg_code in order to report halogen derivatives other than alkyl halides,
       aryl halides and acyl halides (e.g., vinyl halides) and C=C compounds other
       than "true" alkenes, enols, enediols, enamines, etc. (e.g. vinyl halides);
       added "strict mode" (option -s) for matching: this uses a more thorough
       comparison of atom types and bond types (new functions atomtypes_OK_strict,
       bondtypes_OK_strict, several minor modifications in other subroutines).

v0.2g  changed procedure readinputfile (+ minor changes in main routine) in
       order to read very large SD input files without "not enough memory" error
       (the previous version attempted to read the entire SD file into a buffer,
       the new version reads only one molecule at once);
       fixed a minor bug in read_mol2file which led to undefined bond types;
       added definition of "debug" option as a compiler flag.

v0.2h  fixed a small bug which caused program hangs with (e.g.) some metal
       complexes: this was solved just by increasing the number of possible
       neighbor atoms from 8 to 16 (now defined in the constant max_neighbors).

v0.2i  introduced some more plausibility checking (number of directly connected
       atoms), result stored in variable mol_OK (set in count_neighbors);
       completely rewrote function matrix_OK which now can handle up to
       max_neighbors substituents per atom (previous versions were limited to
       4 substituents; larger numbers are required e.g. for ferrocenes and
       similar compounds).

v0.2j  new function raw_hetbond_count: ignores bond order (in contrast to function
       hetbond_count), this is better suitable for molecular statistics.
       ATTENTION: previously generated "fingerprints" in databases must be updated!
       improved recognition of non-conformant SD files (with missing "M  END" line)

v0.2k  changed quick_match routine to compare atom types only by element in
       order to avoid (rare) cases of non-matching identical input files;
       slightly changed error message for non-existant input files

v0.3   changed function update_Htotal in order to distinguish between 3-valent
       and 5-valent phosphorus (thanks to H. Feldman for this suggestion);
       added a table (array ringprop) to store ring sizes and aromaticity for
       faster lookup; changed aromaticity detection (chk_arom) to be fully
       independent of Kekulé structures in condensed ring systems; changed add_ring
       to store new rings in ascending order (with respect to ring size): this
       will cause the aromaticity detection to start with smaller rings;
       added additional calls to chk_arom when in SSR ring search mode (to ensure
       that all aromatic rings are found); speeded up function is_newring;
       added option "-M": this restores the behavior of previous versions,
       i.e. metal atoms will be accepted as ring members (which costs a lot
       of performance for coordinate compounds like ferrocenes but gives
       only little useful information); when used in databases: use the _same_
       mode ("-M" or not) _both_ for checkmol (creation of fingerprints) and matchmol;
       fixed ugly linebreaks in show_usage;

v0.3a  fixed a bug in read_charges (which was introduced in the v0.3 code cleanup)

v0.3b  fixed recognition of ketenes, CO2, CS2, CSO, phosphines, boronic esters

v0.3c  fixed recognition of hydrazines, hydroxylamines, thiocarboxylic acids,
       orthocarboxylic acids; slightly changed textual output for a few functional
       groups

v0.3d  added bond topology feature ("any", "ring", "chain", as defined by MDL specs,
       plus "always_any", "excess-ringcount", "exact-ringcount") to bond properties
       and implemented its use for substructure searches, either as specified in the
       query MDL molfile or by using "strict" mode; added ring_count property to
       "bonds" section of Checkmol-tweaked molfiles (using an unused field in the MDL
       molfile specs);
       added option for E/Z isomer search, either globally (-g option) or per
       double bond: bstereo_xyz property, encoded either by a leading "0" in the
       bond block of a checkmol-tweaked molfile or by using the "stereo care" flag
       as defined in the MDL file specs (both atoms of a double bond must carry
       this flag)

v0.3e  fixed wrong position of "stereo care" flag expected in input molfile;

v0.3f  added option for enantiospecific search (R/S isomer search); this can be
       turned on either globally (-G option or using the "chiral" flag in the
       "counts line" of the query MDL molfile) or per atom (using - or abusing -
       the "stereo care" flag); enantiomers are compared using the 3D coordinates
       of the substituents at a chiral center, so we do not have to rely on the
       presence of correct "up" and "down" bond types in the candidate structures;
       nevertheless, "pseudo-3D" structures (i.e. flat 2D structures with "up"
       and "down" bond notation) can be also used, even in mixed mode both as
       query structure and candidate structure; added support for "up" and "down"
       bond notation in MDL molfile output (if -m option is used and these bond
       types were present in the input file); refined function find_ndl_ref_atom;
       fixed Alchemy atom type misinterpretation (upper/lower case mismatch);

v0.3g  minor fixes: recognition of sp2 nitrogen in N=N structures (function
       update_atypes); extended E/Z geometry check to C=N and N=N bonds; accelerated
       exact search by initial comparison of C,O,N counts; added meaningful error
       message for input structures with 0 (zero) atoms;

v0.3h  added a match matrix clean-up step to function is_matching in order to
       remove "impossible" multiple matches prior to chirality check; added
       function ndl_maybe_chiral as a plausibility check; added support for other
       atoms than carbon as chiral centers; changed function is_cis from a
       distance-based cis/trans check into a torsion-based check;

v0.3i  minor fixes in functions max_neighbors and chk_ammon;
       revert ringsearch_mode to its original value for each new molecule
       (main routine)

v0.3j  various improvements in checkmol: added alkenyl (vinyl) and alkynyl
       residues as possible substituents for amides, esters, ethers, etc.;
       extended recognition of orthocarboxylic acid derivatives; refined
       recognition of oxohetarenes and iminohetareenes, ureas, nitroso
       compounds; fixed reading of "old-style" charges from MDL files;
       refined aromaticity check; fixed a bug in procedure chk_arom and
       renamed function is_exocyclic_methylene_C into
       find_exocyclic_methylene_C; added assignment of pointers to "nil"
       after calling "freemem" in procedure zap_molecule and zap_needle;
       matchmol: refined selection of needle reference atom; added "strict
       chirality check" mode (if -x -s and -G options are used): returns
       "false" if a chiral or pseudochiral atom is matched against its
       counterpart with undefined or "flat" geometry; added an alternative
       selection mechanism for the needle reference atom, based on the
       Morgan algorithm;

v0.3k  improved matching of nitro compounds (ionic vs. non-ionic
       representation); some finetuning in recognition of N-oxides,
       hydroxylamines, etc.; modified get_molstat in order to ignore
       charges in nitro groups, N-oxides, sulfoxides, azides, etc. in
       ionic representation; added function normalize_ionic_bonds;
       fixed a bug in is_alkenyl; fixed a bug (introduced in v0.3j) in
       the treatment of "old-style" MDL charges;

v0.3l  minor adjustments in quick_match (for unknown atom types,
       opposite charges); some performance improvements, e.g. in function
       path_pos, reduced extensive calls to is_heavyatom and is_metal
       by extending atom_rec with boolean fields 'heavy' and 'metal'
       (thanks to E.-G.Schmid); added molstat field n_rBz (number of
       benzene rings; deactivated by default); added -l command-line
       option to checkmol (lists all molstat codes with a short
       explanation); added graceful handling of NoStruct (i.e.,
       zero-atom) molecules

v0.3m  minor bug fixes (chk_imine), let normalize_ionic_bonds return "true"
       if some bond-order change was made; fixed incorrect implementation
       of the Morgan algorithm (thanks to H.Feldman), fixed is_matching in
       order to prevent wrong matches of larger rings on smaller ones; count
       halogens also as type "X" in molstat (as before v0.3l); added some more
       molstat descriptors for most elements (deactivated by default; for
       activation, uncomment the definition of "extended_molstat", see below);
       added support for the generation of binary fingerprints with output
       either as boolean values (-f option) or as decimal values (-F option),
       fingerprint patterns are supplied by the user (in SDF format, max. 62
       structures per file); added a version check ("tweaklevel",
       ringsearch_mode, opt_metalrings) for "tweaked" MDL molfiles

v0.3n  increased max_vringsize in SSR mode from 10 to 12 (now defined in
       ssr_vringsize); optionally changed ring statistics to ignore envelope
       rings (which completely contain another ring, "reduced SAR") in order to
       make molstat in SAR and SSR mode more compatible) - disabled by default;
       introduced a new molstat descriptor: n_br2p (number of bonds which
       belong to two or more rings); changed n_b1 counter to ignore bonds to
       explicit hydrogens; added procedure fix_ssr_ringcounts
       (see comment in the code); added bond property mdl_stereo for preservation
       of original value; slightly changed "get_molstat" and "update_atypes" in
       order to consolidate atom type for various N atoms; remember "chiral" flag
       status in molfile output (-m); modified chk_arom in order to accept rings
       as aromatic if all bonds are of type 'A'; several minor bug fixes

v0.3o  minor changes in update_atypes (nitrogen); changed the matching
       algorithm in order to take care of disconnected molecular graphs
       (e.g., salts); refined matching of atoms with formal charges


Credits:

Rami Jbara (rami_jbara@hotmail.com) designed the 8-character functional
group codes.

I am also very grateful to Ernst-Georg Schmid (Bayer Business Services,
Leverkusen/Germany) and to Howard Feldman (The Blueprint Initiative,
Toronto/Canada) for numerous bug reports and suggestions for improvements.



===============================================================================
DISCLAIMER
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software Foundation,
Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
===============================================================================
*/

/* $DEFINE debug                    uncomment if desired */
/* $DEFINE extended_molstat         uncomment if desired */
/* $DEFINE reduced_SAR              uncomment if desired */

/* uses
  SYSUTILS, MATH; */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "safemem.h"
#ifdef MAKE_SHARED_LIBRARY
#include "barsoi.h"
#endif

#ifdef MAKE_SHARED_LIBRARY
#define version         "0.3p C (pgchem::barsoi)"
#else
#define version         "0.3p C"
#endif

#undef REDUCED_SAR

/* How many recursions before aborting is_ringpath and reverting to SSR */
#define max_ringpath_recursion_depth   100000

/* How many recursions before aborting is_match and reverting to non-exhaustive match */
#define max_match_recursion_depth   100000

#define tweaklevel      2
    /* v0.3m; accept tweaks in MDL molfiles only if same level */
#define max_atoms       1024
#define max_bonds       1024
#define slack           8192
#define max_ringsize    128
#define max_rings       1024
#define max_fg          256
#define max_neighbors   16	/* new in v0.2h */

#define TAB             '\t'

#define max_matchpath_length  256
#define pmCheckMol      1001
#define pmMatchMol      1002

#define rs_sar          2001	/* ring search mode: SAR = set of all rings */
#define rs_ssr          2002
    /*                   SSR = set of small rings */

#define btopo_any       0	/* bond topology, new in v0.3d */
#define btopo_ring      1
#define btopo_chain     2
#define btopo_always_any  3	/* even in "strict mode" */
#define btopo_excess_rc  4
/* bond in candidate must be included in _more_ rings than
                        the matching bond in the query ==> specific search for
                        annulated systems */
#define btopo_exact_rc  5
    /* bond in query and candidate must have same ring count */

#define bstereo_any     0
    /* new in v0.3d, any E/Z isomer matches (for double bonds) */
#define bstereo_xyz     1
    /* E/Z match is checked by using XYZ coordinates of the atoms */
#define bstereo_up      11	/* new in v0.3f, flags for single bonds */
#define bstereo_down    16
#define bstereo_either    14	/* 0.3x */
#define bstereo_double_either    13	/* 0.3x */

#define fpf_boolean     3001
    /* v0.3m; format for fingerprint output (n:T, n:F) */
#define fpf_decimal     3002
    /* v0.3m; format for fingerprint output as decimal value of bitstring */
#define fp_blocksize    62
    /* v0.3m; 1 bit may be used for sign (e.g. by PHP), 1 bit for "exact hit" */
#define ssr_vringsize   12	/* v0.3n; max. ring size in SSR mode */

/* Definitions for functional groups: */
#define fg_cation       1
#define fg_anion        2
#define fg_carbonyl     3
#define fg_aldehyde     4
#define fg_ketone       5
#define fg_thiocarbonyl  6
#define fg_thioaldehyde  7
#define fg_thioketone   8
#define fg_imine        9
#define fg_hydrazone    10
#define fg_semicarbazone  11
#define fg_thiosemicarbazone  12
#define fg_oxime        13
#define fg_oxime_ether  14
#define fg_ketene       15
#define fg_ketene_acetal_deriv  16
#define fg_carbonyl_hydrate  17
#define fg_hemiacetal   18
#define fg_acetal       19
#define fg_hemiaminal   20
#define fg_aminal       21
#define fg_thiohemiaminal  22
#define fg_thioacetal   23
#define fg_enamine      24
#define fg_enol         25
#define fg_enolether    26
#define fg_hydroxy      27
#define fg_alcohol      28
#define fg_prim_alcohol  29
#define fg_sec_alcohol  30
#define fg_tert_alcohol  31
#define fg_1_2_diol     32
#define fg_1_2_aminoalcohol  33
#define fg_phenol       34
#define fg_1_2_diphenol  35
#define fg_enediol      36
#define fg_ether        37
#define fg_dialkylether  38
#define fg_alkylarylether  39
#define fg_diarylether  40
#define fg_thioether    41
#define fg_disulfide    42
#define fg_peroxide     43
#define fg_hydroperoxide  44
#define fg_hydrazine    45
#define fg_hydroxylamine  46
#define fg_amine        47
#define fg_prim_amine   48
#define fg_prim_aliph_amine  49
#define fg_prim_arom_amine  50
#define fg_sec_amine    51
#define fg_sec_aliph_amine  52
#define fg_sec_mixed_amine  53
#define fg_sec_arom_amine  54
#define fg_tert_amine   55
#define fg_tert_aliph_amine  56
#define fg_tert_mixed_amine  57
#define fg_tert_arom_amine  58
#define fg_quart_ammonium  59
#define fg_n_oxide      60
#define fg_halogen_deriv  61
#define fg_alkyl_halide  62
#define fg_alkyl_fluoride  63
#define fg_alkyl_chloride  64
#define fg_alkyl_bromide  65
#define fg_alkyl_iodide  66
#define fg_aryl_halide  67
#define fg_aryl_fluoride  68
#define fg_aryl_chloride  69
#define fg_aryl_bromide  70
#define fg_aryl_iodide  71
#define fg_organometallic  72
#define fg_organolithium  73
#define fg_organomagnesium  74
#define fg_carboxylic_acid_deriv  75
#define fg_carboxylic_acid  76
#define fg_carboxylic_acid_salt  77
#define fg_carboxylic_acid_ester  78
#define fg_lactone      79
#define fg_carboxylic_acid_amide  80
#define fg_carboxylic_acid_prim_amide  81
#define fg_carboxylic_acid_sec_amide  82
#define fg_carboxylic_acid_tert_amide  83
#define fg_lactam       84
#define fg_carboxylic_acid_hydrazide  85
#define fg_carboxylic_acid_azide  86
#define fg_hydroxamic_acid  87
#define fg_carboxylic_acid_amidine  88
#define fg_carboxylic_acid_amidrazone  89
#define fg_nitrile      90
#define fg_acyl_halide  91
#define fg_acyl_fluoride  92
#define fg_acyl_chloride  93
#define fg_acyl_bromide  94
#define fg_acyl_iodide  95
#define fg_acyl_cyanide  96
#define fg_imido_ester  97
#define fg_imidoyl_halide  98
#define fg_thiocarboxylic_acid_deriv  99
#define fg_thiocarboxylic_acid  100
#define fg_thiocarboxylic_acid_ester  101
#define fg_thiolactone  102
#define fg_thiocarboxylic_acid_amide  103
#define fg_thiolactam   104
#define fg_imido_thioester  105
#define fg_oxohetarene  106
#define fg_thioxohetarene  107
#define fg_iminohetarene  108
#define fg_orthocarboxylic_acid_deriv  109
#define fg_carboxylic_acid_orthoester  110
#define fg_carboxylic_acid_amide_acetal  111
#define fg_carboxylic_acid_anhydride  112
#define fg_carboxylic_acid_imide  113
#define fg_carboxylic_acid_unsubst_imide  114
#define fg_carboxylic_acid_subst_imide  115
#define fg_co2_deriv    116
#define fg_carbonic_acid_deriv  117
#define fg_carbonic_acid_monoester  118
#define fg_carbonic_acid_diester  119
#define fg_carbonic_acid_ester_halide  120
#define fg_thiocarbonic_acid_deriv  121
#define fg_thiocarbonic_acid_monoester  122
#define fg_thiocarbonic_acid_diester  123
#define fg_thiocarbonic_acid_ester_halide  124
#define fg_carbamic_acid_deriv  125
#define fg_carbamic_acid  126
#define fg_carbamic_acid_ester  127
#define fg_carbamic_acid_halide  128
#define fg_thiocarbamic_acid_deriv  129
#define fg_thiocarbamic_acid  130
#define fg_thiocarbamic_acid_ester  131
#define fg_thiocarbamic_acid_halide  132
#define fg_urea         133
#define fg_isourea      134
#define fg_thiourea     135
#define fg_isothiourea  136
#define fg_guanidine    137
#define fg_semicarbazide  138
#define fg_thiosemicarbazide  139
#define fg_azide        140
#define fg_azo_compound  141
#define fg_diazonium_salt  142
#define fg_isonitrile   143
#define fg_cyanate      144
#define fg_isocyanate   145
#define fg_thiocyanate  146
#define fg_isothiocyanate  147
#define fg_carbodiimide  148
#define fg_nitroso_compound  149
#define fg_nitro_compound  150
#define fg_nitrite      151
#define fg_nitrate      152
#define fg_sulfuric_acid_deriv  153
#define fg_sulfuric_acid  154
#define fg_sulfuric_acid_monoester  155
#define fg_sulfuric_acid_diester  156
#define fg_sulfuric_acid_amide_ester  157
#define fg_sulfuric_acid_amide  158
#define fg_sulfuric_acid_diamide  159
#define fg_sulfuryl_halide  160
#define fg_sulfonic_acid_deriv  161
#define fg_sulfonic_acid  162
#define fg_sulfonic_acid_ester  163
#define fg_sulfonamide  164
#define fg_sulfonyl_halide  165
#define fg_sulfone      166
#define fg_sulfoxide    167
#define fg_sulfinic_acid_deriv  168
#define fg_sulfinic_acid  169
#define fg_sulfinic_acid_ester  170
#define fg_sulfinic_acid_halide  171
#define fg_sulfinic_acid_amide  172
#define fg_sulfenic_acid_deriv  173
#define fg_sulfenic_acid  174
#define fg_sulfenic_acid_ester  175
#define fg_sulfenic_acid_halide  176
#define fg_sulfenic_acid_amide  177
#define fg_thiol        178
#define fg_alkylthiol   179
#define fg_arylthiol    180
#define fg_phosphoric_acid_deriv  181
#define fg_phosphoric_acid  182
#define fg_phosphoric_acid_ester  183
#define fg_phosphoric_acid_halide  184
#define fg_phosphoric_acid_amide  185
#define fg_thiophosphoric_acid_deriv  186
#define fg_thiophosphoric_acid  187
#define fg_thiophosphoric_acid_ester  188
#define fg_thiophosphoric_acid_halide  189
#define fg_thiophosphoric_acid_amide  190
#define fg_phosphonic_acid_deriv  191
#define fg_phosphonic_acid  192
#define fg_phosphonic_acid_ester  193
#define fg_phosphine    194
#define fg_phosphinoxide  195
#define fg_boronic_acid_deriv  196
#define fg_boronic_acid  197
#define fg_boronic_acid_ester  198
#define fg_alkene       199
#define fg_alkyne       200
#define fg_aromatic     201
#define fg_heterocycle  202
#define fg_alpha_aminoacid  203
#define fg_alpha_hydroxyacid  204

typedef enum
{ false = 0, true } boolean;


typedef char str2[3];

typedef char str3[4];

typedef char str4[5];

typedef char str5[6];

typedef char str8[9];

typedef struct atom_rec
{
  str2 element;
  str3 atype;
  float x, y, z;
  int formal_charge;
  float real_charge;
  int Hexp;			/* explicit H count */
  int Htot;			/* total H count */
  int neighbor_count, ring_count;
  boolean arom, stereo_care;	/* new in v0.3d */
  boolean q_arom;		// v0.3p potentially aromatic in a query structure
  boolean heavy;		/* new in v0.3l */
  boolean metal;		/* new in v0.3l */
  int nvalences;		/* new in v0.3m */
  boolean tag;			/* new in v0.3o */
  int nucleon_number;		/* 0.3.x */
  int radical_type;		/* 0.3.x */
} atom_rec;

typedef struct bond_rec
{
  int a1, a2;
  char btype;
  int ring_count;
  boolean arom;
  boolean q_arom;		// v0.3p potentially aromatic in a query structure
  int topo;			/* new in v0.3d, see MDL file description */
  int stereo;			/* new in v0.3d */
  int mdl_stereo;		/* new in v0.3n */
} bond_rec;

typedef int ringpath_type[max_ringsize];
typedef int matchpath_type[max_matchpath_length];

typedef atom_rec atomlist[max_atoms];
typedef bond_rec bondlist[max_bonds];
typedef ringpath_type ringlist[max_rings];
typedef int neighbor_rec[max_neighbors];	/* new in v0.2h */
typedef boolean fglist[max_fg];

typedef char molbuftype[max_atoms + max_bonds + slack][256];

typedef boolean matchmatrix[max_neighbors][max_neighbors];
    /* new in v0.2i */

typedef struct molstat_rec
{
  int n_QA, n_QB, n_chg;	/* number of query atoms, query bonds, charges */
  int n_C1, n_C2, n_C;
  /* number of sp, sp2 hybridized, and total no. of carbons */
  int n_CHB1p, n_CHB2p, n_CHB3p, n_CHB4;
  /* number of C atoms with at least 1, 2, 3 hetero bonds */
  int n_O2, n_O3;		/* number of sp2 and sp3 oxygens */
  int n_N1, n_N2, n_N3;		/* number of sp, sp2, and sp3 nitrogens */
  int n_S, n_SeTe;
  /* number of sulfur atoms and selenium or tellurium atoms */
  int n_F, n_Cl, n_Br, n_I;
  /* number of fluorine, chlorine, bromine, iodine atoms */
  int n_P, n_B;			/* number of phosphorus and boron atoms */
  int n_Met, n_X;
  /* number of metal and "other" atoms (not listed elsewhere); v0.3l */
  int n_b1, n_b2, n_b3, n_bar;
  /* number single, double, triple, and aromatic bonds */
  int n_C1O, n_C2O, n_CN, n_XY;
  /* number of C-O single bonds, C=O double bonds, CN bonds (any type), hetero/hetero bonds */
  int n_r3, n_r4, n_r5, n_r6, n_r7, n_r8;
  /* number of 3-, 4-, 5-, 6-, 7-, and 8-membered rings */
  int n_r9, n_r10, n_r11, n_r12, n_r13p;
  /* number of 9-, 10-, 11-, 12-, and 13plus-membered rings */
  int n_rN, n_rN1, n_rN2, n_rN3p;
  /* number of rings containing N (any number), 1 N, 2 N, and 3 N or more */
  int n_rO, n_rO1, n_rO2p;
  /* number of rings containing O (any number), 1 O, and 2 O or more */
  int n_rS, n_rX, n_rAr, n_rBz;
  /* number of rings containing S (any number), any heteroatom (any number),  */
  /* number of aromatic rings, number of benzene rings */
  int n_br2p;			/* number of bonds belonging to more than one ring (v0.3n) */
/* p2c: checkmol.pas, line 590:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
  /*$IFDEF extended_molstat */
  int n_psg01, n_psg02, n_psg13, n_psg14;
  /* number of atoms belonging to elements of group 1, 2, etc.  */
  int n_psg15, n_psg16, n_psg17, n_psg18;
  /* number of atoms belonging to elements of group 15, 16, etc.  */
  int n_pstm, n_psla;
  /* number of transition metals, lanthanides/actinides */
  int n_iso, n_rad;
  /* number of isotopes, radicals */
  /*$ENDIF */
} molstat_rec;

typedef struct ringprop_rec
{
  /* new in v0.3 */
  int size;
  boolean arom, envelope;
} ringprop_rec;

typedef ringprop_rec ringprop_type[max_rings];

typedef struct p_3d
{
  /* new in v0.3d */
  double x, y, z;
} p_3d;

typedef int chirpath_type[4];	/* new in v0.3f */

#ifdef MAKE_SHARED_LIBRARY
static boolean yet_initialized = false;
#endif

typedef struct connval_rec
{
  /* new in v0.3j */
  int def;			/* better as longint for large molecules? */
  int tmp;
} connval_rec;

typedef connval_rec connval_type[max_atoms];	/* new in v0.3j  */


static int progmode;
static char progname[256];

#ifndef MAKE_SHARED_LIBRARY
static int i;			/* general purpose index */
#endif

static int li;
static boolean opt_none, opt_verbose, opt_text, opt_text_de, opt_code,
  opt_bin, opt_bitstring, opt_stdin, opt_exact, opt_debug,
  opt_molout, opt_molstat, opt_molstat_X, opt_xmdlout, opt_strict;
    /* new in v0.2f */
static boolean opt_metalrings;	/* new in v0.3 */
static boolean opt_geom;	/* new in v0.3d */
static boolean opt_chiral;	/* new in v0.3f */
static boolean opt_iso;		/* new in v0.3x */
static boolean opt_chg;		/* new in v0.3x */
static boolean opt_rad;		/* new in v0.3x */
static int opt_rs;		/* new in v0.3i */

#ifdef MAKE_SHARED_LIBRARY
static int opt_rs_dll = RPA_DEFAULT;
#endif

static boolean opt_fp;		/* new in v0.3m */
static int fpformat;		/* new in v0.3m */
static char filetype[256];
/*molfile : text; */
static char molfilename[256];
static char ndl_molfilename[256];
static char molname[256];
static char ndl_molname[256];
static char tmp_molname[256];	/* v0.3m */
static char molcomment[256];
static int n_atoms, n_bonds, n_rings;
    /* the number of rings we determined ourselves */
static int n_countablerings;
    /* v0.3n; number of rings which are not envelope rings */
static int n_cmrings;
    /* the number of rings we read from a (CheckMol-tweaked) MDL molfile */
//static int n_charges;         /* number of charges */
static int n_heavyatoms, n_heavybonds, ndl_n_atoms, ndl_n_bonds, ndl_n_rings,
  ndl_n_heavyatoms, ndl_n_heavybonds;
/*cm_mdlmolfile  : boolean; */
static boolean found_arominfo, found_querymol, ndl_querymol;
static int tmp_n_atoms;		/* v0.3m */
static int tmp_n_bonds;		/* v0.3m */
static int tmp_n_rings;		/* v0.3m */
static int tmp_n_heavyatoms;	/* v0.3m */
static int tmp_n_heavybonds;	/* v0.3m */

static atom_rec *atom;
static bond_rec *bond;
static ringpath_type *ring;
static ringprop_rec *ringprop;	/* new in v0.3 */

static atom_rec *ndl_atom;
static bond_rec *ndl_bond;
static ringpath_type *ndl_ring;
static ringprop_rec *ndl_ringprop;	/* new in v0.3 */
static int ndl_ref_atom;	/* since v0.3j as a global variable */
static atom_rec *tmp_atom;	/* v0.3m */
static bond_rec *tmp_bond;	/* v0.3m */
static ringpath_type *tmp_ring;	/* v0.3m */
static ringprop_rec *tmp_ringprop;	/* v0.3m */

static boolean matchresult, matchsummary;	/* v0.3o */
static matchpath_type ndl_matchpath, hst_matchpath;

static fglist fg;
static str4 atomtype;
static str3 newatomtype;

static char (*molbuf)[256];
static int molbufindex;

static boolean mol_in_queue;
static int mol_count;

static molstat_rec molstat, ndl_molstat, tmp_molstat;	/* v0.3m */

static int ringsearch_mode, max_vringsize;	/* for SSR ring search */

static FILE *rfile;
static boolean rfile_is_open, mol_OK;	/* new in v0.2i */

static int n_ar;		/* new in v0.3 */
static int prev_n_ar;		/* new in v0.3 */
static boolean ez_search;	/* new in v0.3d */
static boolean rs_search;	/* new in v0.3f */
static boolean ez_flag;		/* new in v0.3f */
static boolean chir_flag;	/* new in v0.3f */
static boolean rs_strict;	/* new in v0.3j */

static int n_Ctot, n_Otot, n_Ntot;	/* new in v0.3g */
static int ndl_n_Ctot, ndl_n_Otot, ndl_n_Ntot;	/* new in v0.3g */
static int tmp_n_Ctot, tmp_n_Otot, tmp_n_Ntot;	/* new in v0.3m */
static boolean ether_generic;	/* v0.3j */
static boolean amine_generic;	/* v0.3j   */
static boolean hydroxy_generic;	/* v0.3j */
static connval_rec *cv;		/* new in v0.3j */
static long long fpdecimal;	/* v0.3m */
static long long fpincrement;
static int fpindex;
static boolean fp_exacthit, fp_exactblock;
static int tmfcode;
    /* v0.3m; version number for tweaked MDL molfiles (tweaklevel) */
static boolean tmfmismatch;
    /* v0.3m; rely on tweaked MDL molfiles only if same level */
static boolean auto_ssr;
    /* v0.3n; indicates that SAR -> SSR fallback has happened */
static int recursion_depth;
    /* ====================emulated pascal functions============================ */

static boolean
file_exists (const char *fileName)
{
  struct stat filestat;

  return stat (fileName, &filestat) == 0 ? true : false;
}

static void
lblank (int cols, char *nstr)
{
  /* inserts leading blanks up to a given number of total characters */
  char STR1[256];

  if (strlen (nstr) >= cols)
    return;
  while (strlen (nstr) < cols)
    sprintf (nstr, " %s", strcpy (STR1, nstr));
}

static inline double
radtodeg (double rads)
{
  return (180.0 / M_PI) * rads;
}

static void
strdelete (char *s, int pos, int len)
{
  int slen;

  if (--pos < 0)
    return;
  slen = strlen (s) - pos;
  if (slen <= 0)
    return;
  s += pos;
  if (slen <= len)
    {
      *s = 0;
      return;
    }
  while ((*s = s[len]))
    s++;
}

static int
strpos2 (char *s, char *pat, int pos)
{
  char *cp, ch;
  int slen;

  if (--pos < 0)
    return 0;
  slen = strlen (s) - pos;
  cp = s + pos;
  if (!(ch = *pat++))
    return 0;
  pos = strlen (pat);
  slen -= pos;
  while (--slen >= 0)
    {
      if (*cp++ == ch && !strncmp (cp, pat, pos))
	return cp - s;
    }
  return 0;
}

static char *
strsub (char *ret, char *s, int pos, int len)
{
  char *s2;

  if (--pos < 0 || len <= 0)
    {
      *ret = 0;
      return ret;
    }
  while (pos > 0)
    {
      if (!*s++)
	{
	  *ret = 0;
	  return ret;
	}
      pos--;
    }
  s2 = ret;
  while (--len >= 0)
    {
      if (!(*s2++ = *s++))
	return ret;
    }
  *s2 = 0;
  return ret;
}


/*============================= auxiliary functions & procedures */

inline static void
all_lowercase (char *astring)
{
  int i;
  int l = strlen (astring);

  for (i = 0; i < l; i++)
    {
      astring[i] = tolower (astring[i]);
    }
}

static void
init_globals ()
{
  int i;

  opt_verbose = false;
  opt_debug = false;
  opt_exact = false;
  opt_stdin = false;
  opt_text = false;
  opt_code = false;
  opt_bin = false;
  opt_bitstring = false;
  opt_molout = false;
  opt_molstat = false;
  opt_molstat_X = false;
  opt_xmdlout = false;
  opt_strict = false;		/* new in v0.2f */
  opt_metalrings = false;	/* new in v0.3 */
  opt_geom = false;		/* new in v0.3d */
  opt_chiral = false;		/* new in v0.3f */
  opt_fp = false;		/* new in v0.3m */
  opt_iso = false;		/* new in v0.3x */
  opt_chg = false;		/* new in v0.3x */
  opt_rad = false;		/* new in v0.3x */
  /*cm_mdlmolfile   := false; */
  found_arominfo = false;
  found_querymol = false;
  ndl_querymol = false;
  opt_rs = rs_sar;		/* v0.3i */
  /*ringsearch_mode := rs_sar; */
  rfile_is_open = false;	/* new in v0.2g */
  ez_search = false;		/* new in v0.3d */
  rs_search = false;		/* new in v0.3f */
  ez_flag = false;		/* new in v0.3f */
  chir_flag = false;		/* new in v0.3f */
  rs_strict = false;		/* new in v0.3j */
  n_Ctot = 0;
  n_Otot = 0;
  n_Ntot = 0;			/* new in v0.3g */
  ndl_n_Ctot = 0;
  ndl_n_Otot = 0;
  ndl_n_Ntot = 0;		/* new in v0.3g */
  //for (i = 0; i < max_fg; i++)
  //  fg[i] = false;
  memset (fg, 0, sizeof (fglist));
  /*try */
  molbuf = safe_malloc (sizeof (molbuftype));
  /*except
     on e:Eoutofmemory do
     begin
     writeln('Not enough memory');
     halt(4);
     end;
     end; */
  ether_generic = false;	/* v0.3j */
  amine_generic = false;	/* v0.3j */
  hydroxy_generic = false;	/* v0.3j */
  fpformat = fpf_decimal;	/* v0.3m */
  fpindex = 0;			/* v0.3m */
  fp_exacthit = false;		/* v0.3m */
  fp_exactblock = false;	/* v0.3m */
  tmfcode = 0;			/* v0.3m */
  tmfmismatch = false;		/* v0.3m */
  auto_ssr = false;		/* v0.3n */
  recursion_depth = 0;
}


static inline void
init_molstat (mstat)
     molstat_rec *mstat;
{
  /*
     with mstat do
     begin
     n_QA := 0; n_QB := 0; n_chg := 0;
     n_C1 := 0; n_C2 := 0; n_C  := 0;
     n_CHB1p := 0; n_CHB2p := 0; n_CHB3p := 0; n_CHB4 := 0;
     n_O2 := 0; n_O3  := 0;
     n_N1 := 0; n_N2 := 0; n_N3 := 0;
     n_S := 0; n_SeTe := 0;
     n_F := 0; n_Cl := 0; n_Br := 0; n_I := 0;
     n_P := 0; n_B := 0;
     n_Met := 0; n_X := 0;
     n_b1 := 0; n_b2 := 0; n_b3 := 0; n_bar := 0;
     n_C1O := 0; n_C2O := 0; n_CN := 0; n_XY := 0;
     n_r3 := 0; n_r4 := 0; n_r5 := 0; n_r6 := 0; n_r7 := 0;
     n_r8 := 0; n_r9 := 0; n_r10 := 0; n_r11 := 0; n_r12 := 0; n_r13p := 0;
     n_rN := 0; n_rN1 := 0; n_rN2 := 0; n_rN3p := 0;
     n_rO := 0; n_rO1 := 0; n_rO2p := 0;
     n_rS := 0; n_rX := 0;
     n_rAr := 0; n_rBz := 0; n_br2p := 0;
     end;
   */
  memset (mstat, 0, sizeof (molstat_rec));	/* v0.3k */
}


static void
debugoutput (dstr)
     char *dstr;
{
  if (opt_debug)
    printf ("%s\n", dstr);
}


static void
left_trim (trimstr)
     char *trimstr;
{
  while (*trimstr != '\0' && (trimstr[0] == ' ' || trimstr[0] == TAB))
    strdelete (trimstr, 1, 1);
}


static int
left_int (trimstr)
     char *trimstr;
{
  char numstr[256];
  char auxstr[256];
  int auxint = 0;
  int code;
  char STR1[256];

  strcpy (numstr, "-+0123456789");
  *auxstr = '\0';
  while (*trimstr != '\0' && (trimstr[0] == ' ' || trimstr[0] == TAB))
    strdelete (trimstr, 1, 1);
  while ((*trimstr != '\0') &&
	 (strpos2 (numstr, (sprintf (STR1, "%c", trimstr[0]), STR1), 1) > 0))
    {
      sprintf (auxstr + strlen (auxstr), "%c", trimstr[0]);
      strdelete (trimstr, 1, 1);
    }
  code = (sscanf (auxstr, "%ld", &auxint) == 0);
  return auxint;
}

static int
path_pos (int id, int *a_path)
{
  int i = 0;

  for (i = 0; i < max_ringsize; i++)
    {
      if (*(a_path++) == id)
	{
	  return ++i;
	}
    }
  return 0;
}

static int
path_length (int *a_path)
{
  if ((a_path[max_ringsize - 1] != 0) && (path_pos (0, a_path) == 0))
    return max_ringsize;
  else
    return (path_pos (0, a_path) - 1);
}

static int
get_bond (int ba1, int ba2)
{
  int i;
  int b_id = 0;
  int FORLIM;

  if (n_bonds <= 0)
    return b_id;
  FORLIM = n_bonds;
  for (i = 1; i <= FORLIM; i++)
    {
      if (bond[i - 1].a1 == ba1 && bond[i - 1].a2 == ba2 ||
	  bond[i - 1].a1 == ba2 && bond[i - 1].a2 == ba1)
	b_id = i;
    }
  return b_id;
}

static void
clear_atom_tags ()
{
  int i, FORLIM;

  if (n_atoms > 0)
    {
      FORLIM = n_atoms;
      for (i = 0; i < FORLIM; i++)
	atom[i].tag = false;
    }
}

static void
set_atom_tags ()
{
  int i, FORLIM;

  if (n_atoms > 0)
    {
      FORLIM = n_atoms;
      for (i = 0; i < FORLIM; i++)
	atom[i].tag = true;
    }
}

static void
order_ringpath (int *r_path)
{
  /* order should be: array starts with atom of lowest number, followed by neighbor atom with lower number */
  int i, pl, a_ref, a_left, a_right, a_tmp;

  pl = path_length (r_path);
  if (pl < 3)
    return;
  a_ref = n_atoms;
  /* start with highest possible value for an atom number */
  for (i = 0; i < pl; i++)
    {
      if (r_path[i] < a_ref)	/* find the minimum value ==> reference atom */
	a_ref = r_path[i];
    }
  if (a_ref < 1)		/* just to be sure */
    return;
  if (path_pos (a_ref, r_path) < pl)
    a_right = r_path[path_pos (a_ref, r_path)];
  else
    a_right = r_path[0];
  if (path_pos (a_ref, r_path) > 1)
    a_left = r_path[path_pos (a_ref, r_path) - 2];
  else
    a_left = r_path[pl - 1];
  if (a_right == a_left)	/* should never happen */
    return;
  if (a_right < a_left)
    {
      /* correct ring numbering direction, only shift of the reference atom to the left end required */
      while (path_pos (a_ref, r_path) > 1)
	{
	  a_tmp = r_path[0];
	  for (i = 1; i < pl; i++)
	    r_path[i - 1] = r_path[i];
	  r_path[pl - 1] = a_tmp;
	}
      return;
    }
  while (path_pos (a_ref, r_path) < pl)
    {
      /* step one: create "mirrored" ring path with reference atom at right end */
      a_tmp = r_path[pl - 1];
      for (i = pl; i >= 2; i--)
	r_path[i - 1] = r_path[i - 2];
      r_path[0] = a_tmp;
    }
  for (i = 1; i <= pl / 2; i++)
    {				/* one more mirroring */
      a_tmp = r_path[i - 1];
      r_path[i - 1] = r_path[pl - i];
      r_path[pl - i] = a_tmp;
      /* wrong ring numbering direction, two steps required */
    }
}

static void
clear_ndl_atom_tags ()
{
  int i;

  if (ndl_n_atoms > 0)
    {

      for (i = 0; i < ndl_n_atoms; i++)
	ndl_atom[i].tag = false;
    }
}


static void
set_ndl_atom_tags ()
{
  int i;

  if (ndl_n_atoms > 0)
    {

      for (i = 0; i < ndl_n_atoms; i++)
	ndl_atom[i].tag = true;
    }
}


static int
count_tagged_ndl_heavyatoms ()
{
  int i;
  int n = 0;

  if (ndl_n_atoms < 1)
    return n;

  for (i = 0; i < ndl_n_atoms; i++)
    {
      if (ndl_atom[i].heavy && ndl_atom[i].tag)
	n++;
    }
  return n;
}



/*============================= geometry functions ========================== */

static double
dist3d (p1, p2)
     p_3d p1, p2;
{
  double res, TEMP, TEMP1, TEMP2;

  TEMP = p1.x - p2.x;
  TEMP1 = p1.y - p2.y;
  TEMP2 = p1.z - p2.z;
  res = sqrt (TEMP * TEMP + TEMP1 * TEMP1 + TEMP2 * TEMP2);
  return res;
}


/*
function is_cis(p1,p2,p3,p4:p_3d):boolean;  (* new in v0.3d
var                         (* just a simple, distance-based estimation
  total_dist  : double;     (* instead of calculating the dihedral angle
  direct_dist : double;
  res         : boolean;
begin
  res := false;
  total_dist  := dist3d(p1,p2) + dist3d(p2,p3) + dist3d(p3,p4);
  direct_dist := dist3d(p1,p4);
  if (direct_dist < 0.78 * total_dist) then res := true;  (* cutoff value of 0.78 was
  is_cis := res;                                          (* experimentally determined
end;
*/
/* function is_cis was replaced by a new one in v0.3h */


static p_3d
subtract_3d (p1, p2)
     p_3d p1, p2;
{
  p_3d p;

  p.x = p1.x - p2.x;
  p.y = p1.y - p2.y;
  p.z = p1.z - p2.z;
  return p;
}


static p_3d
add_3d (p1, p2)
     p_3d p1, p2;
{
  p_3d p;

  p.x = p1.x + p2.x;
  p.y = p1.y + p2.y;
  p.z = p1.z + p2.z;
  return p;
}


static void
vec2origin (p1, p2)
     p_3d *p1, *p2;
{
  p_3d p;

  p = subtract_3d (*p2, *p1);
  *p2 = p;
  p1->x = 0.0;
  p1->y = 0.0;
  p1->z = 0.0;
}


static double
scalar_prod (p1, p2, p3)
     p_3d p1, p2, p3;
{
  p_3d p;
  double res;

  p = subtract_3d (p2, p1);
  p2 = p;
  p = subtract_3d (p3, p1);
  p3 = p;
  p1.x = 0.0;
  p1.y = 0.0;
  p1.z = 0.0;
  res = p2.x * p3.x + p2.y * p3.y + p2.z * p3.z;
  return res;
}


static p_3d
cross_prod (p1, p2, p3)
     p_3d p1, p2, p3;
{
  p_3d p, orig_p1;

  orig_p1 = p1;
  p = subtract_3d (p2, p1);
  p2 = p;
  p = subtract_3d (p3, p1);
  p3 = p;
  p.x = p2.y * p3.z - p2.z * p3.y;
  p.y = p2.z * p3.x - p2.x * p3.z;
  p.z = p2.x * p3.y - p2.y * p3.x;
  return (add_3d (orig_p1, p));
}


static double
angle_3d (p1, p2, p3)
     p_3d p1, p2, p3;
{
  p_3d lp1, lp2, lp3, p;
  double res = 0.0;
  double magn_1, magn_2, cos_phi;

  lp1 = p1;
  lp2 = p2;
  lp3 = p3;
  p = subtract_3d (lp2, lp1);
  lp2 = p;
  p = subtract_3d (lp3, lp1);
  lp3 = p;
  lp1.x = 0.0;
  lp1.y = 0.0;
  lp1.z = 0.0;
  magn_1 = dist3d (lp1, lp2);
  magn_2 = dist3d (lp1, lp3);
  if (magn_1 * magn_2 == 0)	/* emergency exit */
    return M_PI;
  cos_phi = scalar_prod (lp1, lp2, lp3) / (magn_1 * magn_2);
  if (cos_phi < -1)
    cos_phi = -1.0;
  if (cos_phi > 1)
    cos_phi = 1.0;
  res = acos (cos_phi);
  return res;
}


static double
torsion (p1, p2, p3, p4)
     p_3d p1, p2, p3, p4;
{
  p_3d lp1, lp2, lp3, lp4, d1, c1, c2;
  double res;
  p_3d c1xc2, c2xc1;
  double dist1, dist2, sign;

  /* copy everything into local variables */
  lp1 = p1;
  lp2 = p2;
  lp3 = p3;
  lp4 = p4;
  /* get the vector between the two central atoms */
  d1 = subtract_3d (p3, p2);
  /* shift the first atom parallel to be attached to p3 instead of p2 */
  lp1 = add_3d (p1, d1);
  /* now get the cross product vectors */
  c1 = cross_prod (lp3, lp2, lp1);
  c2 = cross_prod (lp3, lp2, lp4);
  res = angle_3d (p3, c1, c2);
  /*now check if it is clockwise or anticlockwise: */
  /*first, make the cross products of the two cross products c1 and c2 (both ways) */
  c1xc2 = cross_prod (lp3, c1, c2);
  c2xc1 = cross_prod (lp3, c2, c1);
  /*next, get the distances from these points to our refernce point lp2 */
  dist1 = dist3d (lp2, c1xc2);
  dist2 = dist3d (lp2, c2xc1);
  if (dist1 <= dist2)
    sign = 1.0;
  else
    sign = -1.0;
  return (sign * res);
}


static double
ctorsion (p1, p2, p3, p4)
     p_3d p1, p2, p3, p4;
{
  /* calculates "pseudo-torsion" defined by atoms 3 and 4, being both */
  /* attached to atom 2, with respect to axis of atoms 1 and 2 */
  p_3d lp1, lp2, lp3, lp4;
  /*d1 : p_3d; */
  p_3d c1, c2;
  double res;
  p_3d c1xc2, c2xc1;
  double dist1, dist2, sign;

  /* copy everything into local variables */
  lp1 = p1;
  lp2 = p2;
  lp3 = p3;
  lp4 = p4;
  /* get the cross product vectors */
  c1 = cross_prod (lp2, lp1, lp3);
  c2 = cross_prod (lp2, lp1, lp4);
  res = angle_3d (p2, c1, c2);
  /*now check if it is clockwise or anticlockwise: */
  /*first, make the cross products of the two cross products c1 and c2 (both ways) */
  c1xc2 = cross_prod (lp2, c1, c2);
  c2xc1 = cross_prod (lp2, c2, c1);
  /*next, get the distances from these points to our refernce point lp1 */
  dist1 = dist3d (lp1, c1xc2);
  dist2 = dist3d (lp1, c2xc1);
  if (dist1 <= dist2)
    sign = 1.0;
  else
    sign = -1.0;
  return (sign * res);
}


static boolean
is_cis (p1, p2, p3, p4)
     p_3d p1, p2, p3, p4;
{
  /* new in v0.3h, uses the dihedral angle */
  double phi;
  boolean res = false;

  phi = torsion (p1, p2, p3, p4);
  if (fabs (phi) < M_PI / 2)
    res = true;
  return res;
}


/*====================== end of geometry functions ========================== */

static void
show_usage ()
{
  if (progmode == pmMatchMol)
    {
      printf
	("matchmol version %s  N. Haider, University of Vienna, 2003-2007\n",
	 version);
      printf ("Usage: matchmol [options] <needle> <haystack>\n");
      printf
	(" where <needle> and <haystack> are the two molecules to compare\n");
      printf
	(" (supported formats: MDL *.mol or *.sdf, Alchemy *.mol, Sybyl *.mol2)\n");
      printf (" options can be:\n");
      printf ("    -v  verbose output\n");
      printf ("    -x  exact match\n");
      printf
	("    -s  strict comparison of atom and bond types (including ring check)\n");
      /* new in v0.2f, v0.3d */
      printf ("    -r  force SSR (set of small rings) ring search mode\n");
      printf
	("    -m  write matching molecule as MDL molfile to standard output\n");
      printf
	("        (default output: record number + \":T\" for hit  or \":F\" for miss\n");
      printf ("    -M  accept metal atoms as ring members\n");
      printf ("    -g  check geometry of double bonds (E/Z)\n");
      printf ("    -G  check geometry of chiral centers (R/S)\n");
      printf ("    -a  check charges strict\n");	/* 0.3x */
      printf ("    -i  check isotopes strict\n");	/* 0.3x */
      printf ("    -d  check radicals strict\n");	/* 0.3x */
      printf
	("    -f  fingerprint mode (1 haystack, multiple needles) with boolean output\n");
      printf
	("    -F  fingerprint mode (1 haystack, multiple needles) with decimal output\n");
      return;
    }
  printf ("checkmol version %s  N. Haider, University of Vienna, 2003-2007\n",
	  version);
  printf ("Usage: checkmol [options] <filename>\n");
  printf (" where options can be:\n");
  printf
    ("    -l  print a list of fingerprint codes + explanation and exit\n");
  printf ("    -v  verbose output\n");
  printf ("    -r  force SSR (set of small rings) ring search mode\n");
  printf ("    -M  accept metal atoms as ring members\n");
  printf ("  and one of the following:\n");
  printf
    ("    -e  english text (common name of functional group; default)\n");
  printf ("    -d  german text (common name of functional group)\n");
  printf ("    -c  code (acronym-like code for functional group)\n");
  printf
    ("    -b  binary (a bitstring representing absence or presence of each group)\n");
  printf
    ("    -s  the ASCII representation of the above bitstring, i.e. 0s and 1s)\n");
  printf
    ("    -x  print molecular statistics (number of various atom types, bond types,\n");
  printf ("        ring sizes, etc.\n");
  printf
    ("    -X  same as above, listing all records (even if 0) as comma-separated list\n");
  printf ("    -a  count charges in fingerprint\n");	/* 0.3x */
  printf
    ("    -m  write MDL molfile (with special encoding for aromatic atoms/bonds)\n");
  printf (" options can be combined like -vc\n");
  printf (" <filename> specifies any file in the formats supported\n");
  printf
    (" (MDL *.mol, Alchemy *.mol, Sybyl *.mol2), the filename \"-\" (without quotes)\n");
  printf (" specifies standard input\n");
  /* the "debug" option (-D) remains undocumented */
}


static void
list_molstat_codes ()
{
  printf ("n_atoms:     number of heavy atoms\n");
  printf ("n_bonds:     number of bonds between non-H atoms\n");
  printf ("n_rings:     number of rings\n");
  printf ("n_QA:        number of query atoms\n");
  printf ("n_QB:        number of query bonds\n");
  printf ("n_chg:       number of charges\n");
  printf ("n_C1:        number of sp-hybridized carbon atoms\n");
  printf ("n_C2:        number of sp2-hybridized carbon atoms\n");
  printf ("n_C:         total number of carbon atoms\n");
  printf
    ("n_CHB1p:     number of carbon atoms with at least 1 bond to a hetero atom\n");
  printf
    ("n_CHB2p:     number of carbon atoms with at least 2 bonds to a hetero atom\n");
  printf
    ("n_CHB3p:     number of carbon atoms with at least 3 bonds to a hetero atom\n");
  printf
    ("n_CHB4:      number of carbon atoms with 4 bonds to a hetero atom\n");
  printf ("n_O2:        number of sp2-hybridized oxygen atoms\n");
  printf ("n_O3:        number of sp3-hybridized oxygen atoms\n");
  printf ("n_N1:        number of sp-hybridized nitrogen atoms\n");
  printf ("n_N2:        number of sp2-hybridized nitrogen atoms\n");
  printf ("n_N3:        number of sp3-hybridized nitrogen atoms\n");
  printf ("n_S:         number of sulfur atoms\n");
  printf ("n_SeTe:      total number of selenium and tellurium atoms\n");
  printf ("n_F:         number of fluorine atoms\n");
  printf ("n_Cl:        number of chlorine atoms\n");
  printf ("n_Br:        number of bromine atoms\n");
  printf ("n_I:         number of iodine atoms\n");
  printf ("n_P:         number of phosphorus atoms\n");
  printf ("n_B:         number of boron atoms\n");
  printf ("n_Met:       total number of metal atoms\n");
  printf
    ("n_X:         total number of \"other\" atoms (not listed above) and halogens\n");
  printf ("n_b1:        number of single bonds\n");
  printf ("n_b2:        number of double bonds\n");
  printf ("n_b3:        number of triple bonds\n");
  printf ("n_bar:       number of aromatic bonds\n");
  printf ("n_C1O:       number of C-O single bonds\n");
  printf ("n_C2O:       number of C=O double bonds\n");
  printf ("n_CN:        number of C/N bonds (any type)\n");
  printf ("n_XY:        number of heteroatom/heteroatom bonds (any type)\n");
  printf ("n_r3:        number of 3-membered rings\n");
  printf ("n_r4:        number of 4-membered rings\n");
  printf ("n_r5:        number of 5-membered rings\n");
  printf ("n_r6:        number of 6-membered rings\n");
  printf ("n_r7:        number of 7-membered rings\n");
  printf ("n_r8:        number of 8-membered rings\n");
  printf ("n_r9:        number of 9-membered rings\n");
  printf ("n_r10:       number of 10-membered rings\n");
  printf ("n_r11:       number of 11-membered rings\n");
  printf ("n_r12:       number of 12-membered rings\n");
  printf ("n_r13p:      number of 13-membered or larger rings\n");
  printf ("n_rN:        number of rings containing nitrogen (any number)\n");
  printf ("n_rN1:       number of rings containing 1 nitrogen atom\n");
  printf ("n_rN2:       number of rings containing 2 nitrogen atoms\n");
  printf
    ("n_rN3p:      number of rings containing 3 or more nitrogen atoms\n");
  printf ("n_rO:        number of rings containing oxygen (any number)\n");
  printf ("n_rO1:       number of rings containing 1 oxygen atom\n");
  printf ("n_rO2p:      number of rings containing 2 or more oxygen atoms\n");
  printf ("n_rS:        number of rings containing sulfur (any number)\n");
  printf ("n_rX:        number of heterocycles (any type)\n");
  printf ("n_rar:       number of aromatic rings (any type)\n");
/* p2c: checkmol.pas, line 1207:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
  /*$IFDEF extended_molstat */
  printf ("n_rbz:       number of benzene rings\n");
  printf ("n_br2p:      number of bonds belonging to two or more rings\n");
  printf
    ("n_psg01:     number of atoms belonging to group 1 of the periodic system\n");
  printf
    ("n_psg02:     number of atoms belonging to group 2 of the periodic system\n");
  printf
    ("n_psg13:     number of atoms belonging to group 13 of the periodic system\n");
  printf
    ("n_psg14:     number of atoms belonging to group 14 of the periodic system\n");
  printf
    ("n_psg15:     number of atoms belonging to group 15 of the periodic system\n");
  printf
    ("n_psg16:     number of atoms belonging to group 16 of the periodic system\n");
  printf
    ("n_psg17:     number of atoms belonging to group 17 of the periodic system\n");
  printf
    ("n_psg18:     number of atoms belonging to group 18 of the periodic system\n");
  printf
    ("n_pstm:      number of atoms belonging to the transition metals\n");
  printf
    ("n_psla:      number of atoms belonging to the lanthanides or actinides\n");
  printf ("n_iso:      number of isotopes\n");
  printf ("n_rad:      number of radicals\n");
  /*$ENDIF */
}


/*static void parse_args()
{
  int p;
  char parstr[256];
  char tmpstr[256];
  int l;

  *tmpstr = '\0';
  opt_none = true;
  if (progmode == pmCheckMol) {
    for (p = 1; p < P_argc; p++) {
      strcpy(parstr, P_argv[p]);
      if (!strcmp(parstr, "-l")) {   /* new in v0.3l 
	list_molstat_codes();
	_Escape(0);
      }
      if (p < P_argc - 1) {
	if (strpos2(parstr, "-", 1) == 1 && p < P_argc - 1) {
	  strcpy(tmpstr, P_argv[p]);
	  left_trim(tmpstr);
	  l = 0;
	  if (strpos2(tmpstr, "v", 1) > 0)
	    l++;
	  if (strpos2(tmpstr, "D", 1) > 0)
	    l++;
	  if (strpos2(tmpstr, "r", 1) > 0)
	    l++;
	  if (strpos2(tmpstr, "M", 1) > 0)   /* new in v0.3 
	    l++;
	  if (strlen(tmpstr) > l + 2) {
	    show_usage();
	    _Escape(1);
	  }
	  opt_none = false;
	  if (strpos2(tmpstr, "M", 1) > 0)
	    opt_metalrings = true;
	  if (strpos2(tmpstr, "v", 1) > 0)
	    opt_verbose = true;
/* p2c: checkmol.pas, line 1261:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] 
	  /*$IFDEF debug
	  if (strpos2(tmpstr, "D", 1) > 0)
	    opt_debug = true;
	  /*$ENDIF
	  if (strpos2(tmpstr, "e", 1) > 0)
	    opt_text = true;
	  else {
	    if (strpos2(tmpstr, "d", 1) > 0)
	      opt_text_de = true;
	    else {
	      if (strpos2(tmpstr, "c", 1) > 0)
		opt_code = true;
	      else {
		if (strpos2(tmpstr, "b", 1) > 0)
		  opt_bin = true;
		else {
		  if (strpos2(tmpstr, "s", 1) > 0)
		    opt_bitstring = true;
		}
	      }
	    }
	    if (strpos2(tmpstr, "x", 1) > 0)
	      opt_molstat = true;
	    if (strpos2(tmpstr, "r", 1) > 0)
	      opt_rs = rs_ssr;
	    if (strpos2(tmpstr, "X", 1) > 0) {
	      opt_molstat = true;
	      opt_molstat_X = true;
	    }
	    if (strpos2(tmpstr, "m", 1) > 0) {
	      opt_text = false;
	      opt_text_de = false;
	      opt_bin = false;
	      opt_bitstring = false;
	      opt_code = false;
	      opt_molstat = false;
	      opt_xmdlout = true;
	    }
	  }
	  strcpy(molfilename, tmpstr);
	}
      } else {
	if (strpos2(parstr, "-", 1) == 1) {
	  if (strlen(parstr) > 1) {
	    show_usage();
	    _Escape(1);
	  }
	  opt_stdin = true;
	} else {
	  opt_stdin = false;
	  strcpy(molfilename, parstr);
	}
      }
    }
    if (opt_text == false && opt_text_de == false && opt_code == false &&
	opt_bin == false && opt_bitstring == false && opt_molstat == false &&
	opt_molstat_X == false && opt_xmdlout == false)
      opt_none = true;
  }
  if (progmode == pmMatchMol) {
    *ndl_molfilename = '\0';
    *molfilename = '\0';
    for (p = 1; p < P_argc; p++) {
      strcpy(parstr, P_argv[p]);
      if (p == 1) {
	if (strpos2(parstr, "-", 1) == 1) {
	  if (strpos2(parstr, "v", 1) > 1)
	    opt_verbose = true;
/* p2c: checkmol.pas, line 1329:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
	  /*$IFDEF debug
	     if (strpos2(parstr, "D", 1) > 1)
	     opt_debug = true;
	     /*$ENDIF
	     if (strpos2(parstr, "x", 1) > 1)
	     opt_exact = true;
	     if (strpos2(parstr, "s", 1) > 1)   /* new in v0.2f 
	     opt_strict = true;
	     if (strpos2(parstr, "m", 1) > 1)
	     opt_molout = true;
	     if (strpos2(parstr, "r", 1) > 1)
	     opt_rs = rs_ssr;
	     if (strpos2(parstr, "M", 1) > 0)   /* new in v0.3 
	     opt_metalrings = true;
	     if (strpos2(parstr, "g", 1) > 0)   /* new in v0.3d 
	     opt_geom = true;
	     if (strpos2(parstr, "G", 1) > 0)   /* new in v0.3f 
	     opt_chiral = true;
	     if (strpos2(parstr, "f", 1) > 0) {   /* new in v0.3m 
	     opt_fp = true;
	     fpformat = fpf_boolean;
	     }
	     if (strpos2(parstr, "F", 1) > 0) {   /* new in v0.3m 
	     opt_fp = true;
	     fpformat = fpf_decimal;
	     }
	     if (strpos2(parstr, "h", 1) > 1) {
	     show_usage();
	     _Escape(0);
	     }
	     } else
	     strcpy(ndl_molfilename, parstr);
	     }
	     if (p == P_argc - 2) {
	     if (strpos2(parstr, "-", 1) != 1)
	     strcpy(ndl_molfilename, parstr);
	     }
	     if (p == P_argc - 1) {
	     if (strcmp(parstr, "-"))
	     strcpy(molfilename, parstr);
	     else
	     opt_stdin = true;
	     }
	     }
	     if (opt_geom)   /* v0.3d 
	     ez_search = true;
	     if (opt_chiral)   /* v0.3f 
	     rs_search = true;
	     if (opt_chiral && opt_strict && (opt_exact || opt_fp))
	     /* new in v0.3j, v0.3m 
	     rs_strict = true;
	     if (opt_fp) {   /* v0.3m 
	     opt_molout = false;
	     opt_exact = false;
	     }
	     }  /* progmode = pmMatchMol 
	     ringsearch_mode = opt_rs;   /* v0.3i 
	     }
	   */

static void
parse_args (int argc, char *argv[])
{
  short p;
  char parstr[256];
  char tmpstr[256];
  short l;

  *tmpstr = '\0';
  opt_none = true;
  *molfilename = '\0';
  *ndl_molfilename = '\0';
  if (progmode == pmCheckMol)
    {
      for (p = 1; p <= argc - 1; p++)
	{
	  strcpy (parstr, argv[p]);
	  if (!strcmp (parstr, "-l"))
	    {			/* new in v0.3l */
	      list_molstat_codes ();
	      exit (0);
	    }
	  if (p < argc - 1)
	    {
	      if (strpos2 (parstr, "-", 1) == 1 && p < argc - 1)
		{
		  strcpy (tmpstr, argv[p]);
		  left_trim (tmpstr);
		  l = 0;
		  if (strpos2 (tmpstr, "v", 1) > 0)
		    l++;
		  if (strpos2 (tmpstr, "D", 1) > 0)
		    l++;
		  if (strpos2 (tmpstr, "r", 1) > 0)
		    l++;
		  /*if (strpos2 (tmpstr, "a", 1) > 0)   /* 0.3x 
		     l++; */
		  if (strpos2 (tmpstr, "M", 1) > 0)	/* new in v0.3 */
		    l++;
		  if (strlen (tmpstr) > l + 2)
		    {
		      show_usage ();
		      exit (1);
		    }
		  opt_none = false;
		  if (strpos2 (tmpstr, "M", 1) > 0)
		    opt_metalrings = true;
		  if (strpos2 (tmpstr, "v", 1) > 0)
		    opt_verbose = true;
		  /*{$IFDEF debug
		     if pos('D',tmpstr)>0 then opt_debug       := true;
		     {$ENDIF */
		  if (strpos2 (tmpstr, "e", 1) > 0)
		    opt_text = true;
		  else
		    {
		      if (strpos2 (tmpstr, "d", 1) > 0)
			opt_text_de = true;
		      else
			{
			  if (strpos2 (tmpstr, "c", 1) > 0)
			    opt_code = true;
			  else
			    {
			      if (strpos2 (tmpstr, "b", 1) > 0)
				opt_bin = true;
			      else
				{
				  if (strpos2 (tmpstr, "s", 1) > 0)
				    opt_bitstring = true;
				}
			    }
			}
		      if (strpos2 (tmpstr, "x", 1) > 0)
			opt_molstat = true;
		      if (strpos2 (tmpstr, "r", 1) > 0)
			opt_rs = rs_ssr;
		      /* if (strpos2 (tmpstr, "a", 1) > 0)
		         opt_chg = true; /* 0.3x  */
		      if (strpos2 (tmpstr, "X", 1) > 0)
			{
			  opt_molstat = true;
			  opt_molstat_X = true;
			}
		      if (strpos2 (tmpstr, "m", 1) > 0)
			{
			  opt_text = false;
			  opt_text_de = false;
			  opt_bin = false;
			  opt_bitstring = false;
			  opt_code = false;
			  opt_molstat = false;
			  opt_xmdlout = true;
			}
		    }
		  strcpy (molfilename, tmpstr);
		}
	    }
	  else
	    {
	      if (strpos2 (parstr, "-", 1) == 1)
		{
		  if (strlen (parstr) > 1)
		    {
		      show_usage ();
		      exit (1);
		    }
		  opt_stdin = true;
		}
	      else
		{
		  opt_stdin = false;
		  strcpy (molfilename, parstr);
		}
	    }
	}
      if (opt_text == false && opt_text_de == false && opt_code == false &&
	  opt_bin == false && opt_bitstring == false && opt_molstat == false
	  && opt_molstat_X == false && opt_xmdlout == false
	  && opt_chg == false)
	opt_none = true;	/* 0.3x */
    }
  if (progmode == pmMatchMol)
    {

      for (p = 1; p <= argc - 1; p++)
	{
	  strcpy (parstr, argv[p]);
	  if (p == 1)
	    {
	      if (strpos2 (parstr, "-", 1) == 1)
		{
		  if (strpos2 (parstr, "v", 1) > 1)
		    opt_verbose = true;
		  /*{$IFDEF debug
		     if pos('D',parstr)>1 then opt_debug       := true;
		     {$ENDIF */
		  if (strpos2 (parstr, "x", 1) > 1)
		    opt_exact = true;
		  if (strpos2 (parstr, "s", 1) > 1)	/* new in v0.2f */
		    opt_strict = true;
		  if (strpos2 (parstr, "m", 1) > 1)
		    opt_molout = true;
		  if (strpos2 (parstr, "r", 1) > 1)
		    opt_rs = rs_ssr;
		  if (strpos2 (parstr, "a", 1) > 0)
		    opt_chg = true;	/* 0.3x */
		  if (strpos2 (parstr, "i", 1) > 0)
		    opt_iso = true;	/* 0.3x */
		  if (strpos2 (parstr, "d", 1) > 0)
		    opt_rad = true;	/* 0.3x */
		  if (strpos2 (parstr, "M", 1) > 0)	/* new in v0.3 */
		    opt_metalrings = true;
		  if (strpos2 (parstr, "g", 1) > 0)	/* new in v0.3d */
		    opt_geom = true;
		  if (strpos2 (parstr, "G", 1) > 0)	/* new in v0.3f */
		    opt_chiral = true;
		  if (strpos2 (parstr, "f", 1) > 0)
		    {		/* new in v0.3m */
		      opt_fp = true;
		      fpformat = fpf_boolean;
		    }
		  if (strpos2 (parstr, "F", 1) > 0)
		    {		/* new in v0.3m */
		      opt_fp = true;
		      fpformat = fpf_decimal;
		    }
		  if (strpos2 (parstr, "h", 1) > 1)
		    {
		      show_usage ();
		      exit (0);
		    }
		}
	      else
		strcpy (ndl_molfilename, parstr);
	    }
	  if (p == argc - 2)
	    {
	      if (strpos2 (parstr, "-", 1) != 1)
		strcpy (ndl_molfilename, parstr);
	    }
	  if (p == argc - 1)
	    {
	      if (strcmp (parstr, "-"))
		strcpy (molfilename, parstr);
	      else
		opt_stdin = true;
	    }
	}
      if (opt_geom)		/* v0.3d */
	ez_search = true;
      if (opt_chiral)		/* v0.3f */
	rs_search = true;
      if (opt_chiral && opt_strict && (opt_exact || opt_fp))
	/* new in v0.3j, v0.3m  */
	rs_strict = true;
      if (opt_fp)
	{			/* v0.3m */
	  opt_molout = false;
	  opt_exact = false;
	}
    }				/* progmode = pmMatchMol */
  ringsearch_mode = opt_rs;	/* v0.3i */
}

/*============== input-related functions & procedures ===================== */

static char *
get_filetype (Result, f)
     char *Result;
     char *f;
{
  char rline[256];
  char auxstr[256];
  int i;
  boolean mdl1 = false;
  int ri;
  int sepcount = 0;
  char STR1[256], STR6[256], STR7[256];

  strcpy (auxstr, "unknown");
  i = li;
  ri = li - 1;
  while (ri < molbufindex && sepcount < 1)
    {
      ri++;
      strcpy (rline, molbuf[ri - 1]);
      if (strpos2 (rline, "$$$$", 1) > 0)
	sepcount++;
      if ((i == li) && (strcmp (strsub (STR1, rline, 7, 5), "ATOMS") == 0) &&
	  (strcmp (strsub (STR6, rline, 20, 5), "BONDS") == 0) &&
	  (strcmp (strsub (STR7, rline, 33, 7), "CHARGES") == 0))
	strcpy (auxstr, "alchemy");
      if ((i == li + 3) && (strcmp (strsub (STR1, rline, 35, 5), "V2000") ==
			    0))
	/* and (copy(rline,31,3)='999') */
	mdl1 = true;
      if ((i == li + 1) && (strcmp (strsub (STR1, rline, 3, 6), "-ISIS-") ==
			    0))
	mdl1 = true;
      if ((i == li + 1) && (strcmp (strsub (STR1, rline, 3, 8), "WLViewer") ==
			    0))
	mdl1 = true;
      if ((i == li + 1) && (strcmp (strsub (STR1, rline, 3, 8), "CheckMol") ==
			    0))
	mdl1 = true;
      if ((i == li + 1) && (strcmp (strsub (STR1, rline, 3, 8), "CATALYST") ==
			    0))
	{
	  mdl1 = true;
	  strcpy (auxstr, "mdl");
	}
      if (strpos2 (rline, "M  END", 1) == 1 || mdl1)
	strcpy (auxstr, "mdl");
      if (strpos2 (rline, "@<TRIPOS>MOLECULE", 1) > 0)
	strcpy (auxstr, "sybyl");
      i++;
    }
  /* new in v0.2j: try to identify non-conformant SD-files */
  if (!strcmp (auxstr, "unknown") && sepcount > 0)
    strcpy (auxstr, "mdl");
  return strcpy (Result, auxstr);
}


static void
zap_molecule ()
{
  /* try */
  if (atom != NULL)
    {
      free (atom);
      atom = NULL;		/* added in v0.3j */
    }
  if (bond != NULL)
    {
      free (bond);
      bond = NULL;		/* added in v0.3j */
    }
  if (ring != NULL)
    {
      free (ring);
      ring = NULL;		/* added in v0.3j */
    }
  if (ringprop != NULL)
    {
      free (ringprop);
      ringprop = NULL;		/* added in v0.3j */
    }
  /* except
     on e:Einvalidpointer do begin end;
     end; */
  n_atoms = 0;
  n_bonds = 0;
  n_rings = 0;
}


static void
zap_needle ()
{
  /* try */
  if (ndl_atom != NULL)
    {
      free (ndl_atom);
      ndl_atom = NULL;		/* added in v0.3j */
    }
  if (ndl_bond != NULL)
    {
      free (ndl_bond);
      ndl_bond = NULL;		/* added in v0.3j */
    }
  if (ndl_ring != NULL)
    {
      free (ndl_ring);
      ndl_ring = NULL;		/* added in v0.3j */
    }
  if (ndl_ringprop != NULL)
    {
      free (ndl_ringprop);	/* fixed in v0.3g */
      ndl_ringprop = NULL;	/* added in v0.3j */
    }
  /* except
     on e:Einvalidpointer do begin end;
     end; */
  ndl_n_atoms = 0;
  ndl_n_bonds = 0;
  ndl_n_rings = 0;
}


static void
zap_tmp ()
{
  /* try */
  if (tmp_atom != NULL)
    {
      free (tmp_atom);
      tmp_atom = NULL;		/* added in v0.3j */
    }
  if (tmp_bond != NULL)
    {
      free (tmp_bond);
      tmp_bond = NULL;		/* added in v0.3j */
    }
  if (tmp_ring != NULL)
    {
      free (tmp_ring);
      tmp_ring = NULL;		/* added in v0.3j */
    }
  if (tmp_ringprop != NULL)
    {
      free (tmp_ringprop);	/* fixed in v0.3g */
      tmp_ringprop = NULL;	/* added in v0.3j */
    }
  /* except
     on e:Einvalidpointer do begin end;
     end; */
  tmp_n_atoms = 0;
  tmp_n_bonds = 0;
  tmp_n_rings = 0;
}


static boolean
is_heavyatom (id)
     int id;
{
  str2 el;

  strcpy (el, atom[id - 1].element);
  
  if (!strcmp (el, "DU") || !strcmp (el, "LP"))
    return false;
  /*if (progmode == pmCheckMol && !strcmp (el, "H ")
     && atom[id - 1].nucleon_number < 2)
     return false;              /* 0.3x  */
  if (!strcmp (el, "H "))	/* 0.3 p */
    {
      if (progmode == pmMatchMol && !opt_iso)
	{
	  return false;
	}
      else
	{
	  if (atom[id - 1].nucleon_number < 2)
	    return false;
	}
    }
  return true;
}


static boolean
ndl_alkene_C (ba)
     int ba;
{
  /* new in v0.3f */
  boolean res = false;
  int i, ba2, FORLIM;

  if (ndl_n_atoms <= 0 || ndl_n_bonds <= 0)
    return false;
  FORLIM = ndl_n_bonds;
  for (i = 0; i < FORLIM; i++)
    {
      if (ndl_bond[i].a1 == ba || ndl_bond[i].a2 == ba)
	{
	  if (ndl_bond[i].a1 == ba)
	    ba2 = ndl_bond[i].a2;
	  else
	    ba2 = ndl_bond[i].a1;
	  if (!strcmp (ndl_atom[ba - 1].atype, "C2 ") &&
	      !strcmp (ndl_atom[ba2 - 1].atype, "C2 ")
	      && ndl_bond[i].btype == 'D' && ndl_bond[i].arom == false)
	    res = true;
	}
    }
  return res;
}


static boolean
is_metal (id)
     int id;
{
  boolean r = false;
  str2 el;

  strcpy (el, atom[id - 1].element);
  if (!strcmp (el, "LI") || !strcmp (el, "NA") || !strcmp (el, "K ") ||
      !strcmp (el, "RB") || !strcmp (el, "CS") || !strcmp (el, "BE") ||
      !strcmp (el, "MG") || !strcmp (el, "CA") || !strcmp (el, "SR") ||
      !strcmp (el, "BA") || !strcmp (el, "TI") || !strcmp (el, "ZR") ||
      !strcmp (el, "CR") || !strcmp (el, "MO") || !strcmp (el, "MN") ||
      !strcmp (el, "FE") || !strcmp (el, "CO") || !strcmp (el, "NI") ||
      !strcmp (el, "PD") || !strcmp (el, "PT") || !strcmp (el, "SN") ||
      !strcmp (el, "CU") || !strcmp (el, "AG") || !strcmp (el, "AU") ||
      !strcmp (el, "ZN") || !strcmp (el, "CD") || !strcmp (el, "HG") ||
      !strcmp (el, "AL") || !strcmp (el, "SN") || !strcmp (el, "PB") ||
      !strcmp (el, "SB") || !strcmp (el, "BI"))
/* p2c: checkmol.pas, line 1577: 
 * Note: Line breaker spent 0.0 seconds, 5000 tries on line 1686 [251] */
    /* etc. etc. */
    r = true;
  return r;
}


static int
get_nvalences (a_el)
     char *a_el;
{
  /* changed name and position in v0.3m */
  /* preliminary version; should be extended to element/atomtype */
  int res = 1;

  if (!strcmp (a_el, "H "))
    res = 1;
  /*if (!strcmp (a_el, "D "))   /* v0.3n 
     res = 1; */
  if (!strcmp (a_el, "C "))
    res = 4;
  if (!strcmp (a_el, "N "))
    res = 3;
  if (!strcmp (a_el, "O "))
    res = 2;
  if (!strcmp (a_el, "S "))
    res = 2;
  if (!strcmp (a_el, "SE"))
    res = 2;
  if (!strcmp (a_el, "TE"))
    res = 2;
  if (!strcmp (a_el, "P "))
    res = 3;
  if (!strcmp (a_el, "F "))
    res = 1;
  if (!strcmp (a_el, "CL"))
    res = 1;
  if (!strcmp (a_el, "BR"))
    res = 1;
  if (!strcmp (a_el, "I "))
    res = 1;
  if (!strcmp (a_el, "AT"))
    res = 1;
  if (!strcmp (a_el, "B "))
    res = 3;
  if (!strcmp (a_el, "LI"))
    res = 1;
  if (!strcmp (a_el, "NA"))
    res = 1;
  if (!strcmp (a_el, "K "))
    res = 1;
  if (!strcmp (a_el, "CA"))
    res = 2;
  if (!strcmp (a_el, "SR"))
    res = 2;
  if (!strcmp (a_el, "MG"))
    res = 2;
  if (!strcmp (a_el, "FE"))
    res = 3;
  if (!strcmp (a_el, "MN"))
    res = 2;
  if (!strcmp (a_el, "HG"))
    res = 2;
  if (!strcmp (a_el, "SI"))
    res = 4;
  if (!strcmp (a_el, "SN"))
    res = 4;
  if (!strcmp (a_el, "ZN"))
    res = 2;
  if (!strcmp (a_el, "CU"))
    res = 2;
  if (!strcmp (a_el, "A "))
    res = 4;
  if (!strcmp (a_el, "Q "))
    res = 4;
  return res;
}


static char *
convert_type (Result, oldtype)
     char *Result;
     char *oldtype;
{
  int i;
  str3 newtype;

  sprintf (newtype, "%.3s", oldtype);
  for (i = 0; i <= 2; i++)
    newtype[i] = toupper (newtype[i]);
  if (newtype[0] == '~')
    strcpy (newtype, "VAL");
  if (newtype[0] == '*')
    strcpy (newtype, "STR");
  return strcpy (Result, newtype);
}


static char *
convert_sybtype (Result, oldtype)
     char *Result;
     char *oldtype;
{
  str3 newtype;

  /*  NewType := Copy(OldType,1,3); */
  /*  For i := 1 To 3 Do NewType[i] := UpCase(NewType[i]); */
  /*  If NewType[1] = '~' Then NewType := 'VAL'; */
  /*  If NewType[1] = '*' Then NewType := 'STR'; */
  strcpy (newtype, "DU ");
  if (!strcmp (oldtype, "H    "))
    strcpy (newtype, "H  ");
  if (!strcmp (oldtype, "C.ar "))
    strcpy (newtype, "CAR");
  if (!strcmp (oldtype, "C.2  "))
    strcpy (newtype, "C2 ");
  if (!strcmp (oldtype, "C.3  "))
    strcpy (newtype, "C3 ");
  if (!strcmp (oldtype, "C.1  "))
    strcpy (newtype, "C1 ");
  if (!strcmp (oldtype, "O.2  "))
    strcpy (newtype, "O2 ");
  if (!strcmp (oldtype, "O.3  "))
    strcpy (newtype, "O3 ");
  if (!strcmp (oldtype, "O.co2"))
    strcpy (newtype, "O2 ");
  if (!strcmp (oldtype, "O.spc"))
    strcpy (newtype, "O3 ");
  if (!strcmp (oldtype, "O.t3p"))
    strcpy (newtype, "O3 ");
  if (!strcmp (oldtype, "N.1  "))
    strcpy (newtype, "N1 ");
  if (!strcmp (oldtype, "N.2  "))
    strcpy (newtype, "N2 ");
  if (!strcmp (oldtype, "N.3  "))
    strcpy (newtype, "N3 ");
  if (!strcmp (oldtype, "N.pl3"))
    strcpy (newtype, "NPL");
  if (!strcmp (oldtype, "N.4  "))
    strcpy (newtype, "N3+");
  if (!strcmp (oldtype, "N.am "))
    strcpy (newtype, "NAM");
  if (!strcmp (oldtype, "N.ar "))
    strcpy (newtype, "NAR");
  if (!strcmp (oldtype, "F    "))
    strcpy (newtype, "F  ");
  if (!strcmp (oldtype, "Cl   "))
    strcpy (newtype, "CL ");
  if (!strcmp (oldtype, "Br   "))
    strcpy (newtype, "BR ");
  if (!strcmp (oldtype, "I    "))
    strcpy (newtype, "I  ");
  if (!strcmp (oldtype, "Al   "))
    strcpy (newtype, "AL ");
  if (!strcmp (oldtype, "ANY  "))
    strcpy (newtype, "A  ");
  if (!strcmp (oldtype, "Ca   "))
    strcpy (newtype, "CA ");
  if (!strcmp (oldtype, "Du   "))
    strcpy (newtype, "DU ");
  if (!strcmp (oldtype, "Du.C "))
    strcpy (newtype, "DU ");
  if (!strcmp (oldtype, "H.spc"))
    strcpy (newtype, "H  ");
  if (!strcmp (oldtype, "H.t3p"))
    strcpy (newtype, "H  ");
  if (!strcmp (oldtype, "HAL  "))
    strcpy (newtype, "Cl ");
  if (!strcmp (oldtype, "HET  "))
    strcpy (newtype, "Q  ");
  if (!strcmp (oldtype, "HEV  "))
    strcpy (newtype, "DU ");
  if (!strcmp (oldtype, "K    "))
    strcpy (newtype, "K  ");
  if (!strcmp (oldtype, "Li   "))
    strcpy (newtype, "LI ");
  if (!strcmp (oldtype, "LP   "))
    strcpy (newtype, "LP ");
  if (!strcmp (oldtype, "Na   "))
    strcpy (newtype, "NA ");
  if (!strcmp (oldtype, "P.3  "))
    strcpy (newtype, "P3 ");
  if (!strcmp (oldtype, "S.2  "))
    strcpy (newtype, "S2 ");
  if (!strcmp (oldtype, "S.3  "))
    strcpy (newtype, "S3 ");
  if (!strcmp (oldtype, "S.o  "))
    strcpy (newtype, "SO ");
  if (!strcmp (oldtype, "S.o2 "))
    strcpy (newtype, "SO2");
  if (!strcmp (oldtype, "Si   "))
    strcpy (newtype, "SI ");
  if (!strcmp (oldtype, "P.4  "))
    strcpy (newtype, "P4 ");
  return strcpy (Result, newtype);
}


static char *
convert_MDLtype (Result, oldtype)
     char *Result, *oldtype;
{
  str3 newtype;

  /*  NewType := Copy(OldType,1,3); */
  /*  For i := 1 To 3 Do NewType[i] := UpCase(NewType[i]); */
  /*  If NewType[1] = '~' Then NewType := 'VAL'; */
  /*  If NewType[1] = '*' Then NewType := 'STR'; */
  strcpy (newtype, "DU ");
  if (!strcmp (oldtype, "H  "))
    strcpy (newtype, "H  ");
  if (!strcmp (oldtype, "C  "))
    strcpy (newtype, "C3 ");
  if (!strcmp (oldtype, "O  "))
    strcpy (newtype, "O2 ");
  if (!strcmp (oldtype, "N  "))
    strcpy (newtype, "N3 ");
  if (!strcmp (oldtype, "F  "))
    strcpy (newtype, "F  ");
  if (!strcmp (oldtype, "Cl "))
    strcpy (newtype, "CL ");
  if (!strcmp (oldtype, "Br "))
    strcpy (newtype, "BR ");
  if (!strcmp (oldtype, "I  "))
    strcpy (newtype, "I  ");
  if (!strcmp (oldtype, "Al "))
    strcpy (newtype, "AL ");
  if (!strcmp (oldtype, "ANY"))
    strcpy (newtype, "A  ");
  if (!strcmp (oldtype, "Ca "))
    strcpy (newtype, "CA ");
  if (!strcmp (oldtype, "Du "))
    strcpy (newtype, "DU ");
  if (!strcmp (oldtype, "K  "))
    strcpy (newtype, "K  ");
  if (!strcmp (oldtype, "Li "))
    strcpy (newtype, "LI ");
  if (!strcmp (oldtype, "LP "))
    strcpy (newtype, "LP ");
  if (!strcmp (oldtype, "Na "))
    strcpy (newtype, "NA ");
  if (!strcmp (oldtype, "P  "))
    strcpy (newtype, "P3 ");
  if (!strcmp (oldtype, "S  "))
    strcpy (newtype, "S3 ");
  if (!strcmp (oldtype, "Si "))
    strcpy (newtype, "SI ");
  if (!strcmp (oldtype, "P  "))
    strcpy (newtype, "P4 ");
  if (!strcmp (oldtype, "A  "))
    strcpy (newtype, "A  ");
  if (!strcmp (oldtype, "Q  "))
    strcpy (newtype, "Q  ");
  return strcpy (Result, newtype);
}


static char *
get_element (Result, oldtype)
     char *Result;
     char *oldtype;
{
  char elemstr[256];

  if (!strcmp (oldtype, "H   "))
    strcpy (elemstr, "H ");
  /* if (!strcmp (oldtype, "D   "))      /* v0.3n 
     strcpy (elemstr, "D "); */
  if (!strcmp (oldtype, "CAR "))
    strcpy (elemstr, "C ");
  if (!strcmp (oldtype, "C2  "))
    strcpy (elemstr, "C ");
  if (!strcmp (oldtype, "C3  "))
    strcpy (elemstr, "C ");
  if (!strcmp (oldtype, "C1  "))
    strcpy (elemstr, "C ");
  if (!strcmp (oldtype, "O2  "))
    strcpy (elemstr, "O ");
  if (!strcmp (oldtype, "O3  "))
    strcpy (elemstr, "O ");
  if (!strcmp (oldtype, "O2  "))
    strcpy (elemstr, "O ");
  if (!strcmp (oldtype, "O3  "))
    strcpy (elemstr, "O ");
  if (!strcmp (oldtype, "O3  "))
    strcpy (elemstr, "O ");
  if (!strcmp (oldtype, "N1  "))
    strcpy (elemstr, "N ");
  if (!strcmp (oldtype, "N2  "))
    strcpy (elemstr, "N ");
  if (!strcmp (oldtype, "N3  "))
    strcpy (elemstr, "N ");
  if (!strcmp (oldtype, "NPL "))
    strcpy (elemstr, "N ");
  if (!strcmp (oldtype, "N3+ "))
    strcpy (elemstr, "N ");
  if (!strcmp (oldtype, "NAM "))
    strcpy (elemstr, "N ");
  if (!strcmp (oldtype, "NAR "))
    strcpy (elemstr, "N ");
  if (!strcmp (oldtype, "F   "))
    strcpy (elemstr, "F ");
  if (!strcmp (oldtype, "CL  "))
    strcpy (elemstr, "CL");
  if (!strcmp (oldtype, "BR  "))
    strcpy (elemstr, "BR");
  if (!strcmp (oldtype, "I   "))
    strcpy (elemstr, "I ");
  if (!strcmp (oldtype, "AT  "))
    strcpy (elemstr, "AT");
  if (!strcmp (oldtype, "AL  "))
    strcpy (elemstr, "AL");
  if (!strcmp (oldtype, "DU  "))
    strcpy (elemstr, "DU");
  if (!strcmp (oldtype, "CA  "))
    strcpy (elemstr, "CA");
  if (!strcmp (oldtype, "DU  "))
    strcpy (elemstr, "DU");
  if (!strcmp (oldtype, "Cl  "))
    strcpy (elemstr, "CL");
  if (!strcmp (oldtype, "K   "))
    strcpy (elemstr, "K ");
  if (!strcmp (oldtype, "LI  "))
    strcpy (elemstr, "LI");
  if (!strcmp (oldtype, "LP  "))
    strcpy (elemstr, "LP");
  if (!strcmp (oldtype, "NA  "))
    strcpy (elemstr, "NA");
  if (!strcmp (oldtype, "P3  "))
    strcpy (elemstr, "P ");
  if (!strcmp (oldtype, "S2  "))
    strcpy (elemstr, "S ");
  if (!strcmp (oldtype, "S3  "))
    strcpy (elemstr, "S ");
  if (!strcmp (oldtype, "SO  "))
    strcpy (elemstr, "S ");
  if (!strcmp (oldtype, "SO2 "))
    strcpy (elemstr, "S ");
  if (!strcmp (oldtype, "SI  "))
    strcpy (elemstr, "SI");
  if (!strcmp (oldtype, "P4  "))
    strcpy (elemstr, "P ");
  if (!strcmp (oldtype, "A   "))
    strcpy (elemstr, "A ");
  if (!strcmp (oldtype, "Q   "))
    strcpy (elemstr, "Q ");
  return strcpy (Result, elemstr);
}


static char *
get_sybelement (Result, oldtype)
     char *Result;
     char *oldtype;
{
  int i;
  str2 elemstr;

  if (strpos2 (oldtype, ".", 1) < 2)
    sprintf (elemstr, "%.2s", oldtype);
  else
    {
      sprintf (elemstr, "%.*s", strpos2 (oldtype, ".", 1) - 1, oldtype);
      if (strlen (elemstr) < 2)
	strcat (elemstr, " ");
    }
  for (i = 0; i <= 1; i++)
    elemstr[i] = toupper (elemstr[i]);
  return strcpy (Result, elemstr);
}


static char *
get_MDLelement (Result, oldtype)
     char *Result;
     char *oldtype;
{
  int i;
  str2 elemstr;

  sprintf (elemstr, "%.2s", oldtype);
  for (i = 0; i <= 1; i++)
    elemstr[i] = toupper (elemstr[i]);
  if (elemstr[0] == '~')
    strcpy (elemstr, "??");
  if (elemstr[0] == '*')
    strcpy (elemstr, "??");
  return strcpy (Result, elemstr);
}


static void
read_molfile (mfilename)
     char *mfilename;
{
  /* reads ALCHEMY mol files */
  int n, code;
  char rline[256], tmpstr[256];
  char xstr[256], ystr[256], zstr[256], chgstr[256];
  float xval, yval, zval, chgval;
  char a1str[256], a2str[256], elemstr[256];
  int a1val, a2val, ri;
  char STR1[256];
  int FORLIM;
  atom_rec *WITH;
  bond_rec *WITH1;

  if (n_atoms > 0)
    zap_molecule ();
  ri = li;
  strcpy (rline, molbuf[ri - 1]);
  sprintf (tmpstr, "%.5s", rline);
  code = (sscanf (tmpstr, "%ld", &n_atoms) == 0);
  strsub (tmpstr, rline, 14, 5);
  code = (sscanf (tmpstr, "%ld", &n_bonds) == 0);
  strsub (molname, rline, 42, (int) (strlen (rline) - 42L));
  /* try */
  atom = safe_calloc (n_atoms, sizeof (atom_rec));
  bond = safe_calloc (n_bonds, sizeof (bond_rec));
  ring = safe_calloc (1, sizeof (ringlist));
  ringprop = safe_calloc (1, sizeof (ringprop_type));
  /* except
     on e:Eoutofmemory do
     begin
     writeln('Not enough memory');
     halt(4);
     end;
     end; */
  n_heavyatoms = 0;
  n_heavybonds = 0;
  n_Ctot = 0;			/* v0.3g */
  n_Otot = 0;			/* v0.3g */
  n_Ntot = 0;			/* v0.3g */
  FORLIM = n_atoms;
  for (n = 1; n <= FORLIM; n++)
    {
      ri++;
      strcpy (rline, molbuf[ri - 1]);
      strsub (atomtype, rline, 7, 4);
      sprintf (STR1, "%c", toupper (*atomtype));
      strcpy (atomtype, STR1);	/* fixed in v0.3f */
      get_element (elemstr, atomtype);
      if (!strcmp (elemstr, "C "))
	n_Ctot++;
      if (!strcmp (elemstr, "O "))
	n_Otot++;
      if (!strcmp (elemstr, "N "))
	n_Ntot++;
      convert_type (newatomtype, atomtype);
      strsub (xstr, rline, 14, 7);
      strsub (ystr, rline, 23, 7);
      strsub (zstr, rline, 32, 7);
      strsub (chgstr, rline, 43, 7);
      code = (sscanf (xstr, "%lg", &xval) == 0);
      code = (sscanf (ystr, "%lg", &yval) == 0);
      code = (sscanf (zstr, "%lg", &zval) == 0);
      code = (sscanf (chgstr, "%lg", &chgval) == 0);
      WITH = &atom[n - 1];
      strcpy (WITH->element, elemstr);
      strcpy (WITH->atype, newatomtype);
      WITH->x = xval;
      WITH->y = yval;
      WITH->z = zval;
      WITH->real_charge = chgval;
      if (is_heavyatom (n))
	{
	  n_heavyatoms++;
	  WITH->heavy = true;
	  if (is_metal (n))
	    WITH->metal = true;
	}
      WITH->nvalences = get_nvalences (WITH->element);	/* v0.3m   */
    }
  /*
     with atom^[n] do
     begin
     x := 0; y := 0; z := 0;  (* v0.3g
     formal_charge  := 0;
     real_charge    := 0;
     Hexp           := 0;
     Htot           := 0;
     neighbor_count := 0;
     ring_count     := 0;
     arom           := false;
     stereo_care    := false;
     heavy          := false;
     metal          := false;
     tag            := false;
     end;
   */
  FORLIM = n_bonds;
  for (n = 0; n < FORLIM; n++)
    {
      ri++;
      strcpy (rline, molbuf[ri - 1]);
      strsub (a1str, rline, 9, 3);
      strsub (a2str, rline, 15, 3);
      code = (sscanf (a1str, "%ld", &a1val) == 0);
      /* if code <> 0 then beep; */
      code = (sscanf (a2str, "%ld", &a2val) == 0);
      /* if code <> 0 then beep; */
      WITH1 = &bond[n];
      WITH1->a1 = a1val;
      WITH1->a2 = a2val;
      WITH1->btype = rline[19];
      WITH1->ring_count = 0;
      WITH1->arom = false;
      WITH1->topo = btopo_any;
      WITH1->stereo = bstereo_any;
      WITH1->mdl_stereo = 0;	/* v0.3n */
      if (atom[a1val - 1].heavy && atom[a2val - 1].heavy)
	n_heavybonds++;
    }
  memset (ring, 0, sizeof (ringlist));
  for (n = 0; n < max_rings; n++)
    {				/* new in v0.3 */
      ringprop[n].size = 0;
      ringprop[n].arom = false;
      ringprop[n].envelope = false;
    }
  li = ri + 1;
}


static void
read_mol2file (mfilename)
     char *mfilename;
{
  /* reads SYBYL mol2 files */
  int n, code;
  char sybatomtype[6];
  char tmpstr[256], rline[256];
  char xstr[256], ystr[256], zstr[256], chgstr[256];
  float xval, yval, zval, chgval;
  char a1str[256], a2str[256], elemstr[256];
  int a1val, a2val, ri, FORLIM;
  atom_rec *WITH;
  bond_rec *WITH1;

  if (n_atoms > 0)
    zap_molecule ();
  *rline = '\0';
  ri = li - 1;
  while ((ri < molbufindex) && (strpos2 (rline, "@<TRIPOS>MOLECULE", 1) == 0))
    {
      ri++;
      strcpy (rline, molbuf[ri - 1]);
    }
  if (ri < molbufindex)
    {
      ri++;
      strcpy (molname, molbuf[ri - 1]);
    }
  if (ri < molbufindex)
    {
      ri++;
      strcpy (rline, molbuf[ri - 1]);
    }
  sprintf (tmpstr, "%.5s", rline);
  sscanf (tmpstr, "%ld", &n_atoms);
/* p2c: checkmol.pas: Note: Eliminated unused assignment statement [338] */
  strsub (tmpstr, rline, 7, 5);
  sscanf (tmpstr, "%ld", &n_bonds);
  /* try */
/* p2c: checkmol.pas: Note: Eliminated unused assignment statement [338] */
  atom = safe_calloc (n_atoms, sizeof (atom_rec));
  bond = safe_calloc (n_bonds, sizeof (bond_rec));
  ring = safe_calloc (1, sizeof (ringlist));
  ringprop = safe_calloc (1, sizeof (ringprop_type));
  /* except
     on e:Eoutofmemory do
     begin
     writeln('Not enough memory');
     halt(4);
     end;
     end; */
  n_heavyatoms = 0;
  n_heavybonds = 0;
  n_Ctot = 0;			/* v0.3g */
  n_Otot = 0;			/* v0.3g */
  n_Ntot = 0;			/* v0.3g */
  while ((ri < molbufindex) && (strpos2 (rline, "@<TRIPOS>ATOM", 1) == 0))
    {
      ri++;
      strcpy (rline, molbuf[ri - 1]);
    }
  FORLIM = n_atoms;
  for (n = 1; n <= FORLIM; n++)
    {
      /*
         with atom^[n] do
         begin
         x := 0; y := 0; z := 0;  (* v0.3g
         formal_charge  := 0;
         real_charge    := 0;
         Hexp           := 0;
         Htot           := 0;
         neighbor_count := 0;
         ring_count     := 0;
         arom           := false;
         stereo_care    := false;
         heavy          := false;
         metal          := false;
         tag            := false;
         end;
       */
      if (ri < molbufindex)
	{
	  ri++;
	  strcpy (rline, molbuf[ri - 1]);
	}
      strsub (sybatomtype, rline, 48, 5);
      get_sybelement (elemstr, sybatomtype);
      if (!strcmp (elemstr, "C "))
	n_Ctot++;
      if (!strcmp (elemstr, "O "))
	n_Otot++;
      if (!strcmp (elemstr, "N "))
	n_Ntot++;
      convert_sybtype (newatomtype, sybatomtype);
      strsub (xstr, rline, 18, 9);
      strsub (ystr, rline, 28, 9);
      strsub (zstr, rline, 38, 9);
      strsub (chgstr, rline, 70, 9);
/* p2c: checkmol.pas: Note: Eliminated unused assignment statement [338] */
/* p2c: checkmol.pas: Note: Eliminated unused assignment statement [338] */
      sscanf (xstr, "%lg", &xval);
      sscanf (ystr, "%lg", &yval);
      sscanf (zstr, "%lg", &zval);
      sscanf (chgstr, "%lg", &chgval);
/* p2c: checkmol.pas: Note: Eliminated unused assignment statement [338] */
/* p2c: checkmol.pas: Note: Eliminated unused assignment statement [338] */
      WITH = &atom[n - 1];
      strcpy (WITH->element, elemstr);
      strcpy (WITH->atype, newatomtype);
      WITH->x = xval;
      WITH->y = yval;
      WITH->z = zval;
      WITH->real_charge = chgval;
      if (is_heavyatom (n))
	{
	  n_heavyatoms++;
	  WITH->heavy = true;
	  if (is_metal (n))
	    WITH->metal = true;
	}
      WITH->nvalences = get_nvalences (WITH->element);	/* v0.3m   */
    }
  while ((ri < molbufindex) && (strpos2 (rline, "@<TRIPOS>BOND", 1) == 0))
    {
      ri++;
      strcpy (rline, molbuf[ri - 1]);
    }
  FORLIM = n_bonds;
  for (n = 0; n < FORLIM; n++)
    {
      if (ri < molbufindex)
	{
	  ri++;
	  strcpy (rline, molbuf[ri - 1]);
	}
      strsub (a1str, rline, 9, 3);
      strsub (a2str, rline, 14, 3);
      code = (sscanf (a1str, "%ld", &a1val) == 0);
      if (code != 0)
	printf ("%s\007\n", rline);
      code = (sscanf (a2str, "%ld", &a2val) == 0);
      if (code != 0)
	printf ("%s\007\n", rline);
      WITH1 = &bond[n];
      WITH1->a1 = a1val;
      WITH1->a2 = a2val;
      if (rline[17] == '1')
	WITH1->btype = 'S';
      if (rline[17] == '2')
	WITH1->btype = 'D';
      if (rline[17] == '3')
	WITH1->btype = 'T';
      if (rline[17] == 'a')
	WITH1->btype = 'A';
      WITH1->ring_count = 0;
      WITH1->arom = false;
      WITH1->topo = btopo_any;
      WITH1->stereo = bstereo_any;
      WITH1->mdl_stereo = 0;	/* v0.3n */
      if (atom[a1val - 1].heavy && atom[a2val - 1].heavy)
	n_heavybonds++;
    }
  memset (ring, 0, sizeof (ringlist));
  for (n = 0; n < max_rings; n++)
    {				/* new in v0.3 */
      ringprop[n].size = 0;
      ringprop[n].arom = false;
      ringprop[n].envelope = false;
    }
  li = ri + 1;
}


static void
read_charges (chgstring_)
     char *chgstring_;
{
  char chgstring[256];
  int a_id, a_chg, n_chrg;

  /* typical example: a molecule with 2 cations + 1 anion */
  /* M  CHG  3   8   1  10   1  11  -1 */
  strcpy (chgstring, chgstring_);
  if (strpos2 (chgstring, "M  CHG", 1) <= 0)
    return;
  strdelete (chgstring, 1, 6);
  left_trim (chgstring);
  n_chrg = left_int (chgstring);
  /* this assignment must be kept also in non-debug mode! */
/* p2c: checkmol.pas, line 2077:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
  /*$IFDEF debug */
  /*if (n_chrg == 0)
     debugoutput ("strange... M  CHG present, but no charges found"); */
  /*$ENDIF */
  while (*chgstring != '\0')
    {
      a_id = left_int (chgstring);
      a_chg = left_int (chgstring);
      if (a_id != 0 && a_chg != 0)
	atom[a_id - 1].formal_charge = a_chg;
      //printf ("CHG %i %i\n", a_id, a_chg);
    }
}

static void
read_isotopes (char *isotopestring_)
{
  char isotopestring[256];
  int a_id, a_nucleon_number, n_isotopes;

  /* typical example: a molecule with 2 cations + 1 anion */
  /* M  CHG  3   8   1  10   1  11  -1 */
  strcpy (isotopestring, isotopestring_);
  if (strpos2 (isotopestring, "M  ISO", 1) <= 0)
    return;
  strdelete (isotopestring, 1, 6);
  left_trim (isotopestring);
  n_isotopes = left_int (isotopestring);
  /* this assignment must be kept also in non-debug mode! */
/* p2c: checkmol.pas, line 2077:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
  /*$IFDEF debug */
  /*if (n_chrg == 0)
     debugoutput ("strange... M  CHG present, but no charges found"); */
  /*$ENDIF */
  while (*isotopestring != '\0')
    {
      a_id = left_int (isotopestring);
      a_nucleon_number = left_int (isotopestring);
      if (a_id != 0 && a_nucleon_number != 0)
	{
	  atom[a_id - 1].nucleon_number = a_nucleon_number;
	  if (!strcmp (atom[a_id - 1].element, "H "))
	    {
	      atom[a_id - 1].heavy = true;
	      n_heavyatoms++;
	      strcpy (atom[a_id - 1].atype, "DU ");
	    }
	}
      //printf ("ISO %i %i\n", a_id, a_nucleon_number);
    }
}

static void
read_radicals (radstring_)
     char *radstring_;
{
  char radstring[256];
  int a_id, a_rad, n_rads;

  /* typical example: a molecule with 2 cations + 1 anion */
  /* M  CHG  3   8   1  10   1  11  -1 */
  strcpy (radstring, radstring_);
  if (strpos2 (radstring, "M  RAD", 1) <= 0)
    return;
  strdelete (radstring, 1, 6);
  left_trim (radstring);
  n_rads = left_int (radstring);
  /* this assignment must be kept also in non-debug mode! */
/* p2c: checkmol.pas, line 2077:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
  /*$IFDEF debug */
  /*if (n_chrg == 0)
     debugoutput ("strange... M  CHG present, but no charges found"); */
  /*$ENDIF */
  while (*radstring != '\0')
    {
      a_id = left_int (radstring);
      a_rad = left_int (radstring);
      if (a_id != 0 && a_rad != 0)
	atom[a_id - 1].radical_type = a_rad;
      //printf ("RAD %i %i\n", a_id, a_rad);
    }
}


static void
read_MDLmolfile (char *mfilename)
{
  /* reads MDL mol files */
  int n, v, tmp_n_atoms, tmp_n_bonds, code;	/* v0.3l */
  char rline[256], tmpstr[256];
  char xstr[256], ystr[256], zstr[256], chgstr[256];
  float xval, yval, zval, chgval;
  char a1str[256], a2str[256], elemstr[256];
  int a1val, a2val, ri, rc, bt, bs;
  int sepcount = 0;
  int i;			/* v0.3j */
  boolean clearcharges = true;	/* v0.3j */
  char STR1[256];
  int FORLIM;
  atom_rec *WITH;
  bond_rec *WITH1;

  /* v0.3j */
  if (n_atoms > 0)
    zap_molecule ();
  /*cm_mdlmolfile := false; */
  *rline = '\0';
  ri = li;
  strcpy (molname, molbuf[ri - 1]);	/* line 1 */
  if (ri < molbufindex)		/* line 2 */
    ri++;
  strcpy (rline, molbuf[ri - 1]);
  if (strpos2 (rline, "CheckMol", 1) == 3)
    {
      /*cm_mdlmolfile := true; */
      found_arominfo = true;
      tmfcode = 1;		/* v0.3m (begin) */
      code = 0;
      if ((strlen (rline) >= 39) && (strpos2 (rline, "TMF", 1) == 35))
	{			/* v0.3m; encoding of tweaklevel */
	  strsub (tmpstr, rline, 38, 2);
	  code = (sscanf (tmpstr, "%d", &tmfcode) == 0);
	}
      if (code != 0 || tmfcode != tweaklevel)
	tmfmismatch = true;
      else
	tmfmismatch = false;
      if ((strpos2 (rline, ":r0", 1) >= 40 && ringsearch_mode != rs_sar) |
	  (strpos2 (rline, ":r1", 1) >= 40 && ringsearch_mode != rs_ssr))
	tmfmismatch = true;
      if ((strpos2 (rline, ":m0", 1) >= 40 && opt_metalrings == true) |
	  (strpos2 (rline, ":m1", 1) >= 40 && opt_metalrings == false))
	tmfmismatch = true;
/* p2c: checkmol.pas, line 2128:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
      /*$IFDEF debug */
      //if (tmfmismatch)
      //printf ("\"tweaked\" molfile: version mismatch!\n");
      // else
      //printf ("\"tweaked\" molfile: version OK");
    }
  /*$ENDIF */
  /* v0.3m (end) */
  if (ri < molbufindex)		/* line 3 */
    ri++;
  strcpy (rline, molbuf[ri - 1]);
  strcpy (molcomment, rline);
  if (ri < molbufindex)		/* line 4 */
    ri++;
  strcpy (rline, molbuf[ri - 1]);
  sprintf (tmpstr, "%.3s", rline);
  sscanf (tmpstr, "%d", &n_atoms);
/* p2c: checkmol.pas: Note: Eliminated unused assignment statement [338] */
  strsub (tmpstr, rline, 4, 3);
  sscanf (tmpstr, "%d", &n_bonds);
/* p2c: checkmol.pas: Note: Eliminated unused assignment statement [338] */
  strsub (tmpstr, rline, 10, 3);
  /* if it is a CheckMol-tweaked molfile, this is the number of rings */
  n_cmrings = 0;
  code = (sscanf (tmpstr, "%d", &n_cmrings) == 0);
  if (code != 0)
    n_cmrings = 0;
  /* do some range checking for n_atoms, n_bonds; new in v0.3l */
  tmp_n_atoms = n_atoms;
  if (n_atoms > max_atoms)
    n_atoms = max_atoms;
  if (n_atoms < 0)
    n_atoms = 0;
  tmp_n_bonds = n_bonds;
  if (n_bonds > max_bonds)
    n_bonds = max_bonds;
  if (n_bonds < 0)
    n_bonds = 0;
  if (n_atoms == 0
#ifndef MAKE_SHARED_LIBRARY
      && opt_verbose
#endif
    )
    {				/* v0.3l */
      printf ("WARNING: Possible NoStruct read!\n");
      printf ("NoStructs are proprietary, obsolete and dangerous.\n");
    }
  /* try */
  atom = safe_calloc (n_atoms, sizeof (atom_rec));
  bond = safe_calloc (n_bonds, sizeof (bond_rec));
  /* this would be only one safe_calloc() in C;  v0.3l */
  ring = safe_calloc (1, sizeof (ringlist));
  ringprop = safe_calloc (1, sizeof (ringprop_type));
  /* except
     on e:Eoutofmemory do
     begin
     writeln('Not enough memory');
     (* close(molfile);
     halt(4);
     exit;
     end;
     end; */
  /* check for the chirality flag */
  if (strlen (rline) > 14 && rline[14] == '1')	/* new in v0.3f */
    chir_flag = true;
  n_heavyatoms = 0;
  n_heavybonds = 0;
  n_Ctot = 0;			/* v0.3g */
  n_Otot = 0;			/* v0.3g */
  n_Ntot = 0;			/* v0.3g */
  if (n_atoms > 0)
    {				/* v0.3l */
      for (n = 1; n <= tmp_n_atoms; n++)
	{
	  if (n <= max_atoms)
	    v = n;
	  else
	    v = max_atoms;
	  /* just for safety; v0.3l */
	  /*
	     with atom^[v] do
	     begin
	     x := 0; y := 0; z := 0;  (* v0.3g
	     formal_charge  := 0;
	     real_charge    := 0;
	     Hexp           := 0;
	     Htot           := 0;
	     neighbor_count := 0;
	     ring_count     := 0;
	     arom           := false;
	     stereo_care    := false;
	     metal          := false;
	     heavy          := false;
	     tag            := false;
	     end;
	   */
	  /* replaced by fillchar() after getmem() (see above); v0.3l */
	  ri++;
	  strcpy (rline, molbuf[ri - 1]);
	  strsub (atomtype, rline, 32, 3);
	  get_MDLelement (elemstr, atomtype);
	  if (!strcmp (elemstr, "C "))
	    n_Ctot++;
	  if (!strcmp (elemstr, "O "))
	    n_Otot++;
	  if (!strcmp (elemstr, "N "))
	    n_Ntot++;

	  convert_MDLtype (newatomtype, atomtype);
	  strsub (xstr, rline, 1, 10);	/* fixed in v0.3k (was: 2,9 etc.) */
	  strsub (ystr, rline, 11, 10);
	  strsub (zstr, rline, 21, 10);
	  /*chgstr := '0'; */
	  strsub (chgstr, rline, 37, 3);	/* new in v0.3j */
	  sscanf (chgstr, "%f", &chgval);
/* p2c: checkmol.pas: Note: Eliminated unused assignment statement [338] */
	  if (chgval != 0)
	    {
	      if (chgval >= 1 && chgval <= 7)
		chgval = 4.0 - chgval;
	      else
		{
		  chgval = 0.0;
/* p2c: checkmol.pas: Note: Eliminated unused assignment statement [338] */
		}
	    }			/* end (v0.3j) */
	  sscanf (xstr, "%f", &xval);
	  sscanf (ystr, "%f", &yval);
	  sscanf (zstr, "%f", &zval);
	  /* v0.3k: removed superfluous val(chgstr,chgval,code) */
/* p2c: checkmol.pas: Note: Eliminated unused assignment statement [338] */
/* p2c: checkmol.pas: Note: Eliminated unused assignment statement [338] */
	  WITH = &atom[v - 1];
	  strcpy (WITH->element, elemstr);
	  if (!strcmp (elemstr, "A ") || !strcmp (elemstr, "Q ") ||
	      !strcmp (elemstr, "X "))
	    /* 'X ' added in v0.3n */
	    found_querymol = true;

	  strcpy (WITH->atype, newatomtype);


	  if (!strcmp (elemstr, "D "))
	    {
	      strcpy (WITH->element, "H ");
	      WITH->nucleon_number = 2;
	    }			/* 0.3x */
	  if (!strcmp (elemstr, "T "))
	    {
	      strcpy (WITH->element, "H ");
	      WITH->nucleon_number = 3;
	    }			/* 0.3x */



	  WITH->x = xval;
	  WITH->y = yval;
	  WITH->z = zval;
	  WITH->formal_charge = (long) floor (chgval + 0.5);
	  WITH->real_charge = 0.0;	/* v0.3j */
	  /* read aromaticity flag from CheckMol-tweaked MDL molfile */
	  if (strlen (rline) > 37 && rline[37] == '0')
	    {
	      WITH->arom = true;
	      found_arominfo = true;
	    }
	  /* new in v0.3d: read stereo care flag */
	  if (strlen (rline) > 47 && rline[47] == '1')
	    WITH->stereo_care = true;
	  if (is_heavyatom (n))
	    {
	      n_heavyatoms++;
	      WITH->heavy = true;
	      if (is_metal (n))
		WITH->metal = true;
	    }
	  WITH->nvalences = get_nvalences (WITH->element);
	  /* v0.3m                 */
	}
    }				/* if (n_atoms > 0)... */
  if (n_bonds > 0)
    {				/* v0.3l */
      for (n = 1; n <= tmp_n_bonds; n++)
	{
	  if (n <= max_bonds)
	    v = n;
	  else
	    v = max_bonds;
	  /* just for safety; v0.3l */
	  ri++;
	  strcpy (rline, molbuf[ri - 1]);
	  sprintf (a1str, "%.3s", rline);
	  strsub (a2str, rline, 4, 3);
	  code = (sscanf (a1str, "%d", &a1val) == 0);
	  if (code != 0)	/* v0.3l */
	    a1val = 1;
	  code = (sscanf (a2str, "%d", &a2val) == 0);
	  if (code != 0)	/* v0.3l */
	    a2val = 1;
	  WITH1 = &bond[v - 1];
	  WITH1->a1 = a1val;
	  WITH1->a2 = a2val;
	  if (rline[8] == '1')	/* single */
	    WITH1->btype = 'S';
	  if (rline[8] == '2')	/* double */
	    WITH1->btype = 'D';
	  if (rline[8] == '3')	/* triple */
	    WITH1->btype = 'T';
	  if (rline[8] == '4')	/* aromatic */
	    WITH1->btype = 'A';
	  if (rline[8] == '5')	/* single or double */
	    WITH1->btype = 'l';
	  if (rline[8] == '6')	/* single or aromatic */
	    WITH1->btype = 's';
	  if (rline[8] == '7')	/* double or aromatic */
	    WITH1->btype = 'd';
	  if (rline[8] == '8')	/* any */
	    WITH1->btype = 'a';
	  sprintf (STR1, "%c", WITH1->btype);
	  if (strpos2 ("lsda", STR1, 1) > 0)
	    found_querymol = true;
	  WITH1->arom = false;
	  WITH1->q_arom = false;	/* 0.3p */
	  /* read aromaticity flag from CheckMol-tweaked MDL molfile */
	  if (WITH1->btype == 'A' || rline[7] == '0')
	    {
	      WITH1->arom = true;
	      if (rline[7] == '0')
		found_arominfo = true;
	    }
	  strsub (tmpstr, rline, 13, 3);
	  /* new in v0.3d: read ring_count from tweaked molfile */
	  code = (sscanf (tmpstr, "%d", &rc) == 0);
	  if (code != 0 || rc < 0 || progmode == pmCheckMol || tmfmismatch)
	    WITH1->ring_count = 0;
	  else
	    WITH1->ring_count = rc;
	  /* v0.3n: added tmfmismatch check */
	  strsub (tmpstr, rline, 16, 3);	/* new in v0.3d: read bond topology; */
	  code = (sscanf (tmpstr, "%d", &bt) == 0);
	  /* extended features are encoded by leading zero */
	  if (code != 0 || (unsigned long) bt > 2)
	    WITH1->topo = btopo_any;
	  else
	    {
	      if (tmpstr[1] == '0')
		WITH1->topo = bt + 3;
	      else
		WITH1->topo = bt;
	    }
	  /* v0.3n changed >5 into >2 */
	  /* new in v0.3d: add stereo property from MDL "stereo care" flag in atom block */
	  WITH1->stereo = bstereo_any;
	  if (WITH1->btype == 'D')
	    {
	      if (atom[WITH1->a1 - 1].stereo_care
		  && atom[WITH1->a2 - 1].stereo_care)
		{		/* this is the MDL-conformant encoding, */
		  WITH1->stereo = bstereo_xyz;	/* for an alternative see below */
		  ez_flag = true;	/* v0.3f */
		}
	      else
		{		/* this extended feature is encoded by a leading zero */
		  strsub (tmpstr, rline, 10, 3);
		  /* new in v0.3d: read bond stereo specification; */
		  code = (sscanf (tmpstr, "%d", &bs) == 0);
		  WITH1->mdl_stereo = bs;	/* v0.3n */
		  if (code != 0 || bs <= 0 || bs > 2)
		    WITH1->stereo = bstereo_any;
		  else
		    WITH1->stereo = bstereo_xyz;
		  if (tmpstr[1] == '0')
		    WITH1->stereo = bstereo_xyz;
		}
	    }
	  /*if stereo <> bstereo_any then ez_search := true; */
	  if (WITH1->stereo != bstereo_any)	/* changed in v0.3f */
	    ez_flag = true;
	  if (WITH1->btype == 'S' && strlen (rline) > 11 && rline[11] == '1')
	    WITH1->stereo = bstereo_up;
	  if (WITH1->btype == 'S' && strlen (rline) > 11 && rline[11] == '6')
	    WITH1->stereo = bstereo_down;
	  if (WITH1->btype == 'S' && strlen (rline) > 11 && rline[11] == '4')
	    WITH1->stereo = bstereo_either;	/* 0.3x */
	  if (WITH1->btype == 'D' && strlen (rline) > 11 && rline[11] == '3')
	    WITH1->stereo = bstereo_double_either;	/* 0.3x */
	  strsub (tmpstr, rline, 10, 3);
	  /* new in v0.3n: save original bond stereo specification; */
	  sscanf (tmpstr, "%d", &bs);
	  /* v0.3n */
/* p2c: checkmol.pas: Note: Eliminated unused assignment statement [338] */
	  WITH1->mdl_stereo = bs;	/* v0.3n */
	  if (atom[a1val - 1].heavy && atom[a2val - 1].heavy)
	    n_heavybonds++;
	}
    }				/* if (n_bonds > 0)... */
  while (ri < molbufindex && sepcount < 1)
    {
      ri++;
      strcpy (rline, molbuf[ri - 1]);
      if (strpos2 (rline, "M  CHG", 1) > 0)
	{			/* new in v0.3j */
	  if (clearcharges)
	    {			/* "M  CHG" supersedes all "old-style" charge values */

	      for (i = 0; i < n_atoms; i++)
		atom[i].formal_charge = 0;
	    }
	  read_charges (rline);
	  clearcharges = false;
	  /* subsequent "M  CHG" lines must not clear previous values */
	}

      if (strpos2 (rline, "M  ISO", 1) > 0)
	read_isotopes (rline);	/* 0.3x */

      if (strpos2 (rline, "M  RAD", 1) > 0)
	read_radicals (rline);	/* 0.3x */

      if (strpos2 (rline, "$$$$", 1) > 0)
	{
	  sepcount++;
	  if (molbufindex > ri + 2)	/* we assume this is an SDF file */
	    mol_in_queue = true;
	}
    }
  memset (ring, 0, sizeof (ringlist));
  for (n = 0; n < max_rings; n++)
    {				/* new in v0.3 */
      ringprop[n].size = 0;
      ringprop[n].arom = false;
      ringprop[n].envelope = false;
    }
  li = ri + 1;
}



static void
write_MDLmolfile ()
{
  int i;
  char tmpstr[256];
  char wline[256];
  int a_chg;
  int a_iso;
  int a_rad;
  char tmflabel[256];		/* v0.3m */
  char STR1[256], STR7[256];
  int FORLIM;

  sprintf (tmflabel, "%d", (int) tweaklevel);	/* v0.3m */
  while (strlen (tmflabel) < 2)	/* v0.3m */
    sprintf (tmflabel, "0%s", strcpy (STR1, tmflabel));
  sprintf (tmflabel, "TMF%s", strcpy (STR1, tmflabel));	/* v0.3m */
  if (strlen (molname) > 80)
    sprintf (molname, "%.80s", strcpy (STR1, molname));
  puts (molname);
  printf ("  CheckMol                        %s", tmflabel);	/* v0.3m */
  if (ringsearch_mode == rs_sar)	/* v0.3m */
    printf (":r0");
  if (ringsearch_mode == rs_ssr)	/* v0.3m */
    printf (":r1");
  if (opt_metalrings)
    printf (":m1");
  else
    printf (":m0");
  /* v0.3m */
  printf ("\n%s\n", molcomment);
  *wline = '\0';
  *tmpstr = '\0';
  sprintf (tmpstr, "%d", n_atoms);
  lblank (3L, tmpstr);
  strcat (wline, tmpstr);
  *tmpstr = '\0';		/* first 3 digits: number of atoms */
  sprintf (tmpstr, "%d", n_bonds);
  lblank (3L, tmpstr);
  strcat (wline, tmpstr);
  *tmpstr = '\0';		/* next 3 digits: number of bonds */
  strcpy (tmpstr, "  0");
  strcat (wline, tmpstr);
  *tmpstr = '\0';		/* next 3 digits: number of atom lists (not used by us) */
/* p2c: checkmol.pas, line 2388:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
#ifdef REDUCED_SAR
  sprintf (tmpstr, "%d", n_countablerings);
  /* v0.3n; changed n_rings into n_countablerings */
#else
  sprintf (tmpstr, "%d", n_rings);
#endif
  lblank (3L, tmpstr);
  strcat (wline, tmpstr);
  *tmpstr = '\0';
  /* officially "obsolete", we use it for the number of rings */
  strcat (wline, "  ");		/* v0.3n: obey chiral flag */
  if (chir_flag)
    strcat (wline, "1");
  else
    strcat (wline, "0");
  /* v0.3n */
  strcat (wline, "               999 V2000");
  /* v0.3n (adjust string length) */
  puts (wline);
  FORLIM = n_atoms;
  for (i = 0; i < FORLIM; i++)
    {
      *wline = '\0';
      sprintf (tmpstr, "%1.4f", atom[i].x);
      lblank (10L, tmpstr);
      strcat (wline, tmpstr);
      sprintf (tmpstr, "%1.4f", atom[i].y);
      lblank (10L, tmpstr);
      strcat (wline, tmpstr);
      sprintf (tmpstr, "%1.4f", atom[i].z);
      lblank (10L, tmpstr);
      strcat (wline, tmpstr);
      strcpy (tmpstr, atom[i].element);
      /* tmpstr := lowercase(tmpstr); REPLACE!!! */
      //tmpstr[0] = toupper (tmpstr[0]);
      all_lowercase (tmpstr);
      tmpstr[0] = toupper (tmpstr[0]);
      /*wline := wline + ' '+atom^[i].element+' '; */
      sprintf (wline + strlen (wline), " %s ", tmpstr);
      strcat (wline, " 0");	/* mass difference (isotopes) */
      /* now we code aromaticity into the old-style charge column (charges are now in the M  CHG line) */
      if (atom[i].arom)
	strcpy (tmpstr, " 00");
      else
	strcpy (tmpstr, "  0");
      strcat (wline, tmpstr);
      strcat (wline, "  0  0  0  0  0  0  0  0  0  0");
      puts (wline);
    }
  FORLIM = n_bonds;
  for (i = 0; i < FORLIM; i++)
    {
      *wline = '\0';
      sprintf (tmpstr, "%d", bond[i].a1);
      lblank (3L, tmpstr);
      strcat (wline, tmpstr);
      sprintf (tmpstr, "%d", bond[i].a2);
      lblank (3L, tmpstr);
      strcat (wline, tmpstr);
      if (bond[i].btype == 'S')
	strcpy (tmpstr, "  1");
      if (bond[i].btype == 'D')
	strcpy (tmpstr, "  2");
      if (bond[i].btype == 'T')
	strcpy (tmpstr, "  3");
      if (bond[i].btype == 'A')
	strcpy (tmpstr, "  4");
      if (bond[i].btype == 'l')
	strcpy (tmpstr, "  5");
      if (bond[i].btype == 's')
	strcpy (tmpstr, "  6");
      if (bond[i].btype == 'd')
	strcpy (tmpstr, "  7");
      if (bond[i].btype == 'a')
	strcpy (tmpstr, "  8");
      /* now encode our own aromaticity information */
      if (bond[i].arom)
	tmpstr[1] = '0';
      strcat (wline, tmpstr);	/* next, encode bond stereo property (v0.3f) */
      /*if (bond^[i].stereo = bstereo_up) then wline := wline + '  1' else */
      /*  if (bond^[i].stereo = bstereo_down) then wline := wline + '  6' else */
      /*    wline := wline + '  0'; */
      /* restore original value from MDL molfile (v0.3n) */
      /* wline := wline + '  ' + inttostr(bond^[i].mdl_stereo);    REPLACE!!! */
      *tmpstr = '\0';
      sprintf (tmpstr, "%i", bond[i].mdl_stereo);
      strcat (wline, "  ");
      strcat (wline, tmpstr);
      *tmpstr = '\0';
      /* now encode the ring_count of this bond (using a field which officially is "not used") */
      /* tmpstr := inttostr(bond^[i].ring_count); REPLACE!!! */
      sprintf (tmpstr, "%i", bond[i].ring_count);
      while (strlen (tmpstr) < 3)
	sprintf (tmpstr, " %s", strcpy (STR1, tmpstr));
      sprintf (wline + strlen (wline), "%s  0  0", tmpstr);
      puts (wline);
    }
  FORLIM = n_atoms;
  for (i = 1; i <= FORLIM; i++)
    {
      a_chg = atom[i - 1].formal_charge;
      if (a_chg != 0)
	{
	  strcpy (wline, "M  CHG  1 ");
	  sprintf (tmpstr, "%d", i);
	  lblank (3L, tmpstr);
	  sprintf (wline + strlen (wline), "%s ", tmpstr);
	  sprintf (tmpstr, "%d", a_chg);
	  lblank (3L, tmpstr);
	  strcat (wline, tmpstr);
	  puts (wline);
	}
    }
  for (i = 1; i <= FORLIM; i++)	/* 0.3x */
    {
      a_iso = atom[i - 1].nucleon_number;
      if (a_iso != 0)
	{
	  strcpy (wline, "M  ISO  1 ");
	  sprintf (tmpstr, "%d", i);
	  lblank (3L, tmpstr);
	  sprintf (wline + strlen (wline), "%s ", tmpstr);
	  sprintf (tmpstr, "%d", a_iso);
	  lblank (3L, tmpstr);
	  strcat (wline, tmpstr);
	  puts (wline);
	}
    }
  for (i = 1; i <= FORLIM; i++)	/* 0.3x */
    {
      a_rad = atom[i - 1].radical_type;
      if (a_rad != 0)
	{
	  strcpy (wline, "M  RAD  1 ");
	  sprintf (tmpstr, "%d", i);
	  lblank (3L, tmpstr);
	  sprintf (wline + strlen (wline), "%s ", tmpstr);
	  sprintf (tmpstr, "%d", a_rad);
	  lblank (3L, tmpstr);
	  strcat (wline, tmpstr);
	  puts (wline);
	}
    }
  printf ("M  END\n");
}


/*============= chemical processing functions && procedures ============ */

static boolean
is_electroneg (a_el)
     char *a_el;
{
  /* new in v0.3j */
  boolean res = false;

  if (!strcmp (a_el, "N "))
    res = true;
  if (!strcmp (a_el, "P "))
    res = true;
  if (!strcmp (a_el, "O "))
    res = true;
  if (!strcmp (a_el, "S "))
    res = true;
  if (!strcmp (a_el, "SE"))
    res = true;
  if (!strcmp (a_el, "TE"))
    res = true;
  if (!strcmp (a_el, "F "))
    res = true;
  if (!strcmp (a_el, "CL"))
    res = true;
  if (!strcmp (a_el, "BR"))
    res = true;
  if (!strcmp (a_el, "I "))
    res = true;
  if (!strcmp (a_el, "AT"))
    res = true;
  return res;
}


static void
count_neighbors ()
{
  /* counts heavy-atom neighbors and explicit hydrogens */
  int i, FORLIM;

  if (n_atoms < 1 || n_bonds < 1)
    return;
  FORLIM = n_bonds;
  for (i = 0; i < FORLIM; i++)
    {
      if (atom[bond[i].a1 - 1].heavy)
	atom[bond[i].a2 - 1].neighbor_count++;
      if (atom[bond[i].a2 - 1].heavy)
	atom[bond[i].a1 - 1].neighbor_count++;
      if (!strcmp (atom[bond[i].a1 - 1].element, "H "))
	atom[bond[i].a2 - 1].Hexp++;
      if (!strcmp (atom[bond[i].a2 - 1].element, "H "))
	atom[bond[i].a1 - 1].Hexp++;
      /* plausibility check (new in v02.i) */
      if (atom[bond[i].a1 - 1].neighbor_count > max_neighbors ||
	  atom[bond[i].a2 - 1].neighbor_count > max_neighbors)
	{
	  mol_OK = false;
	  /*writeln('invalid molecule!'); */
	}
    }
}


static void
get_neighbors (Result, id)
     int *Result;
     int id;
{
  int i = 0;
  //neighbor_rec nb_tmp;
  int nb_count = 0;
  //int FORLIM = n_bonds;

  //memset (Result, 0, sizeof (neighbor_rec));

  for (i; i < n_bonds; i++)
    {
      if (bond[i].a1 == id && atom[bond[i].a2 - 1].heavy
	  && nb_count < max_neighbors)
	{
	  Result[nb_count++] = bond[i].a2;
	}
      if (bond[i].a2 == id && atom[bond[i].a1 - 1].heavy
	  && nb_count < max_neighbors)
	{
	  Result[nb_count++] = bond[i].a1;
	}
    }
  //return memcpy (Result, nb_tmp, sizeof (neighbor_rec));
}


static void
get_ndl_neighbors (int *Result, int id)
{
  int i;
  //neighbor_rec nb_tmp;
  int nb_count = 0;
  /* v0.3i: use max_neighbors instead of a fixed value of 8 */
  //int FORLIM = ndl_n_bonds;

  //memset (nb_tmp, 0, sizeof (neighbor_rec));

  for (i = 0; i < ndl_n_bonds; i++)
    {
      if (ndl_bond[i].a1 == id && nb_count < max_neighbors &&
	  ndl_atom[ndl_bond[i].a2 - 1].heavy)
	{
	  nb_count++;
	  Result[nb_count - 1] = ndl_bond[i].a2;
	}
      if (ndl_bond[i].a2 == id && nb_count < max_neighbors &&
	  ndl_atom[ndl_bond[i].a1 - 1].heavy)
	{
	  nb_count++;
	  Result[nb_count - 1] = ndl_bond[i].a1;
	}
    }
  //return memcpy (Result, nb_tmp, sizeof (neighbor_rec));
}


static int *
get_nextneighbors (int *Result, int id, int prev_id)
{
  int i;
  neighbor_rec nb_tmp;
  int nb_count = 0;
  int FORLIM;

  /* gets all neighbors except prev_id (usually the atom where we came from */
  memset (nb_tmp, 0, sizeof (neighbor_rec));
  FORLIM = n_bonds;
  for (i = 0; i < FORLIM; i++)
    {
      if (bond[i].a1 == id && bond[i].a2 != prev_id &&
	  nb_count < max_neighbors && atom[bond[i].a2 - 1].heavy)
	{
	  nb_count++;
	  nb_tmp[nb_count - 1] = bond[i].a2;
	}
      if (bond[i].a2 == id && bond[i].a1 != prev_id &&
	  nb_count < max_neighbors && atom[bond[i].a1 - 1].heavy)
	{
	  nb_count++;
	  nb_tmp[nb_count - 1] = bond[i].a1;
	}
    }
  return memcpy (Result, nb_tmp, sizeof (neighbor_rec));
}


static int *
get_ndl_nextneighbors (int *Result, int id, int prev_id)
{
  int i;
  neighbor_rec nb_tmp;
  int nb_count = 0;
  int FORLIM;

  /* gets all neighbors except prev_id (usually the atom where we came from */
  memset (nb_tmp, 0, sizeof (neighbor_rec));
  FORLIM = ndl_n_bonds;
  for (i = 0; i < FORLIM; i++)
    {
      if (ndl_bond[i].a1 == id && ndl_bond[i].a2 != prev_id &&
	  nb_count < max_neighbors && ndl_atom[ndl_bond[i].a2 - 1].heavy)
	{
	  nb_count++;
	  nb_tmp[nb_count - 1] = ndl_bond[i].a2;
	}
      if (ndl_bond[i].a2 == id && ndl_bond[i].a1 != prev_id &&
	  nb_count < max_neighbors && ndl_atom[ndl_bond[i].a1 - 1].heavy)
	{
	  nb_count++;
	  nb_tmp[nb_count - 1] = ndl_bond[i].a1;
	}
    }
  return memcpy (Result, nb_tmp, sizeof (neighbor_rec));
}


  /*static int path_pos (id, a_path) int id;
     int *a_path;
     {
     /* new version in v0.3l 
     int i;
     int pp = 0;

     for (i = 1; i <= max_ringsize; i++)
     {
     if (a_path[i - 1] == id)
     {
     pp = i;
     /* p2c: checkmol.pas, line 2620:
     * Warning: Expected a '(', found a semicolon [227] 
     /* p2c: checkmol.pas, line 2620:
     * Warning: Expected an expression, found a semicolon [227] 
     fflush (0);
     P_ioresult = 0;
     }
     }
     return pp;
     } */

static int
matchpath_pos (int id, int *a_path)
{
  int i;
  int pp = 0;

  for (i = max_matchpath_length; i >= 1; i--)
    {
      if (a_path[i - 1] == id)
	pp = i;
    }
  return pp;
}


static int
matchpath_length (int *a_path)
{
  if (a_path[max_matchpath_length - 1] != 0)
    return max_matchpath_length;
  else
    return (matchpath_pos (0L, a_path) - 1);
}

static int
get_ndl_bond (int ba1, int ba2)
{
  int i;
  int b_id = 0;
  int FORLIM;

  if (ndl_n_bonds <= 0)
    return b_id;
  FORLIM = ndl_n_bonds;
  for (i = 1; i <= FORLIM; i++)
    {
      if (ndl_bond[i - 1].a1 == ba1 && ndl_bond[i - 1].a2 == ba2 ||
	  ndl_bond[i - 1].a1 == ba2 && ndl_bond[i - 1].a2 == ba1)
	b_id = i;
    }
  return b_id;
}

static int
ringcompare (int *rp1, int *rp2)
{
  int i, j;
  int rc = 0;
  int rs1, rs2;
  int n_common = 0;
  int max_cra;

  rs1 = path_length (rp1);
  rs2 = path_length (rp2);
  if (rs1 < rs2)
    max_cra = rs1;
  else
    max_cra = rs2;
  for (i = 0; i < rs1; i++)
    {
      for (j = 0; j < rs2; j++)
	{
	  if (rp1[i] == rp2[j])
	    n_common++;
	}
    }
  if (rs1 == rs2 && n_common == max_cra)
    {
      rc = 0;
      return rc;
    }
  if (n_common == 0)
    rc += 8;
  if (n_common < max_cra)
    {
      rc += 4;
      return rc;
    }
  if (rs1 < rs2)
    rc++;
  else
    rc += 2;
  return rc;
}


static boolean
rc_identical (int rc_int)
{
  if (rc_int == 0)
    return true;
  else
    return false;
}


static boolean
rc_1in2 (int rc_int)
{
  if (rc_int & 1)
    return true;
  else
    return false;
}


static boolean
rc_2in1 (int rc_int)
{
  rc_int /= 2;
  if (rc_int & 1)
    return true;
  else
    return false;
}


static boolean
rc_different (int rc_int)
{
  rc_int /= 4;
  if (rc_int & 1)
    return true;
  else
    return false;
}


static boolean
rc_independent (int rc_int)
{
  rc_int /= 8;
  if (rc_int & 1)
    return true;
  else
    return false;
}


static boolean
is_newring (int *n_path)
{
  int i, j;
  boolean nr = true;
  boolean same_ring;
  ringpath_type tmp_path;
  int rc_result;
  boolean found_ring;
  int pl;			/* new in v0.3 */
  int FORLIM;

  pl = path_length (n_path);	/* new in v0.3 */
  if (n_rings <= 0)
    return true;
  switch (ringsearch_mode)
    {

    case rs_sar:
      found_ring = false;
      i = 0;
      while (i < n_rings && !found_ring)
	{
	  i++;
	  if (pl != ringprop[i - 1].size)	/* compare only rings of same size */
	    continue;
	  same_ring = true;
	  for (j = 0; j < max_ringsize; j++)
	    {
	      if (ring[i - 1][j] != n_path[j])
		same_ring = false;
	    }
	  if (same_ring)
	    {
	      nr = false;
	      found_ring = true;
	    }
	}
      break;

    case rs_ssr:
      FORLIM = n_rings;
      for (i = 0; i < FORLIM; i++)
	{
	  for (j = 0; j < max_ringsize; j++)
	    tmp_path[j] = ring[i][j];
	  rc_result = ringcompare (n_path, tmp_path);
	  if (rc_identical (rc_result))
	    nr = false;
	  if (rc_1in2 (rc_result))
	    {
	      /* exchange existing ring by smaller one */
	      for (j = 0; j < max_ringsize; j++)
		ring[i][j] = n_path[j];
	      /* update ring property record (new in v0.3) */
	      ringprop[i].size = pl;
	      nr = false;
/* p2c: checkmol.pas, line 2841:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
	      /*$IFDEF debug */
	      /* debugoutput('replacing ring '+inttostr(i)+' by smaller one (ringsize: '+inttostr(path_length(n_path))+')'); */
	      /*$ENDIF */
	    }
	  if (rc_2in1 (rc_result))
	    {
	      /* new ring contains existing one, but is larger ==> discard */
	      nr = false;
	    }
	}
      break;
    }
  return nr;
}


static void
add_ring (int *n_path)
{
  /* new in v0.3: store rings in an ordered way (with ascending ring size) */
  int i;
  int j = 0;
  int k, s, pl;
/* p2c: checkmol.pas, line 2862:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
  /*$IFDEF debug */
  char dstr[256];
  int FORLIM;

  /*$ENDIF */
  pl = path_length (n_path);

  if (pl < 1)
    return;

  if (n_rings >= max_rings)
    {

/* p2c: checkmol.pas, line 2906:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
      /*$IFDEF debug 
         debugoutput ("max_rings exceeded!");
         /*$ENDIF */
      return;
    }
  n_rings++;

/* p2c: checkmol.pas, line 2871:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
  /*$IFDEF debug */
  /* dstr := '';
     for j := 1 to pl do dstr := dstr + inttostr((n_path[j])) + ' ';
     debugoutput('adding ring '+inttostr(n_rings)+':  '+dstr); */
  /*$ENDIF */
  if (n_rings > 1)
    {
      FORLIM = n_rings;
      for (i = 1; i < FORLIM; i++)
	{
	  s = ringprop[i - 1].size;
	  if (pl >= s)
	    j = i;
	}
    }
  j++;				/* the next position is ours */
  if (j < n_rings)
    {				/* push up the remaining rings by one position */
      for (k = n_rings; k > j; k--)
	{
	  ringprop[k - 1].size = ringprop[k - 2].size;
	  /*ringprop^[k].arom := ringprop^[(k-1)].arom; */
	  for (i = 0; i < max_ringsize; i++)
	    ring[k - 1][i] = ring[k - 2][i];
	}
    }
  ringprop[j - 1].size = pl;
  for (i = 0; i < max_ringsize; i++)
    {
      ring[j - 1][i] = n_path[i];
      /*inc(atom^[(n_path[i])].ring_count); */
    }
}


  /* static boolean is_ringpath (s_path) int *s_path;
     {
     boolean Result;
     int i, j;
     neighbor_rec nb;
     boolean rp = false;
     boolean new_atom;
     int a_last, pl;
     ringpath_type l_path;
     int FORLIM;

     /* p2c: checkmol.pas: Note: Eliminated unused assignment statement [338] 
     memset (nb, 0, sizeof (neighbor_rec));
     memset (l_path, 0, sizeof (ringpath_type));
     pl = path_length (s_path);
     if (pl < 1)
     {
     /* p2c: checkmol.pas, line 2928:
     * Note: Turbo Pascal conditional compilation directive was ignored [218] 
     /*$IFDEF debug 
     debugoutput ("Oops! Got zero-length s_path!");
     /*$ENDIF 
     return Result;
     }
     for (i = 0; i < pl; i++)
     l_path[i] = s_path[i];
     /* check if the last atom is a metal and stop if opt_metalrings is not set (v0.3) 
     if (opt_metalrings == false)
     {
     if (atom[l_path[pl - 1] - 1].metal)
     {
     /* p2c: checkmol.pas, line 2942:
     * Note: Turbo Pascal conditional compilation directive was ignored [218] 
     /*$IFDEF debug 
     debugoutput ("skipping metal in ring search");
     /*$ENDIF 
     return false;
     }
     }
     /* check if ring is already closed 
     if (pl > 2 && l_path[pl - 1] == l_path[0])
     {
     l_path[pl - 1] = 0;        /* remove last entry (redundant!) 
     order_ringpath (l_path);
     if (is_newring (l_path))
     {
     if (n_rings >= max_rings)
     {
     /* p2c: checkmol.pas, line 2958:
     * Note: Turbo Pascal conditional compilation directive was ignored [218] 
     /*$IFDEF debug 
     debugoutput ("maximum number of rings exceeded!");
     /*$ENDIF 
     return false;
     }
     add_ring (l_path);
     }
     /* p2c: checkmol.pas: Note: Eliminated unused assignment statement [338] 
     return true;
     }
     /* any other case: ring is not (yet) closed 
     a_last = l_path[pl - 1];
     get_neighbors (nb, a_last);
     if (atom[a_last - 1].neighbor_count <= 1)
     return false;
     if (n_rings >= max_rings)
     /* added in v0.2: check if max_rings is reached *
     {                          /* if ring is not closed, continue searching 
     return false;
     }
     FORLIM = atom[a_last - 1].neighbor_count;
     for (i = 0; i < FORLIM; i++)
     {
     new_atom = true;
     for (j = 1; j < pl; j++)
     {
     if (nb[i] == l_path[j])
     {                    /* v0.3k 
     new_atom = false;
     /* p2c: checkmol.pas, line 2982:
     * Warning: Expected a '(', found a semicolon [227] *
     /* p2c: checkmol.pas, line 2982:
     * Warning: Expected an expression, found a semicolon [227] 
     fflush (0);
     P_ioresult = 0;      /* v0.3k *
     }
     }

     /* added in v0.1a: check if max_rings not yet reached 
     /* added in v0.2:  limit ring size to max_vringsize instead of max_ringsize 
     if (new_atom && pl < max_vringsize && n_rings < max_rings)
     {
     l_path[pl] = nb[i];
     if (pl < max_ringsize - 1)   /* just to be sure 
     l_path[pl + 1] = 0;
     if (is_ringpath (l_path))
     rp = true;
     }
     }
     return rp;
     } */

static boolean
is_ringpath (int *s_path)
{
  int i, j;
  neighbor_rec nb;
  boolean rp = false;
  boolean new_atom;
  int a_last, pl;
  ringpath_type l_path;
  int FORLIM, pl_prev, pl_next, max_ringsize_dec;

/* p2c: checkmol.pas: Note: Eliminated unused assignment statement [338] */
  memset ((void *) nb, 0, sizeof (neighbor_rec));
  memset ((void *) l_path, 0, sizeof (ringpath_type));
  pl = path_length (s_path);
  if (pl < 1)
    {
/* p2c: checkmol.pas, line 2524:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */

      return false;
    }

  pl_prev = pl - 1;
  //memcpy (l_path, s_path, sizeof (ringpath_type));
  memcpy (l_path, s_path, pl * sizeof (int));

  /* check if the last atom is a metal and stop if opt_metalrings is not set (v0.3) */
  if (opt_metalrings == false)
    {
      if (atom[l_path[pl_prev] - 1].metal)
	{
/* p2c: checkmol.pas, line 2538:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
	  return false;
	}
    }
  /* check if ring is already closed */
  if (pl > 2 && l_path[pl_prev] == l_path[0])
    {
      l_path[pl_prev] = 0;	/* remove last entry (redundant!) */
      order_ringpath (l_path);
      if (is_newring (l_path))
	{
	  if (n_rings >= max_rings)
	    {
/* p2c: checkmol.pas, line 2554:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */

	      return false;
	    }
	  add_ring (l_path);
	}
/* p2c: checkmol.pas: Note: Eliminated unused assignment statement [338] */
      return true;
    }
  /* any other case: ring is not (yet) closed */
  a_last = l_path[pl_prev];
  get_neighbors (nb, a_last);
  if (atom[a_last - 1].neighbor_count <= 1)
    return false;
  if (n_rings >= max_rings)
    /* added in v0.2: check if max_rings is reached */
    {				/* if ring is not closed, continue searching */
      return false;
    }
  FORLIM = atom[a_last - 1].neighbor_count;
  pl_next = pl + 1;
  max_ringsize_dec = max_ringsize - 1;
  for (i = 0; i < FORLIM; i++)
    {
      new_atom = true;
      for (j = 1; j < pl; j++)
	{
	  if (nb[i] == l_path[j])
	    {
	      new_atom = false;
	      break;
	    }
	}
      /* added in v0.1a: check if max_rings not yet reached */
      /* added in v0.2:  limit ring size to max_vringsize instead of max_ringsize */
      if (new_atom && pl < max_vringsize && n_rings < max_rings)
	{
	  l_path[pl] = nb[i];
	  if (pl < max_ringsize_dec)	/* just to be sure */
	    l_path[pl_next] = 0;

	  if (max_ringpath_recursion_depth != 0
	      && ++recursion_depth > max_ringpath_recursion_depth)
	    {
#ifndef MAKE_SHARED_LIBRARY
	      if (opt_verbose)
#endif
		printf
		  ("Warning: max. number of ringpath recursions (%i) reached\n",
		   max_ringpath_recursion_depth);
	      n_rings = max_rings;
	      return false;
	    }

	  //printf("%i\n",recursion_depth);
	  //fflush(stdout);

	  if (is_ringpath (l_path))
	    rp = true;
	  //return true;
	}
    }
  return rp;
}

static boolean
is_ringbond (int b_id)
{
  int i, ra1, ra2;
  neighbor_rec nb;
  ringpath_type search_path;
  boolean rb = false;
  int FORLIM;

  recursion_depth = 0;

  ra1 = bond[b_id - 1].a1;
  ra2 = bond[b_id - 1].a2;
  memset (nb, 0, sizeof (neighbor_rec));
  memset (search_path, 0, sizeof (ringpath_type));
  get_neighbors (nb, ra2);
  if (atom[ra2 - 1].neighbor_count <= 1 || atom[ra1 - 1].neighbor_count <= 1)
    return false;
  search_path[0] = ra1;
  search_path[1] = ra2;
  FORLIM = atom[ra2 - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      if (nb[i] != ra1 && atom[nb[i] - 1].heavy)
	{
	  search_path[2] = nb[i];
	  if (is_ringpath (search_path))
	    rb = true;
	  //return true;
	}
    }
  return rb;
}


static void
chk_ringbonds ()
{
  int i, a1rc, a2rc;

  if (n_bonds < 1)
    return;

  for (i = 0; i < n_bonds; i++)
    {
      a1rc = atom[bond[i].a1 - 1].ring_count;
      a2rc = atom[bond[i].a2 - 1].ring_count;
      if (n_rings == 0 || a1rc < n_rings && a2rc < n_rings)
	{
	  is_ringbond (i + 1);
	  /*inc(bond^[i].ring_count); */
	}
    }
}


/* v0.3d: moved procedure update_ringcount a bit down */


static int
raw_hetbond_count (int a)
{
  /* new in v0.2j, ignores bond order */
  int i;
  neighbor_rec nb;
  str2 nb_el;
  int hbc = 0;
  int FORLIM;

  if (a <= 0 || a > n_atoms)
    return hbc;
  memset (nb, 0, sizeof (neighbor_rec));
  get_neighbors (nb, a);
  if (atom[a - 1].neighbor_count <= 0)
    return hbc;
  FORLIM = atom[a - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      strcpy (nb_el, atom[nb[i] - 1].element);
      if (strcmp (nb_el, "C ") && strcmp (nb_el, "A ") && strcmp (nb_el, "H ")	/* &&  strcmp (nb_el, "D ") */
	  && strcmp (nb_el, "LP") && strcmp (nb_el, "DU"))
	/* added 'D ' in v0.3n */
	hbc++;
    }
  return hbc;
}


static int
hetbond_count (int a)
{
  int i;
  neighbor_rec nb;
  str2 nb_el;
  float hbc = 0.0;
  int FORLIM;

  if (a <= 0 || a > n_atoms)
    return ((int) floor (hbc + 0.5));
  memset (nb, 0, sizeof (neighbor_rec));
  get_neighbors (nb, a);
  if (atom[a - 1].neighbor_count <= 0)
    return ((int) floor (hbc + 0.5));
  FORLIM = atom[a - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      strcpy (nb_el, atom[nb[i] - 1].element);
      if (strcmp (nb_el, "C ") && strcmp (nb_el, "A ") && strcmp (nb_el, "H ")	/*&& strcmp (nb_el, "D ") */
	  && strcmp (nb_el, "LP") && strcmp (nb_el, "DU"))
	{			/* added 'D ' in v0.3n */
	  if (bond[get_bond (a, nb[i]) - 1].btype == 'S')
	    hbc += 1.0;
	  if (bond[get_bond (a, nb[i]) - 1].btype == 'A')
	    hbc += 1.5;
	  if (bond[get_bond (a, nb[i]) - 1].btype == 'D')
	    hbc += 2.0;
	  if (bond[get_bond (a, nb[i]) - 1].btype == 'T')
	    hbc += 3.0;
	}
    }
  return ((int) floor (hbc + 0.5));
}


static int
hetatom_count (int a)
{
  int i;
  neighbor_rec nb;
  str2 nb_el;
  int hac = 0;
  int FORLIM;

  if (a <= 0 || a > n_atoms)
    return hac;
  memset (nb, 0, sizeof (neighbor_rec));
  get_neighbors (nb, a);
  if (atom[a - 1].neighbor_count <= 0)
    return hac;
  FORLIM = atom[a - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      strcpy (nb_el, atom[nb[i] - 1].element);
      if (strcmp (nb_el, "C ") && strcmp (nb_el, "H ")
	  /*&&  strcmp (nb_el, "D ") */  && strcmp (nb_el, "LP")
	  && strcmp (nb_el, "DU"))
	/* added 'D ' in v0.3n */
	hac++;
    }
  return hac;
}


static int
ndl_hetbond_count (int a)
{
  int i;
  neighbor_rec nb;
  str2 nb_el;
  float hbc = 0.0;
  int FORLIM;

  if (a <= 0 || a > n_atoms)
    return ((int) floor (hbc + 0.5));
  memset (nb, 0, sizeof (neighbor_rec));
  get_ndl_neighbors (nb, a);
  if (ndl_atom[a - 1].neighbor_count <= 0)
    return ((int) floor (hbc + 0.5));
  FORLIM = ndl_atom[a - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      strcpy (nb_el, ndl_atom[nb[i] - 1].element);
      if (strcmp (nb_el, "C ") && strcmp (nb_el, "H ")
	  /*&&  strcmp (nb_el, "D ") */  && strcmp (nb_el, "LP")
	  && strcmp (nb_el, "DU"))
	{			/* added 'D ' in v0.3n */
	  if (ndl_bond[get_ndl_bond (a, nb[i]) - 1].btype == 'S')
	    hbc += 1.0;
	  if (ndl_bond[get_ndl_bond (a, nb[i]) - 1].btype == 'A')
	    hbc += 1.5;
	  if (ndl_bond[get_ndl_bond (a, nb[i]) - 1].btype == 'D')
	    hbc += 2.0;
	  if (ndl_bond[get_ndl_bond (a, nb[i]) - 1].btype == 'T')
	    hbc += 3.0;
	}
    }
  return ((int) floor (hbc + 0.5));
}


static int
ndl_hetatom_count (int a)
{
  int i;
  neighbor_rec nb;
  str2 nb_el;
  int hac = 0;
  int FORLIM;

  if (a <= 0 || a > ndl_n_atoms)
    return hac;
  memset (nb, 0, sizeof (neighbor_rec));
  get_ndl_neighbors (nb, a);
  if (ndl_atom[a - 1].neighbor_count <= 0)
    return hac;
  FORLIM = ndl_atom[a - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {				/* note: query atoms like 'A' should be present only in the needle */
      strcpy (nb_el, ndl_atom[nb[i] - 1].element);
      if (strcmp (nb_el, "C ") && strcmp (nb_el, "A ") && strcmp (nb_el, "H ")	/*&&  strcmp (nb_el, "D ") */
	  && strcmp (nb_el, "LP") && strcmp (nb_el, "DU"))
	/* added 'D ' in v0.3n */
	hac++;
    }
  return hac;
}


static boolean
is_oxo_C (int id)
{
  boolean Result;
  int i;
  boolean r = false;
  neighbor_rec nb;
  int FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  if (id < 1 || id > n_atoms)
    return Result;
  get_neighbors (nb, id);
  if (strcmp (atom[id - 1].element, "C ") || atom[id - 1].neighbor_count <= 0)
    return false;
  FORLIM = atom[id - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      if (bond[get_bond (id, nb[i]) - 1].btype == 'D' &&
	  !strcmp (atom[nb[i] - 1].element, "O "))
	/* no N, amidines are different... */
	r = true;
      /* or
         (atom^[(nb[i])].element = 'S ')  or
         (atom^[(nb[i])].element = 'SE') */
    }
  return r;
}


static boolean
is_thioxo_C (int id)
{
  boolean Result;
  int i;
  boolean r = false;
  neighbor_rec nb;
  int FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  if (id < 1 || id > n_atoms)
    return Result;
  get_neighbors (nb, id);
  if (strcmp (atom[id - 1].element, "C ") || atom[id - 1].neighbor_count <= 0)
    return false;
  FORLIM = atom[id - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      if (bond[get_bond (id, nb[i]) - 1].btype == 'D' &&
	  (!strcmp (atom[nb[i] - 1].element, "S ") ||
	   !strcmp (atom[nb[i] - 1].element, "SE")))
	/* no N, amidines are different... */
	r = true;
    }
  return r;
}


static boolean
is_imino_C (int id)
{
  boolean Result;
  int i;
  boolean r = false;
  neighbor_rec nb;
  int FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  if (id < 1 || id > n_atoms)
    return Result;
  get_neighbors (nb, id);
  if (strcmp (atom[id - 1].element, "C ") || atom[id - 1].neighbor_count <= 0)
    return false;
  FORLIM = atom[id - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      if (bond[get_bond (id, nb[i]) - 1].btype == 'D' &&
	  !strcmp (atom[nb[i] - 1].element, "N "))
	r = true;
    }
  return r;
}


static boolean
is_true_imino_C (int id)
{
  boolean Result;
  int i;
  boolean r = true;
  neighbor_rec nb;
  str2 nb_el;
  int a_n = 0;
  int b;			/* v0.3j */
  int FORLIM;

/* p2c: checkmol.pas: Note: Eliminated unused assignment statement [338] */
  memset (nb, 0, sizeof (neighbor_rec));
  if (id < 1 || id > n_atoms)
    return Result;
  get_neighbors (nb, id);
  if (strcmp (atom[id - 1].element, "C ") || atom[id - 1].neighbor_count <= 0)
    return false;
  FORLIM = atom[id - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      b = get_bond (id, nb[i]);	/* v0.3j */
      if (bond[b - 1].btype == 'D' && bond[b - 1].arom == false &&
	  !strcmp (atom[nb[i] - 1].element, "N "))
	/* v0.3j */
	a_n = nb[i];
    }
  if (a_n <= 0)
    return false;
  memset (nb, 0, sizeof (neighbor_rec));
  get_neighbors (nb, a_n);
  FORLIM = atom[a_n - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      strcpy (nb_el, atom[nb[i] - 1].element);
      if (strcmp (nb_el, "C ") && strcmp (nb_el, "H ")
	  /*&&  strcmp (nb_el, "D ") */ )
	/* v0.3n: D */
	r = false;
    }
  return r;
}


static boolean
is_true_exocyclic_imino_C (int id, int r_id)
{
  /* v0.3j */
  boolean Result;
  int i, j;
  boolean r = false;
  neighbor_rec nb;
  ringpath_type testring;
  int ring_size, b, FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  if (id < 1 || id > n_atoms)
    return Result;
  get_neighbors (nb, id);
  memset (testring, 0, sizeof (ringpath_type));
  ring_size = ringprop[r_id - 1].size;	/* v0.3j */
  for (j = 0; j < ring_size; j++)	/* v0.3j */
    testring[j] = ring[r_id - 1][j];
  if (strcmp (atom[id - 1].element, "C ") || atom[id - 1].neighbor_count <= 0)
    return false;
  FORLIM = atom[id - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      b = get_bond (id, nb[i]);
      if (bond[b - 1].btype == 'D' && bond[b - 1].arom == false &&
	  !strcmp (atom[nb[i] - 1].element, "N "))
	{
	  r = true;
	  for (j = 0; j < ring_size; j++)
	    {
	      if (nb[i] == ring[r_id - 1][j])
		r = false;
	    }
	}
    }
  return r;
}


static boolean
is_exocyclic_imino_C (int id, int r_id)
{
  boolean Result;
  int i, j;
  boolean r = false;
  neighbor_rec nb;
  ringpath_type testring;
  int ring_size, FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  if (id < 1 || id > n_atoms)
    return Result;
  get_neighbors (nb, id);
  memset (testring, 0, sizeof (ringpath_type));
  ring_size = ringprop[r_id - 1].size;	/* v0.3j */
  for (j = 0; j < ring_size; j++)	/* v0.3j */
    testring[j] = ring[r_id - 1][j];
  if (strcmp (atom[id - 1].element, "C ") || atom[id - 1].neighbor_count <= 0)
    return false;
  FORLIM = atom[id - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      if (bond[get_bond (id, nb[i]) - 1].btype == 'D' &&
	  !strcmp (atom[nb[i] - 1].element, "N "))
	{
	  r = true;
	  for (j = 0; j < ring_size; j++)
	    {
	      if (nb[i] == ring[r_id - 1][j])
		r = false;
	    }
	}
    }
  return r;
}


static int
find_exocyclic_methylene_C (int id, int r_id)
{
  /* renamed and rewritten in v0.3j */
  int i, j;
  int r = 0;
  neighbor_rec nb;
  ringpath_type testring;
  int ring_size, FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  if (id < 1 || id > n_atoms)
    return 0;
  get_neighbors (nb, id);
  memset (testring, 0, sizeof (ringpath_type));
  ring_size = ringprop[r_id - 1].size;	/* v0.3j */
  for (j = 0; j < ring_size; j++)	/* v0.3j */
    testring[j] = ring[r_id - 1][j];
  if (strcmp (atom[id - 1].element, "C ") || atom[id - 1].neighbor_count <= 0)
    return r;
  FORLIM = atom[id - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      if (bond[get_bond (id, nb[i]) - 1].btype == 'D' &&
	  !strcmp (atom[nb[i] - 1].element, "C "))
	{
	  r = nb[i];
	  for (j = 0; j < ring_size; j++)
	    {
	      if (nb[i] == ring[r_id - 1][j])
		r = 0;
	    }
	}
    }
  return r;
}


static boolean
is_hydroxy (int a_view, int a_ref)
{
  boolean r = false;

  if (atom[a_view - 1].
      heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S'))
    {
      if (!strcmp (atom[a_ref - 1].atype, "O3 ") &&
	  atom[a_ref - 1].neighbor_count == 1)
	r = true;
    }
  return r;
}


static boolean
is_sulfanyl (int a_view, int a_ref)
{
  boolean r = false;

  if (atom[a_view - 1].
      heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S'))
    {
      if (!strcmp (atom[a_ref - 1].atype, "S3 ") &&
	  atom[a_ref - 1].neighbor_count == 1)
	r = true;
    }
  return r;
}


static boolean
is_amino (int a_view, int a_ref)
{
  boolean r = false;

  if (atom[a_view - 1].
      heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S'))
    {
      if ((!strcmp (atom[a_ref - 1].atype, "N3 ") ||
	   !strcmp (atom[a_ref - 1].atype, "N3+")) &&
	  atom[a_ref - 1].neighbor_count == 1)
	r = true;
    }
  return r;
}


static boolean
is_alkyl (int a_view, int a_ref)
{
  int i;
  boolean r = false;
  neighbor_rec nb;
  str2 nb_el;
  int het_count = 0;
  int FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S'))
      || strcmp (atom[a_ref - 1].atype, "C3 ")
      || atom[a_ref - 1].arom != false)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  FORLIM = atom[a_ref - 1].neighbor_count - 2;
  for (i = 0; i <= FORLIM; i++)
    {
      strcpy (nb_el, atom[nb[i] - 1].element);
      if (strcmp (nb_el, "C ") && strcmp (nb_el, "H ")
	  /*&&  strcmp (nb_el, "D ") */  && strcmp (nb_el, "DU")
	  && strcmp (nb_el, "LP"))
	/* added 'D ' in v0.3n */
	het_count++;
    }
  if (het_count <= 1)		/* we consider (e.g.) alkoxyalkyl groups as alkyl */
    r = true;
  return r;
}


static boolean
is_true_alkyl (int a_view, int a_ref)
{
  int i;
  boolean r = false;
  neighbor_rec nb;
  str2 nb_el;
  int het_count = 0;
  int FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S'))
      || strcmp (atom[a_ref - 1].atype, "C3 ")
      || atom[a_ref - 1].arom != false)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  FORLIM = atom[a_ref - 1].neighbor_count - 2;
  for (i = 0; i <= FORLIM; i++)
    {
      strcpy (nb_el, atom[nb[i] - 1].element);
      if (strcmp (nb_el, "C ") && strcmp (nb_el, "H ")
	  /*&&  strcmp (nb_el, "D ") */  && strcmp (nb_el, "DU"))
	/* added 'D ' in v0.3n */
	het_count++;
    }
  if (het_count == 0)		/* */
    r = true;
  return r;
}


static boolean
is_alkenyl (int a_view, int a_ref)
{
  /* new in v0.3j */
  int i;
  boolean r = false;
  neighbor_rec nb;
  str2 nb_el;
  str3 nb_at;
  int c2_count = 0, het_count = 0;
  int FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S'))
      || strcmp (atom[a_ref - 1].atype, "C2 ")
      || atom[a_ref - 1].arom != false)
    {
      return false;
    }				/* v0.3k: changed c2_count = 1 into c2_count >= 1 */
  get_nextneighbors (nb, a_ref, a_view);
  FORLIM = atom[a_ref - 1].neighbor_count - 2;
  for (i = 0; i <= FORLIM; i++)
    {
      strcpy (nb_el, atom[nb[i] - 1].element);
      strcpy (nb_at, atom[nb[i] - 1].atype);
      if (strcmp (nb_el, "C ") && strcmp (nb_el, "H ")
	  /*&&  strcmp (nb_el, "D ") */  && strcmp (nb_el, "DU")
	  && strcmp (nb_el, "LP"))
	/* added 'D ' in v0.3n */
	het_count++;
      if (!strcmp (nb_at, "C2 "))
	c2_count++;
    }
  if (c2_count >= 1 && het_count <= 1)
    /* we consider (e.g.) alkoxyalkenyl groups as alkenyl */
    r = true;
  return r;
}


static boolean
is_alkynyl (int a_view, int a_ref)
{
  /* new in v0.3j */
  int i;
  boolean r = false;
  neighbor_rec nb;
  str3 nb_at;
  int c1_count = 0;
  int FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S'))
      || strcmp (atom[a_ref - 1].atype, "C1 ")
      || atom[a_ref - 1].arom != false)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  FORLIM = atom[a_ref - 1].neighbor_count - 2;
  for (i = 0; i <= FORLIM; i++)
    {
      strcpy (nb_at, atom[nb[i] - 1].atype);
      if (!strcmp (nb_at, "C1 "))
	c1_count++;
    }
  if (c1_count == 1)
    r = true;
  return r;
}


static boolean
is_aryl (int a_view, int a_ref)
{
  boolean r = false;

  if ((atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S'))
      && !strcmp (atom[a_ref - 1].element, "C ")
      && atom[a_ref - 1].arom == true)
    r = true;
  return r;
}


static boolean
is_alkoxy (int a_view, int a_ref)
{
  boolean r = false;
  neighbor_rec nb;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].atype, "O3 ")
      || atom[a_ref - 1].neighbor_count != 2)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  if (is_alkyl (a_ref, nb[0]))
    r = true;
  return r;
}


static boolean
is_siloxy (int a_view, int a_ref)
{
  boolean r = false;
  neighbor_rec nb;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].atype, "O3 ")
      || atom[a_ref - 1].neighbor_count != 2)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  if (!strcmp (atom[nb[0] - 1].element, "SI"))
    r = true;
  return r;
}


static boolean
is_true_alkoxy (int a_view, int a_ref)
{
  boolean r = false;
  neighbor_rec nb;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].atype, "O3 ")
      || atom[a_ref - 1].neighbor_count != 2)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  if (is_true_alkyl (a_ref, nb[0]))
    r = true;
  return r;
}


static boolean
is_aryloxy (int a_view, int a_ref)
{
  boolean r = false;
  neighbor_rec nb;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].atype, "O3 ")
      || atom[a_ref - 1].neighbor_count != 2)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  if (is_aryl (a_ref, nb[0]))
    r = true;
  return r;
}


static boolean
is_alkenyloxy (int a_view, int a_ref)
{
  /* v0.3j */
  boolean r = false;
  neighbor_rec nb;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].atype, "O3 ")
      || atom[a_ref - 1].neighbor_count != 2)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  if (is_alkenyl (a_ref, nb[0]))
    r = true;
  return r;
}


static boolean
is_alkynyloxy (int a_view, int a_ref)
{
  /* v0.3j */
  boolean r = false;
  neighbor_rec nb;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].atype, "O3 ")
      || atom[a_ref - 1].neighbor_count != 2)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  if (is_alkynyl (a_ref, nb[0]))
    r = true;
  return r;
}


static boolean
is_alkylsulfanyl (int a_view, int a_ref)
{
  boolean r = false;
  neighbor_rec nb;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].atype, "S3 ")
      || atom[a_ref - 1].neighbor_count != 2)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  if (is_alkyl (a_ref, nb[0]))
    r = true;
  return r;
}


static boolean
is_true_alkylsulfanyl (int a_view, int a_ref)
{
  boolean r = false;
  neighbor_rec nb;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].atype, "S3 ")
      || atom[a_ref - 1].neighbor_count != 2)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  if (is_true_alkyl (a_ref, nb[0]))
    r = true;
  return r;
}


static boolean
is_arylsulfanyl (int a_view, int a_ref)
{
  boolean r = false;
  neighbor_rec nb;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].atype, "S3 ")
      || atom[a_ref - 1].neighbor_count != 2)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  if (is_aryl (a_ref, nb[0]))
    r = true;
  return r;
}


static boolean
is_alkenylsulfanyl (a_view, a_ref)
     int a_view, a_ref;
{
  /* v0.3j */
  boolean r = false;
  neighbor_rec nb;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].atype, "S3 ")
      || atom[a_ref - 1].neighbor_count != 2)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  if (is_alkenyl (a_ref, nb[0]))
    r = true;
  return r;
}


static boolean
is_alkynylsulfanyl (a_view, a_ref)
     int a_view, a_ref;
{
  /* v0.3j */
  boolean r = false;
  neighbor_rec nb;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].atype, "S3 ")
      || atom[a_ref - 1].neighbor_count != 2)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  if (is_alkynyl (a_ref, nb[0]))
    r = true;
  return r;
}


static boolean
is_alkylamino (a_view, a_ref)
     int a_view, a_ref;
{
  boolean r = false;
  neighbor_rec nb;
  int alkyl_count = 0;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].element, "N ")
      || atom[a_ref - 1].neighbor_count != 2)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  if (is_alkyl (a_ref, nb[0]))
    alkyl_count++;
  if (alkyl_count == 1)
    r = true;
  return r;
}


static boolean
is_dialkylamino (a_view, a_ref)
     int a_view, a_ref;
{
  int i;
  boolean r = false;
  neighbor_rec nb;
  int alkyl_count = 0;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].element, "N ")
      || atom[a_ref - 1].neighbor_count != 3)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  for (i = 0; i <= 1; i++)
    {
      if (is_alkyl (a_ref, nb[i]))
	alkyl_count++;
    }
  if (alkyl_count == 2)
    r = true;
  return r;
}


static boolean
is_arylamino (a_view, a_ref)
     int a_view, a_ref;
{
  boolean r = false;
  neighbor_rec nb;
  int aryl_count = 0;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].element, "N ")
      || atom[a_ref - 1].neighbor_count != 2)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  if (is_aryl (a_ref, nb[0]))
    aryl_count++;
  if (aryl_count == 1)
    r = true;
  return r;
}


static boolean
is_diarylamino (a_view, a_ref)
     int a_view, a_ref;
{
  int i;
  boolean r = false;
  neighbor_rec nb;
  int aryl_count = 0;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].element, "N ")
      || atom[a_ref - 1].neighbor_count != 3)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  for (i = 0; i <= 1; i++)
    {
      if (is_aryl (a_ref, nb[i]))
	aryl_count++;
    }
  if (aryl_count == 2)
    r = true;
  return r;
}


static boolean
is_alkylarylamino (a_view, a_ref)
     int a_view, a_ref;
{
  int i;
  boolean r = false;
  neighbor_rec nb;
  int alkyl_count = 0, aryl_count = 0;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].element, "N ")
      || atom[a_ref - 1].neighbor_count != 3)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  for (i = 0; i <= 1; i++)
    {
      if (is_alkyl (a_ref, nb[i]))
	alkyl_count++;
      if (is_aryl (a_ref, nb[i]))
	aryl_count++;
    }
  if (alkyl_count == 1 && aryl_count == 1)
    r = true;
  return r;
}


static boolean
is_C_monosubst_amino (a_view, a_ref)
     int a_view, a_ref;
{
  /* new in v0.3j */
  boolean r = false;
  neighbor_rec nb;
  int c_count = 0;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].atype, "N3 ")
      && strcmp (atom[a_ref - 1].atype, "NAM")
      || atom[a_ref - 1].neighbor_count != 2)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  if (!strcmp (atom[nb[0] - 1].element, "C "))
    c_count++;
  if (c_count == 1)
    r = true;
  return r;
}


static boolean
is_C_disubst_amino (a_view, a_ref)
     int a_view, a_ref;
{
  /* new in v0.3j */
  int i;
  boolean r = false;
  neighbor_rec nb;
  int b;
  int c_count = 0;

  b = get_bond (a_view, a_ref);
  memset (nb, 0, sizeof (neighbor_rec));
  if (!(atom[a_view - 1].heavy && bond[b - 1].btype == 'S' &&
	bond[b - 1].arom == false))
    return false;
  if (strcmp (atom[a_ref - 1].atype, "N3 ")
      && strcmp (atom[a_ref - 1].atype, "NAM")
      || atom[a_ref - 1].neighbor_count != 3)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  for (i = 0; i <= 1; i++)
    {
      if (!strcmp (atom[nb[i] - 1].element, "C "))
	c_count++;
    }
  if (c_count == 2)
    r = true;
  return r;
}


static boolean
is_subst_amino (a_view, a_ref)
     int a_view, a_ref;
{
  boolean r = false;

  if (is_amino (a_view, a_ref) || is_alkylamino (a_view, a_ref) |
      is_arylamino (a_view, a_ref) || is_dialkylamino (a_view, a_ref) |
      is_alkylarylamino (a_view, a_ref) || is_diarylamino (a_view, a_ref))
    r = true;
  return r;
}


static boolean
is_true_alkylamino (a_view, a_ref)
     int a_view, a_ref;
{
  boolean r = false;
  neighbor_rec nb;
  int alkyl_count = 0;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].atype, "N3 ")
      && strcmp (atom[a_ref - 1].atype, "N3+")
      || atom[a_ref - 1].neighbor_count != 2)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  if (is_true_alkyl (a_ref, nb[0]))
    alkyl_count++;
  if (alkyl_count == 1)
    r = true;
  return r;
}


static boolean
is_true_dialkylamino (a_view, a_ref)
     int a_view, a_ref;
{
  int i;
  boolean r = false;
  neighbor_rec nb;
  int alkyl_count = 0;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].atype, "N3 ")
      && strcmp (atom[a_ref - 1].atype, "N3+")
      || atom[a_ref - 1].neighbor_count != 3)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  for (i = 0; i <= 1; i++)
    {
      if (is_true_alkyl (a_ref, nb[i]))
	alkyl_count++;
    }
  if (alkyl_count == 2)
    r = true;
  return r;
}


static boolean
is_true_alkylarylamino (a_view, a_ref)
     int a_view, a_ref;
{
  int i;
  boolean r = false;
  neighbor_rec nb;
  int alkyl_count = 0, aryl_count = 0;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].atype, "N3 ")
      && strcmp (atom[a_ref - 1].atype, "N3+")
      || atom[a_ref - 1].neighbor_count != 3)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  for (i = 0; i <= 1; i++)
    {
      if (is_true_alkyl (a_ref, nb[i]))
	alkyl_count++;
      if (is_aryl (a_ref, nb[i]))
	aryl_count++;
    }
  if (alkyl_count == 1 && aryl_count == 1)
    r = true;
  return r;
}


static boolean
is_hydroxylamino (a_view, a_ref)
     int a_view, a_ref;
{
  int i;
  boolean r = false;
  neighbor_rec nb;
  int oh_count = 0, het_count = 0;	/* v0.3k */
  str2 nb_el;			/* v0.3k */
  int FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].element, "N ")
      || atom[a_ref - 1].neighbor_count < 2)
    /* v0.3c */
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  FORLIM = atom[a_ref - 1].neighbor_count - 2;
  for (i = 0; i <= FORLIM; i++)
    {				/* v0.3c */
      if (is_hydroxy (a_ref, nb[i]))
	oh_count++;
      strcpy (nb_el, atom[nb[i] - 1].element);	/* v0.3k */
      if (strcmp (nb_el, "C ") && strcmp (nb_el, "H ")
	  /*&&  strcmp (nb_el, "D ") */  && strcmp (nb_el, "DU")
	  && strcmp (nb_el, "LP"))
	/* v0.3k */
	het_count++;
      /* v0.3n: D */
    }
  if (oh_count == 1 && het_count == 1)
    r = true;
  return r;
}


static boolean
is_nitro (a_view, a_ref)
     int a_view, a_ref;
{
  int i;
  boolean r = false;
  neighbor_rec nb;
  int o_count = 0, bond_count = 0;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].element, "N ")
      || atom[a_ref - 1].neighbor_count != 3)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  for (i = 0; i <= 1; i++)
    {
      if (!strcmp (atom[nb[i] - 1].element, "O "))
	o_count++;
      if (bond[get_bond (a_ref, nb[i]) - 1].btype == 'S')
	bond_count++;
      if (bond[get_bond (a_ref, nb[i]) - 1].btype == 'D')
	bond_count += 2;
    }
  if (o_count == 2 && bond_count >= 3)
    r = true;
  return r;
}


static boolean
is_azido (a_view, a_ref)
     int a_view, a_ref;
{
  boolean r = false;
  neighbor_rec nb;
  int bond_count = 0, n1 = 0, n2 = 0, n3 = 0;

  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].element, "N ")
      || atom[a_ref - 1].neighbor_count != 2)
    return false;
  n1 = a_ref;
  memset (nb, 0, sizeof (neighbor_rec));
  get_nextneighbors (nb, n1, a_view);
  if (!strcmp (atom[nb[0] - 1].element, "N "))
    {
      n2 = nb[0];
      if (bond[get_bond (n1, n2) - 1].btype == 'S')
	bond_count++;
      if (bond[get_bond (n1, n2) - 1].btype == 'D')
	bond_count += 2;
      if (bond[get_bond (n1, n2) - 1].btype == 'T')
	bond_count += 3;
    }
  if (n2 > 0 && atom[n2 - 1].neighbor_count == 2)
    {
      memset (nb, 0, sizeof (neighbor_rec));
      get_nextneighbors (nb, n2, n1);
      if (!strcmp (atom[nb[0] - 1].element, "N "))
	{
	  n3 = nb[0];
	  if (bond[get_bond (n2, n3) - 1].btype == 'S')
	    bond_count++;
	  if (bond[get_bond (n2, n3) - 1].btype == 'D')
	    bond_count += 2;
	  if (bond[get_bond (n2, n3) - 1].btype == 'T')
	    bond_count += 3;
	}
    }
  if (n1 > 0 && n2 > 0 && n3 > 0 && atom[n3 - 1].neighbor_count == 1 &&
      bond_count > 3)
    r = true;
  return r;
}


static boolean
is_diazonium (a_view, a_ref)
     int a_view, a_ref;
{
  boolean r = false;
  neighbor_rec nb;
  int bond_count = 0, chg_count = 0, n1 = 0, n2 = 0;

  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].element, "N ")
      || atom[a_ref - 1].neighbor_count != 2)
    return false;
  n1 = a_ref;
  chg_count = atom[n1 - 1].formal_charge;
  memset (nb, 0, sizeof (neighbor_rec));
  get_nextneighbors (nb, n1, a_view);
  if (!strcmp (atom[nb[0] - 1].element, "N "))
    {
      n2 = nb[0];
      chg_count += atom[n2 - 1].formal_charge;
      if (bond[get_bond (n1, n2) - 1].btype == 'S')
	bond_count++;
      if (bond[get_bond (n1, n2) - 1].btype == 'D')
	bond_count += 2;
      if (bond[get_bond (n1, n2) - 1].btype == 'T')
	bond_count += 3;
    }
  if (n1 > 0 && n2 > 0 && atom[n2 - 1].neighbor_count == 1
      && bond_count >= 2 && chg_count > 0)
    r = true;
  return r;
}


static boolean
is_hydroximino_C (id)
     int id;
{
  boolean Result;
  int i;
  boolean r = false;
  neighbor_rec nb;
  int a_het = 0;
  int FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  if (id < 1 || id > n_atoms)
    return Result;
  get_neighbors (nb, id);
  if (strcmp (atom[id - 1].element, "C ") || atom[id - 1].neighbor_count <= 0)
    return false;
  FORLIM = atom[id - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      if ((bond[get_bond (id, nb[i]) - 1].btype == 'D' &&
	   !strcmp (atom[nb[i] - 1].element,
		    "N ")) && (hetbond_count (nb[i]) == 3))
	a_het = nb[i];
    }
  if (a_het <= 0)
    return false;
  memset (nb, 0, sizeof (neighbor_rec));
  get_neighbors (nb, a_het);
  if (strcmp (atom[a_het - 1].element, "N ")
      || atom[a_het - 1].neighbor_count <= 0)
    return false;
  FORLIM = atom[a_het - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      if (is_hydroxy (a_het, nb[i]))
	r = true;
    }
  return r;
}


static boolean
is_hydrazono_C (id)
     int id;
{
  boolean Result;
  int i;
  boolean r = false;
  neighbor_rec nb;
  int a_het = 0;
  int FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  if (id < 1 || id > n_atoms)
    return Result;
  get_neighbors (nb, id);
  if (strcmp (atom[id - 1].element, "C ") || atom[id - 1].neighbor_count <= 0)
    return false;
  FORLIM = atom[id - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      if (bond[get_bond (id, nb[i]) - 1].btype == 'D' &&
	  !strcmp (atom[nb[i] - 1].element, "N "))
	{
	  /* and
	     (hetbond_count(nb[i]) = 3)  */
	  a_het = nb[i];
	}
    }
  if (a_het <= 0)
    return false;
  memset (nb, 0, sizeof (neighbor_rec));
  get_neighbors (nb, a_het);
  if (strcmp (atom[a_het - 1].element, "N ")
      || atom[a_het - 1].neighbor_count <= 0)
    return false;
  FORLIM = atom[a_het - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      if (is_amino (a_het, nb[i]) || is_alkylamino (a_het, nb[i]) |
	  is_alkylarylamino (a_het, nb[i]) || is_arylamino (a_het, nb[i]) |
	  is_dialkylamino (a_het, nb[i]) || is_diarylamino (a_het, nb[i]))
	r = true;
    }
  return r;
}


static boolean
is_alkoxycarbonyl (a_view, a_ref)
     int a_view, a_ref;
{
  int i;
  boolean r = false;
  neighbor_rec nb;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (!(is_oxo_C (a_ref) && atom[a_ref - 1].neighbor_count == 3))
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  for (i = 0; i <= 1; i++)
    {
      if (is_alkoxy (a_ref, nb[i]))
	r = true;
    }
  return r;
}


static boolean
is_aryloxycarbonyl (a_view, a_ref)
     int a_view, a_ref;
{
  int i;
  boolean r = false;
  neighbor_rec nb;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (!(is_oxo_C (a_ref) && atom[a_ref - 1].neighbor_count == 3))
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  for (i = 0; i <= 1; i++)
    {
      if (is_aryloxy (a_ref, nb[i]))
	r = true;
    }
  return r;
}


static boolean
is_carbamoyl (a_view, a_ref)
     int a_view, a_ref;
{
  int i;
  boolean r = false;
  neighbor_rec nb;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (!(is_oxo_C (a_ref) && atom[a_ref - 1].neighbor_count == 3))
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  for (i = 0; i <= 1; i++)
    {
      if (!strcmp (atom[nb[i] - 1].atype, "N3 ") ||
	  !strcmp (atom[nb[i] - 1].atype, "NAM"))
	r = true;
    }
  return r;
}


static boolean
is_alkoxythiocarbonyl (a_view, a_ref)
     int a_view, a_ref;
{
  int i;
  boolean r = false;
  neighbor_rec nb;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (!(is_thioxo_C (a_ref) && atom[a_ref - 1].neighbor_count == 3))
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  for (i = 0; i <= 1; i++)
    {
      if (is_alkoxy (a_ref, nb[i]))
	r = true;
    }
  return r;
}


static boolean
is_aryloxythiocarbonyl (a_view, a_ref)
     int a_view, a_ref;
{
  int i;
  boolean r = false;
  neighbor_rec nb;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (!(is_thioxo_C (a_ref) && atom[a_ref - 1].neighbor_count == 3))
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  for (i = 0; i <= 1; i++)
    {
      if (is_aryloxy (a_ref, nb[i]))
	r = true;
    }
  return r;
}


static boolean
is_thiocarbamoyl (a_view, a_ref)
     int a_view, a_ref;
{
  int i;
  boolean r = false;
  neighbor_rec nb;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (!(is_thioxo_C (a_ref) && atom[a_ref - 1].neighbor_count == 3))
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  for (i = 0; i <= 1; i++)
    {
      if (!strcmp (atom[nb[i] - 1].atype, "N3 ") ||
	  !strcmp (atom[nb[i] - 1].atype, "NAM"))
	r = true;
    }
  return r;
}


static boolean
is_alkanoyl (a_view, a_ref)
     int a_view, a_ref;
{
  int i;
  boolean r = false;
  neighbor_rec nb;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (!(is_oxo_C (a_ref) && atom[a_ref - 1].neighbor_count == 3))
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  for (i = 0; i <= 1; i++)
    {
      if (is_alkyl (a_ref, nb[i]))
	r = true;
    }
  return r;
}


static boolean
is_aroyl (a_view, a_ref)
     int a_view, a_ref;
{
  int i;
  boolean r = false;
  neighbor_rec nb;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (!(is_oxo_C (a_ref) && atom[a_ref - 1].neighbor_count == 3))
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  for (i = 0; i <= 1; i++)
    {
      if (is_aryl (a_ref, nb[i]))
	r = true;
    }
  return r;
}


static boolean
is_acyl (a_view, a_ref)
     int a_view, a_ref;
{
  boolean r = false;

  if (is_alkanoyl (a_view, a_ref) || is_aroyl (a_view, a_ref))
    r = true;
  return r;
}


static boolean
is_acyl_gen (a_view, a_ref)
     int a_view, a_ref;
{
  /* new in v0.3j */
  boolean r = false;

  if (is_oxo_C (a_ref))
    r = true;
  return r;
}


static boolean
is_acylamino (a_view, a_ref)
     int a_view, a_ref;
{
  boolean r = false;
  neighbor_rec nb;
  int acyl_count = 0;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].element, "N ")
      || atom[a_ref - 1].neighbor_count != 2)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  if (is_acyl (a_ref, nb[0]))
    acyl_count++;
  if (acyl_count == 1)
    r = true;
  return r;
}


static boolean
is_subst_acylamino (a_view, a_ref)
     int a_view, a_ref;
{
  /* may be substituted _or_ unsubstituted acylamino group! */
  int i;
  boolean r = false;
  neighbor_rec nb;
  int acyl_count = 0;
  int FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].element, "N ")
      || atom[a_ref - 1].neighbor_count < 2)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  FORLIM = atom[a_ref - 1].neighbor_count - 2;
  for (i = 0; i <= FORLIM; i++)
    {
      if (is_acyl_gen (a_ref, nb[i]))	/* v0.3j */
	acyl_count++;
    }
  if (acyl_count > 0)
    r = true;
  return r;
}


static boolean
is_hydrazino (a_view, a_ref)
     int a_view, a_ref;
{
  int i;
  boolean r = false;
  neighbor_rec nb;
  int nr_count = 0;
  int FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].element, "N ")
      || atom[a_ref - 1].neighbor_count < 2)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  FORLIM = atom[a_ref - 1].neighbor_count - 2;
  for (i = 0; i <= FORLIM; i++)
    {				/* fixed in v0.3c */
      if (is_amino (a_ref, nb[i]) || is_subst_amino (a_ref, nb[i]))
	nr_count++;
    }
  if (nr_count == 1)
    r = true;
  return r;
}


static boolean
is_nitroso (a_view, a_ref)
     int a_view, a_ref;
{
  /* new in v0.3j */
  boolean r = false;
  neighbor_rec nb;
  int o_count = 0;
  int a2;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].element, "N ")
      || atom[a_ref - 1].neighbor_count != 2)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  a2 = nb[0];
  if ((strcmp (atom[a2 - 1].element, "O ") == 0) &
      (bond[get_bond (a_ref, a2) - 1].btype == 'D'))
    o_count++;
  if (o_count == 1)
    r = true;
  return r;
}


static boolean
is_subst_hydrazino (a_view, a_ref)
     int a_view, a_ref;
{
  int i;
  boolean r = false;
  neighbor_rec nb;
  int nr_count = 0;
  int a2, FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  if (!
      (atom[a_view - 1].
       heavy && (bond[get_bond (a_view, a_ref) - 1].btype == 'S')))
    return false;
  if (strcmp (atom[a_ref - 1].element, "N ")
      || atom[a_ref - 1].neighbor_count < 2)
    return false;
  get_nextneighbors (nb, a_ref, a_view);
  FORLIM = atom[a_ref - 1].neighbor_count - 2;
  for (i = 0; i <= FORLIM; i++)
    {
      a2 = nb[i];
      if ((strcmp (atom[a2 - 1].element, "N ") ==
	   0) && (!is_nitroso (a_ref, a2)))
	/* v0.3j */
	nr_count++;
    }
  if (nr_count == 1)
    r = true;
  return r;
}


static boolean
is_cyano (a_view, a_ref)
     int a_view, a_ref;
{
  boolean r = false;

  if (((strcmp (atom[a_view - 1].atype, "C1 ") == 0) &
       (bond[get_bond (a_view, a_ref) - 1].btype == 'T')) &&
      !strcmp (atom[a_ref - 1].atype, "N1 ") &&
      atom[a_ref - 1].neighbor_count == 1)
    r = true;
  return r;
}


static boolean
is_cyano_c (a_ref)
     int a_ref;
{
  int i;
  boolean r = false;
  neighbor_rec nb;
  int FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  if (strcmp (atom[a_ref - 1].atype, "C1 ")
      || atom[a_ref - 1].neighbor_count <= 0)
    return false;
  get_neighbors (nb, a_ref);
  FORLIM = atom[a_ref - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      if (is_cyano (a_ref, nb[i]))
	r = true;
    }
  return r;
}


static boolean
is_nitrile (a_view, a_ref)
     int a_view, a_ref;
{
  boolean r = false;
  neighbor_rec nb;
  str2 nb_el;

  if (!is_cyano (a_view, a_ref))
    return false;
  if (atom[a_view - 1].neighbor_count == 1
      && atom[a_view - 1].formal_charge == 0)
    return true;
/* p2c: checkmol.pas: Note: Eliminated unused assignment statement [338] */
  get_nextneighbors (nb, a_view, a_ref);
  strcpy (nb_el, atom[nb[0] - 1].element);
  if (!strcmp (nb_el, "C ")
      || !strcmp (nb_el, "H ") /*|| !strcmp (nb_el, "D ") */ )
    /* v0.3n: D */
    r = true;
  /* HCN is also a nitrile! */
  return r;
}


static boolean
is_isonitrile (a_view, a_ref)
     int a_view, a_ref;
{
  /* only recognized with CN triple bond! */
  boolean r = false;

  if (((strcmp (atom[a_view - 1].atype, "C1 ") == 0) &
       (bond[get_bond (a_view, a_ref) - 1].btype == 'T')) &&
      !strcmp (atom[a_ref - 1].atype, "N1 ") &&
      atom[a_ref - 1].neighbor_count == 2
      && atom[a_view - 1].neighbor_count == 1)
    r = true;
  return r;
}


static boolean
is_cyanate (a_view, a_ref)
     int a_view, a_ref;
{
  boolean r = false;
  neighbor_rec nb;

  if (!is_cyano (a_view, a_ref))
    return false;
  if (atom[a_view - 1].neighbor_count != 2)
    return false;
  get_nextneighbors (nb, a_view, a_ref);
  if (is_alkoxy (a_view, nb[0]) || is_aryloxy (a_view, nb[0]))
    r = true;
  return r;
}


static boolean
is_thiocyanate (a_view, a_ref)
     int a_view, a_ref;
{
  boolean r = false;
  neighbor_rec nb;

  if (!is_cyano (a_view, a_ref))
    return false;
  if (atom[a_view - 1].neighbor_count != 2)
    return false;
  get_nextneighbors (nb, a_view, a_ref);
  if (is_alkylsulfanyl (a_view, nb[0]) || is_arylsulfanyl (a_view, nb[0]))
    r = true;
  return r;
}


static void
update_Htotal ()
{
  int i, j, b_id;
  neighbor_rec nb;
  int single_count, double_count, triple_count, arom_count, total_bonds,
    Htotal, nval;
  /* new in v0.3 */
  boolean diazon = false;	/* new in v0.3j */
  neighbor_rec nb2;		/* new in v0.3j */
  int a1, a2, a3;		/* new in v0.3j */
  int FORLIM, FORLIM1;

  if (n_atoms < 1)
    return;
  memset (nb, 0, sizeof (neighbor_rec));
  FORLIM = n_atoms;
  for (i = 1; i <= FORLIM; i++)
    {
      single_count = 0;
      double_count = 0;
      triple_count = 0;
      arom_count = 0;
      total_bonds = 0;
      Htotal = 0;
      get_neighbors (nb, i);
      if (atom[i - 1].neighbor_count > 0)
	{
	  /* count single, double, triple, and aromatic bonds to all neighbor atoms */
	  FORLIM1 = atom[i - 1].neighbor_count;
	  for (j = 0; j < FORLIM1; j++)
	    {
	      b_id = get_bond (i, nb[j]);
	      if (b_id > 0)
		{
		  if (bond[b_id - 1].btype == 'S')
		    single_count++;
		  if (bond[b_id - 1].btype == 'D')
		    double_count++;
		  if (bond[b_id - 1].btype == 'T')
		    triple_count++;
		  if (bond[b_id - 1].btype == 'A')
		    arom_count++;
		}
	    }
	  /*check for diazonium salts */
	  a1 = i;
	  a2 = nb[0];
	  if (!strcmp (atom[a1 - 1].element, "N ") &&
	      !strcmp (atom[a2 - 1].element, "N "))
	    {
	      if (atom[a2 - 1].neighbor_count == 2)
		{
		  get_nextneighbors (nb2, a2, a1);
		  a3 = nb2[0];
		  if ((strcmp (atom[a3 - 1].element, "C ") ==
		       0) && is_diazonium (a3, a2))
		    diazon = true;
		}
	    }
	}
      total_bonds = single_count + double_count * 2 + triple_count * 3 +
	(int) (1.5 * arom_count);
      /* calculate number of total hydrogens per atom */
      /*nval := nvalences(atom^[i].element);    (* new in v0.3 */
      nval = atom[i - 1].nvalences;	/* new in v0.3m */
      if (!strcmp (atom[i - 1].element, "P "))
	{
	  if (total_bonds - atom[i - 1].formal_charge > 3)	/* refined in v0.3n */
	    nval = 5;
	}			/*  */
      if (!strcmp (atom[i - 1].element, "S "))
	{			/* v0.3h */
	  if (total_bonds > 2 && atom[i - 1].formal_charge < 1)
	    /* updated in v0.3j */
	    nval = 4;
	  if (total_bonds > 4)	/* this will need some refinement... */
	    nval = 6;
	}			/*  */
      Htotal = nval - total_bonds + atom[i - 1].formal_charge;
      if (diazon)		/* v0.3j */
	Htotal = 0;
      if (Htotal < 0)		/* e.g., N in nitro group */
	Htotal = 0;
      atom[i - 1].Htot = Htotal;
      if (atom[i - 1].Hexp > atom[i - 1].Htot)	/* v0.3n; just to be sure... */
	atom[i - 1].Htot = atom[i - 1].Hexp;
    }
}


static void
update_atypes ()
{
  int i, j, b_id;
  neighbor_rec nb;
  int single_count, double_count, triple_count, arom_count, acyl_count,
    C_count, O_count, total_bonds, NdO_count, NdC_count, Htotal, FORLIM,
    FORLIM1;

  if (n_atoms < 1)
    return;
  memset (nb, 0, sizeof (neighbor_rec));
  FORLIM = n_atoms;
  for (i = 0; i < FORLIM; i++)
    {
      single_count = 0;
      double_count = 0;
      triple_count = 0;
      arom_count = 0;
      total_bonds = 0;
      acyl_count = 0;
      C_count = 0;
      O_count = 0;
      NdO_count = 0;
      NdC_count = 0;
      Htotal = 0;
      get_neighbors (nb, i + 1);
      if (atom[i].neighbor_count > 0)
	{
	  /* count single, double, triple, and aromatic bonds to all neighbor atoms */
	  FORLIM1 = atom[i].neighbor_count;
	  for (j = 0; j < FORLIM1; j++)
	    {
	      if (is_oxo_C (nb[j]) || is_thioxo_C (nb[j]))
		acyl_count++;
	      if (!strcmp (atom[nb[j] - 1].element, "C "))
		C_count++;
	      if (!strcmp (atom[nb[j] - 1].element, "O "))
		O_count++;
	      b_id = get_bond (i + 1, nb[j]);
	      if (b_id > 0)
		{
		  if (bond[b_id - 1].btype == 'S')
		    single_count++;
		  if (bond[b_id - 1].btype == 'D')
		    double_count++;
		  if (bond[b_id - 1].btype == 'T')
		    triple_count++;
		  if (bond[b_id - 1].btype == 'A')
		    /* v0.3n: special treatment for acyclic bonds */
		    {		/* flagged as "aromatic" (in query structures) */
		      if (bond[b_id - 1].ring_count > 0)
			arom_count++;
		      else
			double_count++;
		    }
		  if ((!strcmp (atom[i].element, "N ") &&
		       !strcmp (atom[nb[j] - 1].element, "O ")) ||
		      (!strcmp (atom[i].element, "O ") &&
		       !strcmp (atom[nb[j] - 1].element, "N ")))
		    {
		      /* check if it is an N-oxide drawn with a double bond ==> should be N3 */
		      if (bond[b_id - 1].btype == 'D')
			NdO_count++;
		    }
		  if ((!strcmp (atom[i].element, "N ") &&
		       !strcmp (atom[nb[j] - 1].element, "C ")) ||
		      (!strcmp (atom[i].element, "C ") &&
		       !strcmp (atom[nb[j] - 1].element, "N ")))
		    {
		      if (bond[b_id - 1].btype == 'D')
			NdC_count++;
		    }
		}
	    }
	  total_bonds = single_count + double_count * 2 + triple_count * 3 +
	    (int) (1.5 * arom_count);
	  /* calculate number of total hydrogens per atom */
	  /*Htotal := nvalences(atom^[i].element) - total_bonds + atom^[i].formal_charge; */
	  Htotal = atom[i].nvalences - total_bonds + atom[i].formal_charge;
	  if (Htotal < 0)	/* e.g., N in nitro group */
	    Htotal = 0;
	  atom[i].Htot = Htotal;
	  /* refine atom types, based on bond types */
	  if (!strcmp (atom[i].element, "C "))
	    {
	      if (arom_count > 1)
		strcpy (atom[i].atype, "CAR");
	      if (triple_count == 1 || double_count == 2)
		strcpy (atom[i].atype, "C1 ");
	      if (double_count == 1)
		strcpy (atom[i].atype, "C2 ");
	      if (triple_count == 0 && double_count == 0 && arom_count < 2)
		strcpy (atom[i].atype, "C3 ");
	    }
	  if (!strcmp (atom[i].element, "O "))
	    {
	      if (double_count == 1)
		strcpy (atom[i].atype, "O2 ");
	      if (double_count == 0)
		strcpy (atom[i].atype, "O3 ");
	    }
	  if (!strcmp (atom[i].element, "N "))
	    {
	      if (total_bonds > 3)
		{
		  if (O_count == 0)
		    {
		      if (single_count > 3 ||
			  single_count == 2 && double_count == 1
			  && C_count >= 2)
			atom[i].formal_charge = 1;
		    }
		  else
		    {
		      if (O_count == 1 && atom[i].formal_charge == 0)	/* v0.3m */
			strcpy (atom[i].atype, "N3 ");
		      if (O_count == 2 && atom[i].formal_charge == 0)
			{
			  if (atom[i].neighbor_count > 2)	/* nitro v0.3o */
			    strcpy (atom[i].atype, "N2 ");
			  if (atom[i].neighbor_count == 2)	/* NO2   v0.3o */
			    strcpy (atom[i].atype, "N1 ");
			}
		      /* the rest is left empty, so far.... */
		    }
		}
	      /* could be an N-oxide -> should be found elsewhere  */
	      if (triple_count == 1 ||
		  double_count == 2 && atom[i].neighbor_count == 2)
		/* v0.3n */
		strcpy (atom[i].atype, "N1 ");
	      if (double_count == 1)
		{
		  /*if NdC_count > 0 then atom^[i].atype := 'N2 '; */
		  if (NdC_count == 0 && NdO_count > 0 && C_count >= 2)
		    strcpy (atom[i].atype, "N3 ");
		  /* N-oxide is N3 except in hetarene etc. */
		  else
		    strcpy (atom[i].atype, "N2 ");
		}
	      /* fallback, added in v0.3g  */
	      if (arom_count > 1 || atom[i].arom == true)	/* v0.3n */
		strcpy (atom[i].atype, "NAR");
	      if (triple_count == 0 && double_count == 0)
		{
		  if (atom[i].formal_charge == 0)
		    {
		      if (acyl_count == 0)
			strcpy (atom[i].atype, "N3 ");
		      if (acyl_count > 0)
			strcpy (atom[i].atype, "NAM");
		    }
		  if (atom[i].formal_charge == 1)
		    strcpy (atom[i].atype, "N3+");
		}
	    }
	  if (!strcmp (atom[i].element, "P "))
	    {
	      if (single_count > 4)
		strcpy (atom[i].atype, "P4 ");
	      if (single_count <= 4 && double_count == 0)
		strcpy (atom[i].atype, "P3 ");
	      if (double_count == 2)
		strcpy (atom[i].atype, "P3D");
	    }
	  if (!strcmp (atom[i].element, "S "))
	    {
	      if (double_count == 1 && single_count == 0)
		strcpy (atom[i].atype, "S2 ");
	      if (double_count == 0)
		strcpy (atom[i].atype, "S3 ");
	      if (double_count == 1 && single_count > 0)
		strcpy (atom[i].atype, "SO ");
	      if (double_count == 2 && single_count > 0)
		strcpy (atom[i].atype, "SO2");
	    }
	  /* further atom types should go here */
	}
    }
}


static void
chk_arom ()
{
  int i, j, pi_count, ring_size, b, a1, a2;	/* v0.3n */
  ringpath_type testring;
  int a_ref, a_prev, a_next, b_bk, b_fw, b_exo;
  char bt_bk, bt_fw;
  boolean ar_bk, ar_fw, ar_exo;	/* new in v0.3 */
  boolean conj_intr, ko, aromatic, aromatic_bt;	/* v0.3n */
  int n_db, n_sb, n_ar;
  boolean cumul;
  int exo_mC;			/* v0.3j */
  int arom_pi_diff;		/* v0.3j */
  int FORLIM;

  if (n_rings < 1)
    return;
  FORLIM = n_rings;
  /* first, do a very quick check for benzene, pyridine, etc. */
  for (i = 0; i < FORLIM; i++)
    {
      ring_size = ringprop[i].size;
      if (ring_size == 6)
	{
	  memset (testring, 0, sizeof (ringpath_type));
	  for (j = 0; j < ring_size; j++)
	    testring[j] = ring[i][j];
	  cumul = false;
	  n_sb = 0;
	  n_db = 0;
	  n_ar = 0;
	  a_prev = testring[ring_size - 1];
	  for (j = 1; j <= ring_size; j++)
	    {
	      a_ref = testring[j - 1];
	      if (j < ring_size)
		a_next = testring[j];
	      else
		a_next = testring[0];
	      b_bk = get_bond (a_prev, a_ref);
	      b_fw = get_bond (a_ref, a_next);
	      bt_bk = bond[b_bk - 1].btype;
	      bt_fw = bond[b_fw - 1].btype;
	      if (bt_fw == 'S')
		n_sb++;
	      if (bt_fw == 'D')
		n_db++;
	      if (bt_fw == 'A')
		n_ar++;
	      if (bt_fw != 'A' && bt_bk == bt_fw)
		cumul = true;
	      a_prev = a_ref;
	    }
	  if (n_ar == 6 || n_sb == 3 && n_db == 3 && cumul == false)
	    {			/* this ring is aromatic */
	      a_prev = testring[ring_size - 1];
	      for (j = 0; j < ring_size; j++)
		{
		  a_ref = testring[j];
		  b_bk = get_bond (a_prev, a_ref);
		  bond[b_bk - 1].arom = true;
		  a_prev = a_ref;
		}
	      ringprop[i].arom = true;
	    }
	}
    }
  FORLIM = n_rings;
  for (i = 1; i <= FORLIM; i++)
    {
      if (ringprop[i - 1].arom == false)
	{
	  /* do the hard work only for those rings which are not yet flagged aromatic */
	  memset (testring, 0, sizeof (ringpath_type));
	  ring_size = ringprop[i - 1].size;	/* v0.3j */
	  for (j = 0; j < ring_size; j++)	/* v0.3j */
	    testring[j] = ring[i - 1][j];
	  pi_count = 0;
	  arom_pi_diff = 0;	/* v0.3j */
/* p2c: checkmol.pas: Note: Eliminated unused assignment statement [338] */
	  ko = false;
	  a_prev = testring[ring_size - 1];
	  for (j = 1; j <= ring_size; j++)
	    {
	      a_ref = testring[j - 1];
	      if (j < ring_size)
		a_next = testring[j];
	      else
		a_next = testring[0];
	      b_bk = get_bond (a_prev, a_ref);
	      b_fw = get_bond (a_ref, a_next);
	      bt_bk = bond[b_bk - 1].btype;
	      bt_fw = bond[b_fw - 1].btype;
	      ar_bk = bond[b_bk - 1].arom;
	      ar_fw = bond[b_fw - 1].arom;
	      if (bt_bk == 'S' && bt_fw == 'S' && ar_bk == false
		  && ar_fw == false)
		{
		  /* first, assume the worst case (interrupted conjugation) */
		  conj_intr = true;
		  /* conjugation can be restored by hetero atoms */
		  if (!strcmp (atom[a_ref - 1].atype, "O3 ") ||
		      !strcmp (atom[a_ref - 1].atype, "S3 ") ||
		      !strcmp (atom[a_ref - 1].element, "N ") ||
		      !strcmp (atom[a_ref - 1].element, "SE"))
		    {
		      conj_intr = false;
		      pi_count += 2;	/* lone pair adds for 2 pi electrons */
		    }
		  /* conjugation can be restored by a formal charge at a methylene group */
		  if (!strcmp (atom[a_ref - 1].element, "C ") &&
		      atom[a_ref - 1].formal_charge != 0)
		    {
		      conj_intr = false;
		      pi_count -= atom[a_ref - 1].formal_charge;
		      /* neg. charge increases pi_count! */
		    }
		  /* conjugation can be restored by carbonyl groups etc. */
		  if (is_oxo_C (a_ref) || is_thioxo_C (a_ref) |
		      is_exocyclic_imino_C (a_ref, i))
		    conj_intr = false;
		  /* conjugation can be restored by exocyclic C=C double bond, */
		  /* adds 2 pi electrons to 5-membered rings, not to 7-membered rings (CAUTION!) */
		  /* apply only to non-aromatic exocyclic C=C bonds */
		  exo_mC = find_exocyclic_methylene_C (a_ref, i);	/* v0.3j */
		  if (exo_mC > 0 && (ring_size & 1))
		    {		/* v0.3j */
		      b_exo = get_bond (a_ref, exo_mC);	/* v0.3j  */
		      ar_exo = bond[b_exo - 1].arom;
		      if (((ring_size - 1) & 3) == 0)
			{	/* 5-membered rings and related */
			  conj_intr = false;
			  pi_count += 2;
			}
		      else
			{
			  if (!ar_exo)
			    conj_intr = false;
			}
		    }
		  /* 7-membered rings and related */
		  /* if conjugation is still interrupted ==> knock-out */
		  if (conj_intr)
		    ko = true;
		}
	      else
		{
		  if (bt_bk == 'S' && bt_fw == 'S' && ar_bk == true
		      && ar_fw == true)
		    {
		      if (!strcmp (atom[a_ref - 1].atype, "O3 ") ||
			  !strcmp (atom[a_ref - 1].atype, "S3 ") ||
			  !strcmp (atom[a_ref - 1].element, "N ") ||
			  !strcmp (atom[a_ref - 1].element, "SE"))
			pi_count += 2;	/* lone pair adds for 2 pi electrons */
		      if (!strcmp (atom[a_ref - 1].element, "C ") &&
			  atom[a_ref - 1].formal_charge != 0)
			pi_count -= atom[a_ref - 1].formal_charge;
		      /* neg. charge increases pi_count! */
		      exo_mC = find_exocyclic_methylene_C (a_ref, i);	/* v0.3j */
		      if (exo_mC > 0 && (ring_size & 1))
			{	/* v0.3j */
			  b_exo = get_bond (a_ref, exo_mC);	/* v0.3j */
/* p2c: checkmol.pas: Note: Eliminated unused assignment statement [338] */
			  if (((ring_size - 1) & 3) == 0)
			    /* 5-membered rings and related */
			    pi_count += 2;
			}
		    }
		  else
		    {
		      pi_count++;	/* v0.3j; adjustment for bridgehead N: see below */
		      if (bt_bk == 'S' && bt_fw == 'S' &&
			  (ar_bk == true && ar_fw == false ||
			   ar_bk == false && ar_fw == true))
			{
			  /* v0.3j; if a bridgehead N were not aromatic, it could  */
			  /* contribute 2 pi electrons --> try also this variant */
			  /* (example: CAS 32278-54-9) */
			  if (!strcmp (atom[a_ref - 1].element, "N "))
			    {
			      arom_pi_diff++;
			      /* any other case: increase pi count by one electron */
			    }
			}
		    }
		}
	      /* last command: */
	      a_prev = a_ref;
	    }			/* for j := 1 to ring_size */
	  /* now we can draw our conclusion */
	  /*if not ((ko) or (odd(pi_count))) then */
	  if (!ko)
	    /* v0.3j; odd pi_count might be compensated by arom_pi_diff */
	    {			/* apply Hueckel's rule */
	      if (labs (ring_size - pi_count) < 2 &&
		  (((pi_count - 2) & 3) == 0 ||
		   ((pi_count + arom_pi_diff - 2) & 3) == 0))
		{
		  /* this ring is aromatic */
		  ringprop[i - 1].arom = true;
		  /* now mark _all_ bonds in the ring as aromatic */
		  a_prev = testring[ring_size - 1];
		  for (j = 0; j < ring_size; j++)
		    {
		      a_ref = testring[j];
		      bond[get_bond (a_prev, a_ref) - 1].arom = true;
		      a_prev = a_ref;
		    }
		}
	    }
	}
    }				/* (for i := 1 to n_rings) */
  FORLIM = n_bonds;
  /* finally, mark all involved atoms as aromatic */
  for (i = 0; i < FORLIM; i++)
    {
      if (bond[i].arom)
	{
	  a1 = bond[i].a1;	/* v0.3n */
	  a2 = bond[i].a2;	/* v0.3n */
	  atom[a1 - 1].arom = true;
	  atom[a2 - 1].arom = true;
	  /* v0.3n: update atom types if applicable (C and N) */
	  if (!strcmp (atom[a1 - 1].element, "C "))
	    strcpy (atom[a1 - 1].atype, "CAR");
	  if (!strcmp (atom[a2 - 1].element, "C "))
	    strcpy (atom[a2 - 1].atype, "CAR");
	  if (!strcmp (atom[a1 - 1].element, "N "))
	    strcpy (atom[a1 - 1].atype, "NAR");
	  if (!strcmp (atom[a2 - 1].element, "N "))
	    strcpy (atom[a2 - 1].atype, "NAR");
	}
    }
  FORLIM = n_rings;
  /* update aromaticity information in ringprop */
  /* new in v0.3n: accept rings as aromatic if all bonds are of type 'A' */
  for (i = 0; i < FORLIM; i++)
    {
      memcpy (testring, ring[i], sizeof (ringpath_type));
      /*ring_size := path_length(testring); */
      ring_size = ringprop[i].size;	/* v0.3j */
      aromatic = true;
      aromatic_bt = true;	/* v0.3n */
      a_prev = testring[ring_size - 1];
      for (j = 0; j < ring_size; j++)
	{
	  a_ref = testring[j];
	  b = get_bond (a_prev, a_ref);	/* v0.3n */
	  if (!bond[b - 1].arom)
	    aromatic = false;
	  if (bond[b - 1].btype != 'A')	/* v0.3n */
	    aromatic_bt = false;
	  a_prev = a_ref;
	}
      if (aromatic_bt && !aromatic)
	{			/* v0.3n: update aromaticity flag */
	  a_prev = testring[ring_size - 1];
	  for (j = 0; j < ring_size; j++)
	    {
	      a_ref = testring[j];
	      b = get_bond (a_prev, a_ref);
	      bond[b - 1].arom = true;
	      if (!strcmp (atom[a_ref - 1].element, "C "))
		strcpy (atom[a_ref - 1].atype, "CAR");
	      if (!strcmp (atom[a_ref - 1].element, "N "))
		strcpy (atom[a_ref - 1].atype, "NAR");
	      a_prev = a_ref;
	    }
	  aromatic = true;
	}			/* end v0.3n block   */
      if (aromatic)
	ringprop[i].arom = true;
      else
	ringprop[i].arom = false;
    }
}


static void
write_mol ()
{
  int i, j;
  ringpath_type testring;
  int ring_size, FORLIM;

  /*aromatic : boolean; */
  /*a_prev, a_ref : integer; */
  if (progmode == pmCheckMol)
    printf ("Molecule name: %s\n", molname);
  else
    printf ("Molecule name (haystack): %s\n", molname);
  printf ("atoms: %ld  bonds: %ld  rings: %ld\n", n_atoms, n_bonds, n_rings);
  if (n_atoms < 1)
    return;
  if (n_bonds < 1)
    return;
  FORLIM = n_atoms;
  for (i = 1; i <= FORLIM; i++)
    {
      if (i < 10)
	putchar (' ');
      if (i < 100)
	putchar (' ');
      if (i < 1000)
	putchar (' ');
      printf ("%ld %s %s %f %f ",
	      i, atom[i - 1].element, atom[i - 1].atype, atom[i - 1].x,
	      atom[i - 1].y);
      printf ("%f", atom[i - 1].z);
      printf ("  (%ld heavy-atom neighbors, Hexp: %ld Htot: %ld)",
	      atom[i - 1].neighbor_count, atom[i - 1].Hexp, atom[i - 1].Htot);
      if (atom[i - 1].formal_charge != 0)
	printf ("  charge: %ld", atom[i - 1].formal_charge);
      putchar ('\n');
    }
  FORLIM = n_bonds;
  for (i = 1; i <= FORLIM; i++)
    {
      if (i < 10)
	putchar (' ');
      if (i < 100)
	putchar (' ');
      if (i < 1000)
	putchar (' ');
      printf ("%ld %ld %ld %c",
	      i, bond[i - 1].a1, bond[i - 1].a2, bond[i - 1].btype);
      if (bond[i - 1].ring_count > 0)
	printf (", contained in %ld ring(s)", bond[i - 1].ring_count);
      if (bond[i - 1].arom)
	printf (" (aromatic) ");
      putchar ('\n');
    }
  if (n_rings <= 0)
    return;
  FORLIM = n_rings;
  for (i = 0; i < FORLIM; i++)
    {
      printf ("ring %ld: ", i + 1);
      /*aromatic := true; */
      memset (testring, 0, sizeof (ringpath_type));
      ring_size = ringprop[i].size;	/* v0.3j */
      /*for j := 1 to max_ringsize do if ring^[i,j] > 0 then testring[j] := ring^[i,j]; */
      for (j = 0; j < ring_size; j++)	/* v0.3j */
	testring[j] = ring[i][j];
      /*ring_size := path_length(testring); */
      /*a_prev := testring[ring_size]; */
      for (j = 0; j < ring_size; j++)
	{
	  printf ("%ld ", testring[j]);
	  /*a_ref := testring[j]; */
	  /*if (not bond^[get_bond(a_prev,a_ref)].arom) then aromatic := false; */
	  /*a_prev := a_ref; */
	}
      /*if aromatic then write(' (aromatic)'); */
      if (ringprop[i].arom)
	printf (" (aromatic)");
      if (ringprop[i].envelope)
	printf (" (env)");
      putchar ('\n');
    }
}


static void
write_needle_mol ()
{
  int i, j;
  ringpath_type testring;
  int ring_size;
  boolean aromatic;
  int a_prev, a_ref, FORLIM;

  printf ("Molecule name (needle): %s\n", ndl_molname);
  printf ("atoms: %ld  bonds: %ld  rings: %ld\n",
	  ndl_n_atoms, ndl_n_bonds, ndl_n_rings);
  if (ndl_n_atoms < 1)
    return;
  if (ndl_n_bonds < 1)
    return;
  FORLIM = ndl_n_atoms;
  for (i = 1; i <= FORLIM; i++)
    {
      if (i < 10)
	putchar (' ');
      if (i < 100)
	putchar (' ');
      if (i < 1000)
	putchar (' ');
      printf ("%ld %s %s %f %f ",
	      i, ndl_atom[i - 1].element, ndl_atom[i - 1].atype,
	      ndl_atom[i - 1].x, atom[i - 1].y);
      printf ("%f", ndl_atom[i - 1].z);
      printf ("  (%ld heavy-atom neighbors, Hexp: %ld Htot: %ld)",
	      ndl_atom[i - 1].neighbor_count, ndl_atom[i - 1].Hexp,
	      ndl_atom[i - 1].Htot);
      if (ndl_atom[i - 1].formal_charge != 0)
	printf ("  charge: %ld", ndl_atom[i - 1].formal_charge);
      putchar ('\n');
    }
  FORLIM = ndl_n_bonds;
  for (i = 1; i <= FORLIM; i++)
    {
      if (i < 10)
	putchar (' ');
      if (i < 100)
	putchar (' ');
      if (i < 1000)
	putchar (' ');
      printf ("%ld %ld %ld %c",
	      i, ndl_bond[i - 1].a1, ndl_bond[i - 1].a2,
	      ndl_bond[i - 1].btype);
      if (ndl_bond[i - 1].ring_count > 0)
	printf (", contained in %ld ring(s)", ndl_bond[i - 1].ring_count);
      if (ndl_bond[i - 1].arom)
	printf (" (aromatic) ");
      putchar ('\n');
    }
  if (ndl_n_rings <= 0)
    return;
  FORLIM = ndl_n_rings;
  for (i = 0; i < FORLIM; i++)
    {
      aromatic = true;
      memset (testring, 0, sizeof (ringpath_type));
      for (j = 0; j < max_ringsize; j++)
	{
	  if (ndl_ring[i][j] > 0)
	    testring[j] = ndl_ring[i][j];
	}
      ring_size = path_length (testring);
      printf ("ring %ld: ", i + 1);
      a_prev = testring[ring_size - 1];
      for (j = 0; j < ring_size; j++)
	{
	  printf ("%ld ", testring[j]);
	  a_ref = testring[j];
	  if (!ndl_bond[get_ndl_bond (a_prev, a_ref) - 1].arom)	/* v0.3k */
	    aromatic = false;
	  a_prev = a_ref;
	}
      if (aromatic)
	printf (" (aromatic)");
      putchar ('\n');
    }
}


static void
chk_so2_deriv (a_ref)
     int a_ref;
{
  int i;
  neighbor_rec nb;
  str2 nb_el;
  int het_count = 0, o_count = 0, or_count = 0, hal_count = 0, n_count = 0,
    c_count = 0;
  int FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  if (strcmp (atom[a_ref - 1].atype, "SO2"))
    return;
  get_neighbors (nb, a_ref);
  FORLIM = atom[a_ref - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      if (bond[get_bond (a_ref, nb[i]) - 1].btype == 'S')
	{
	  strcpy (nb_el, atom[nb[i] - 1].element);
	  if (strcmp (nb_el, "C ") && strcmp (nb_el, "H ")
	      /*&&  strcmp (nb_el, "D ") */  && strcmp (nb_el, "DU")
	      && strcmp (nb_el, "LP"))
	    /* added 'D ' in v0.3n */
	    het_count++;
	  if (!strcmp (nb_el, "O "))
	    {
	      o_count++;
	      if (is_alkoxy (a_ref, nb[i]) || is_aryloxy (a_ref, nb[i]))
		or_count++;
	    }
	  if (!strcmp (nb_el, "N "))
	    n_count++;
	  if (!strcmp (nb_el, "C "))
	    c_count++;
	  if (!strcmp (nb_el, "F ") || !strcmp (nb_el, "CL") ||
	      !strcmp (nb_el, "BR") || !strcmp (nb_el, "I ")
	      || !strcmp (nb_el, "AT"))
	    hal_count++;
	}
    }
  if (het_count == 2)
    {				/* sulfuric acid derivative */
      fg[fg_sulfuric_acid_deriv - 1] = true;
      if (o_count == 2)
	{
	  if (or_count == 0)
	    fg[fg_sulfuric_acid - 1] = true;
	  if (or_count == 1)
	    fg[fg_sulfuric_acid_monoester - 1] = true;
	  if (or_count == 2)
	    fg[fg_sulfuric_acid_diester - 1] = true;
	}
      if (o_count == 1)
	{
	  if (or_count == 1 && n_count == 1)
	    fg[fg_sulfuric_acid_amide_ester - 1] = true;
	  if (or_count == 0 && n_count == 1)
	    fg[fg_sulfuric_acid_amide - 1] = true;
	}
      if (n_count == 2)
	fg[fg_sulfuric_acid_diamide - 1] = true;
      if (hal_count > 0)
	fg[fg_sulfuryl_halide - 1] = true;
    }
  if (het_count == 1 && c_count == 1)
    {				/* sulfonic acid derivative */
      fg[fg_sulfonic_acid_deriv - 1] = true;
      if (o_count == 1 && or_count == 0)
	fg[fg_sulfonic_acid - 1] = true;
      if (o_count == 1 && or_count == 1)
	fg[fg_sulfonic_acid_ester - 1] = true;
      if (n_count == 1)
	fg[fg_sulfonamide - 1] = true;
      if (hal_count == 1)
	fg[fg_sulfonyl_halide - 1] = true;
    }
  if (het_count == 0 && c_count == 2)	/* sulfone */
    fg[fg_sulfone - 1] = true;
}


static void
chk_p_deriv (a_ref)
     int a_ref;
{
  int i;
  neighbor_rec nb;
  str2 nb_el, dbl_het;
  int het_count;
  int oh_count = 0, or_count = 0, hal_count = 0, n_count = 0, c_count = 0;
  int FORLIM;

  if (strcmp (atom[a_ref - 1].element, "P "))
    return;
  memset (nb, 0, sizeof (neighbor_rec));
  get_neighbors (nb, a_ref);
  *dbl_het = '\0';
/* p2c: checkmol.pas: Note: Eliminated unused assignment statement [338] */
  FORLIM = atom[a_ref - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      if (bond[get_bond (a_ref, nb[i]) - 1].btype == 'D')
	strcpy (dbl_het, atom[nb[i] - 1].element);
      if (bond[get_bond (a_ref, nb[i]) - 1].btype == 'S')
	{
	  strcpy (nb_el, atom[nb[i] - 1].element);
	  if (!strcmp (nb_el, "C "))
	    c_count++;
	  if (is_hydroxy (a_ref, nb[i]))
	    oh_count++;
	  if (is_alkoxy (a_ref, nb[i]) || is_aryloxy (a_ref, nb[i]))
	    or_count++;
	  if (!strcmp (nb_el, "N "))
	    n_count++;
	  if (!strcmp (nb_el, "F ") || !strcmp (nb_el, "CL") ||
	      !strcmp (nb_el, "BR") || !strcmp (nb_el, "I ")
	      || !strcmp (nb_el, "AT"))
	    hal_count++;
	}
    }
  het_count = oh_count + or_count + hal_count + n_count;
  if (!strcmp (atom[a_ref - 1].atype, "P3D") ||
      !strcmp (atom[a_ref - 1].atype, "P4 "))
    {
      if (!strcmp (dbl_het, "O "))
	{
	  if (c_count == 0)
	    {
	      fg[fg_phosphoric_acid_deriv - 1] = true;
	      if (oh_count == 3)
		fg[fg_phosphoric_acid - 1] = true;
	      if (or_count > 0)
		fg[fg_phosphoric_acid_ester - 1] = true;
	      if (hal_count > 0)
		fg[fg_phosphoric_acid_halide - 1] = true;
	      if (n_count > 0)
		fg[fg_phosphoric_acid_amide - 1] = true;
	    }
	  if (c_count == 1)
	    {
	      fg[fg_phosphonic_acid_deriv - 1] = true;
	      if (oh_count == 2)
		fg[fg_phosphonic_acid - 1] = true;
	      if (or_count > 0)
		fg[fg_phosphonic_acid_ester - 1] = true;
	      /*if (hal_count > 0)  then fg[fg_phosphonic_acid_halide] := true;             */
	      /*if (n_count > 0)    then fg[fg_phosphonic_acid_amide]  := true; */
	    }
	  if (c_count == 3)
	    fg[fg_phosphinoxide - 1] = true;
	}
      if (!strcmp (dbl_het, "S "))
	{
	  if (c_count == 0)
	    {
	      fg[fg_thiophosphoric_acid_deriv - 1] = true;
	      if (oh_count == 3)
		fg[fg_thiophosphoric_acid - 1] = true;
	      if (or_count > 0)
		fg[fg_thiophosphoric_acid_ester - 1] = true;
	      if (hal_count > 0)
		fg[fg_thiophosphoric_acid_halide - 1] = true;
	      if (n_count > 0)
		fg[fg_thiophosphoric_acid_amide - 1] = true;
	    }
	}
    }
  /*  if (atom^[a_ref].atype = 'P4 ') then fg[fg_phosphoric_acid_deriv] := true; */
  if (strcmp (atom[a_ref - 1].atype, "P3 "))	/* changed P3D into P3 in v0.3b */
    return;
  if (c_count == 3 && het_count == 0)
    fg[fg_phosphine - 1] = true;
  if (c_count == 3 && oh_count == 1)
    fg[fg_phosphinoxide - 1] = true;
}


static void
chk_b_deriv (a_ref)
     int a_ref;
{
  int i;
  neighbor_rec nb;
  str2 nb_el;
  int het_count = 0, oh_count = 0, or_count = 0, hal_count = 0, n_count = 0,
    c_count = 0;
  int FORLIM;

  if (strcmp (atom[a_ref - 1].element, "B "))
    return;
  memset (nb, 0, sizeof (neighbor_rec));
  get_neighbors (nb, a_ref);
  FORLIM = atom[a_ref - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      if (bond[get_bond (a_ref, nb[i]) - 1].btype == 'S')
	{
	  strcpy (nb_el, atom[nb[i] - 1].element);
	  if (!strcmp (nb_el, "C "))
	    c_count++;
	  else if (strcmp (nb_el, "H ") /*&& strcmp (nb_el, "D ") */  &&
		   strcmp (nb_el, "LP"))
	    /* v0.3n: D */
	    het_count++;
	  if (is_hydroxy (a_ref, nb[i]))
	    oh_count++;
	  if (is_alkoxy (a_ref, nb[i]) || is_aryloxy (a_ref, nb[i]))
	    /* fixed in v0.3b */
	    or_count++;
	  if (!strcmp (nb_el, "N "))
	    n_count++;
	  if (!strcmp (nb_el, "F ") || !strcmp (nb_el, "CL") ||
	      !strcmp (nb_el, "BR") || !strcmp (nb_el, "I ")
	      || !strcmp (nb_el, "AT"))
	    hal_count++;
	}
    }
  het_count = oh_count + or_count + hal_count + n_count;
  /* fixed in v0.3b */
  if (c_count != 1 || het_count != 2)
    return;
  fg[fg_boronic_acid_deriv - 1] = true;
  if (oh_count == 2)
    fg[fg_boronic_acid - 1] = true;
  if (or_count > 0)
    fg[fg_boronic_acid_ester - 1] = true;
}


static void
chk_ammon (a_ref)
     int a_ref;
{
  int i;
  neighbor_rec nb;
  str2 nb_el;
  int het_count = 0, o_count = 0, or_count = 0, r_count = 0;
  char bt;			/* v0.3k */
  float bo_sum = 0.0;
  boolean ha;
  int FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  if (strcmp (atom[a_ref - 1].atype, "N3+")
      && atom[a_ref - 1].formal_charge == 0)
    return;
  if (strcmp (atom[a_ref - 1].element, "N "))	/* just to be sure;  v0.3i */
    return;
  get_neighbors (nb, a_ref);
  FORLIM = atom[a_ref - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      bt = bond[get_bond (a_ref, nb[i]) - 1].btype;	/* v0.3k */
      strcpy (nb_el, atom[nb[i] - 1].element);	/* v0.3k */
      ha = atom[nb[i] - 1].heavy;	/* v0.3k */
      if (bt == 'S')
	{
	  if (ha)
	    bo_sum += 1.0;
	  if (strcmp (nb_el, "C ") && strcmp (nb_el, "H ")
	      /*&&  strcmp (nb_el, "D ") */  && strcmp (nb_el, "DU"))
	    {			/* added 'D ' in v0.3n */
	      het_count++;
	      if (!strcmp (nb_el, "O "))
		{
		  o_count++;
		  if (atom[nb[i] - 1].neighbor_count > 1)
		    or_count++;
		}
	    }
	  if (is_alkyl (a_ref, nb[i]) || is_aryl (a_ref, nb[i]) |
	      is_alkenyl (a_ref, nb[i]) || is_alkynyl (a_ref, nb[i]))
	    /* v0.3k */
	    r_count++;
	}
      if (bt == 'D')
	{
	  if (ha)
	    bo_sum += 2.0;
	  if (strcmp (nb_el, "C "))
	    {
	      het_count += 2;
	      if (!strcmp (nb_el, "O "))
		o_count += 2;
	    }
	  if (!strcmp (nb_el, "C "))
	    r_count++;
	}
      if (bt == 'A' && ha)
	bo_sum += 1.5;
    }				/* v0.3k: corrected end of "for ..." loop */
  if (het_count == 0 && r_count == 4)
    fg[fg_quart_ammonium - 1] = true;
  if (het_count != 1 || atom[a_ref - 1].neighbor_count < 3)
    return;
  if (o_count == 1 && or_count == 0 && bo_sum > 3)
    fg[fg_n_oxide - 1] = true;	/* finds only aliphatic N-oxides! */
  if ((o_count == 1 && or_count == 1 || o_count == 0) &&
      atom[a_ref - 1].arom == true)
    fg[fg_quart_ammonium - 1] = true;
}


static void
swap_atoms (a1, a2)
     int *a1, *a2;
{
  int a_tmp;

  a_tmp = *a1;
  *a1 = *a2;
  *a2 = a_tmp;
}


static void
orient_bond (a1, a2)
     int *a1, *a2;
{
  str2 a1_el, a2_el;

  strcpy (a1_el, atom[*a1 - 1].element);
  strcpy (a2_el, atom[*a2 - 1].element);
  if (!strcmp (a1_el, "H ") || !strcmp (a2_el, "H ")
      || !strcmp (a1_el, "D ") || !strcmp (a2_el, "D "))
    /* v0.3n: D */
    return;
  if (!strcmp (a2_el, "C ") && strcmp (a1_el, "C "))
    swap_atoms (a1, a2);
  if (!strcmp (a2_el, a1_el))
    {
      if (hetbond_count (*a1) > hetbond_count (*a2))
	swap_atoms (a1, a2);
    }
  if (strcmp (a2_el, "C ") && strcmp (a1_el, "C ") && strcmp (a1_el, a2_el))
    {
      if (!strcmp (a1_el, "O ") || !strcmp (a2_el, "O "))
	{
	  if (!strcmp (a1_el, "O "))
	    swap_atoms (a1, a2);
	}
    }
  if (strcmp (a2_el, "C ") && strcmp (a1_el, "C ") && !strcmp (a1_el, a2_el))
    {
      if (atom[*a2 - 1].neighbor_count - hetbond_count (*a2) >
	  atom[*a1 - 1].neighbor_count - hetbond_count (*a1))
	swap_atoms (a1, a2);
    }
}


static void
chk_imine (a_ref, a_view)
     int a_ref, a_view;
{
  /* a_ref = C, a_view = N */
  int i;
  neighbor_rec nb;
  str2 nb_el;
  int a_het, a_c;
  int het_count = 0, c_count = 0, o_count = 0;	/* v0.3k */
  int FORLIM;

  /* v0.3k */
  if (atom[a_view - 1].neighbor_count == 1)
    {
      if (atom[a_ref - 1].arom == false)
	fg[fg_imine - 1] = true;
      return;
    }
  memset (nb, 0, sizeof (neighbor_rec));
  get_neighbors (nb, a_view);
  if (atom[a_view - 1].neighbor_count <= 1)
    return;
  FORLIM = atom[a_view - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      if ((nb[i] != a_ref) && (bond[get_bond (a_view, nb[i]) - 1].btype ==
			       'S'))
	{
	  strcpy (nb_el, atom[nb[i] - 1].element);
	  if (!strcmp (nb_el, "C "))
	    {
	      a_c = nb[i];
	      c_count++;
	    }
	  if (!strcmp (nb_el, "O ") || !strcmp (nb_el, "N "))
	    {
	      a_het = nb[i];
	      het_count++;
	    }
	  if ((!strcmp (nb_el, "O ")
	       && atom[nb[i] - 1].neighbor_count ==
	       1) && (bond[get_bond (a_view, nb[i]) - 1].arom == false))
	    /* v0.3k */
	    o_count++;
	}
      if ((nb[i] != a_ref) && (bond[get_bond (a_view, nb[i]) - 1].btype ==
			       'D'))
	{			/* v0.3k; make sure we do not count nitro groups in "azi" form etc. */
	  strcpy (nb_el, atom[nb[i] - 1].element);
	  if (!strcmp (nb_el, "O ") || !strcmp (nb_el, "N ")
	      || !strcmp (nb_el, "S "))
	    {
	      a_het = nb[i];	/* v0.3m */
	      het_count++;
	    }
	  if ((!strcmp (nb_el, "O ")
	       && atom[nb[i] - 1].neighbor_count ==
	       1) && (bond[get_bond (a_view, nb[i]) - 1].arom == false))
	    /* v0.3k */
	    o_count++;
	}
    }
  if (c_count == 1)
    {
      if ((is_alkyl (a_view, a_c) || is_aryl (a_view, a_c) |
	   is_alkenyl (a_view, a_c) || is_alkynyl (a_view, a_c))
	  && atom[a_ref - 1].arom == false && het_count == 0)
	/* v0.3k */
	fg[fg_imine - 1] = true;
    }
  if (het_count == 1)
    {
      strcpy (nb_el, atom[a_het - 1].element);
      if (!strcmp (nb_el, "O "))
	{
	  if (is_hydroxy (a_view, a_het))
	    fg[fg_oxime - 1] = true;
	  if (is_alkoxy (a_view, a_het) || is_aryloxy (a_view, a_het) |
	      is_alkenyloxy (a_view, a_het) || is_alkynyloxy (a_view, a_het))
	    fg[fg_oxime_ether - 1] = true;
	}
      if (!strcmp (nb_el, "N "))
	{
	  if (is_amino (a_view, a_het) || is_alkylamino (a_view, a_het) |
	      is_dialkylamino (a_view, a_het) || is_alkylarylamino (a_view,
								    a_het) |
	      is_arylamino (a_view, a_het) || is_diarylamino (a_view, a_het))
	    fg[fg_hydrazone - 1] = true;
	  else
	    {
	      memset (nb, 0, sizeof (neighbor_rec));
	      get_neighbors (nb, a_het);
	      if (atom[a_het - 1].neighbor_count > 1)
		{
		  FORLIM = atom[a_het - 1].neighbor_count;
		  for (i = 0; i < FORLIM; i++)
		    {
		      if (nb[i] != a_view)
			{
			  if (is_carbamoyl (a_het, nb[i]))
			    fg[fg_semicarbazone - 1] = true;
			  if (is_thiocarbamoyl (a_het, nb[i]))
			    fg[fg_thiosemicarbazone - 1] = true;
			}
		    }
		}
	    }
	}
    }				/* v0.3k: nitro groups in "azi" form */
  /* check for semicarbazone or thiosemicarbazone */
  if (het_count == 2 && o_count == 2)
    fg[fg_nitro_compound - 1] = true;
}


static void
chk_carbonyl_deriv (a_view, a_ref)
     int a_view, a_ref;
{
  /* a_view = C */
  int i;
  neighbor_rec nb;
  str2 nb_el;
  int c_count = 0, cn_count = 0;
  char bt;			/* new in v0.3b */
  int n_db = 0;			/* new in v0.3b */
  int FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  get_neighbors (nb, a_view);
  FORLIM = atom[a_view - 1].neighbor_count;
  /* new in v0.3b */
  for (i = 0; i < FORLIM; i++)
    {
      bt = bond[get_bond (a_view, nb[i]) - 1].btype;
      if (bt == 'S')
	{
	  strcpy (nb_el, atom[nb[i] - 1].element);
	  if (!strcmp (nb_el, "C "))
	    {
	      if (is_cyano_c (nb[i]))
		cn_count++;
	      else
		c_count++;
	    }
	}
      else
	{
	  if (bt == 'D')
	    n_db++;
	}
    }
  /* new in v0.3b */
  if (is_oxo_C (a_view))
    {
      fg[fg_carbonyl - 1] = true;
      if (c_count + cn_count < 2)
	{			/* new in v0.3b (detection of ketenes) */
	  if (n_db <= 1)
	    fg[fg_aldehyde - 1] = true;
	  else
	    fg[fg_ketene - 1] = true;
	}
      if (c_count == 2)
	{
	  if (atom[a_view - 1].arom)
	    fg[fg_oxohetarene - 1] = true;
	  else
	    fg[fg_ketone - 1] = true;
	}
      if (cn_count > 0)
	fg[fg_acyl_cyanide - 1] = true;
    }
  if (is_thioxo_C (a_view))
    {
      fg[fg_thiocarbonyl - 1] = true;
      if (c_count < 2)
	fg[fg_thioaldehyde - 1] = true;
      if (c_count == 2)
	{
	  if (atom[a_view - 1].arom)
	    fg[fg_thioxohetarene - 1] = true;
	  else
	    fg[fg_thioketone - 1] = true;
	}
    }
  if (is_imino_C (a_view))
    chk_imine (a_view, a_ref);
}


static void
chk_carboxyl_deriv (a_view, a_ref)
     int a_view, a_ref;
{
  int i;
  neighbor_rec nb;
  str2 nb_el;
  int o_count = 0, n_count = 0, s_count = 0;
  int a_o, a_n, a_s, FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  get_neighbors (nb, a_view);
  FORLIM = atom[a_view - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      if (bond[get_bond (a_view, nb[i]) - 1].btype == 'S')
	{
	  strcpy (nb_el, atom[nb[i] - 1].element);
	  if (strcmp (nb_el, "C "))
	    {
	      if (!strcmp (nb_el, "O "))
		{
		  o_count++;
		  a_o = nb[i];
		}
	      if (!strcmp (nb_el, "N "))
		{
		  n_count++;
		  a_n = nb[i];
		}
	      if (!strcmp (nb_el, "S "))
		{
		  s_count++;
		  a_s = nb[i];
		}
	    }
	}
    }
  if (is_oxo_C (a_view))
    {
      if (o_count == 1)
	{			/* anhydride is checked somewhere else */
	  if (bond[get_bond (a_view, a_o) - 1].arom == false)
	    fg[fg_carboxylic_acid_deriv - 1] = true;
	  if (is_hydroxy (a_view, a_o))
	    {
	      if (atom[a_o - 1].formal_charge == 0)
		fg[fg_carboxylic_acid - 1] = true;
	      if (atom[a_o - 1].formal_charge == -1)
		fg[fg_carboxylic_acid_salt - 1] = true;
	    }
	  if (is_alkoxy (a_view, a_o) || is_aryloxy (a_view, a_o) |
	      is_alkenyloxy (a_view, a_o) || is_alkynyloxy (a_view, a_o))
	    {
	      if (bond[get_bond (a_view, a_o) - 1].arom == false)
		fg[fg_carboxylic_acid_ester - 1] = true;
	      if (bond[get_bond (a_view, a_o) - 1].ring_count > 0)
		{
		  if (bond[get_bond (a_view, a_o) - 1].arom == true)
		    {
		      /*fg[fg_lactone_heteroarom] := true else fg[fg_lactone] := true; */
		      fg[fg_oxohetarene - 1] = true;
		    }
		  else
		    fg[fg_lactone - 1] = true;
		}
	    }
	}
      if (n_count == 1)
	{
	  if (bond[get_bond (a_view, a_n) - 1].arom == false)
	    fg[fg_carboxylic_acid_deriv - 1] = true;
	  else
	    {
	      /*fg[fg_lactam_heteroarom] := true;  (* catches also pyridazines, 1,2,3-triazines, etc. */
	      fg[fg_oxohetarene - 1] = true;
	    }
	  if (is_amino (a_view, a_n)
	      || (!strcmp (atom[a_n - 1].atype, "NAM")
		  && atom[a_n - 1].neighbor_count == 1))
	    {
	      fg[fg_carboxylic_acid_amide - 1] = true;
	      fg[fg_carboxylic_acid_prim_amide - 1] = true;
	    }
	  /*if (is_alkylamino(a_view,a_n)) or (is_arylamino(a_view,a_n)) then  */
	  if (is_C_monosubst_amino (a_view, a_n) &
	      (!is_subst_acylamino (a_view, a_n)))
	    {			/* v0.3j */
	      if (bond[get_bond (a_view, a_n) - 1].arom == false)
		fg[fg_carboxylic_acid_amide - 1] = true;
	      if (bond[get_bond (a_view, a_n) - 1].arom == false)
		fg[fg_carboxylic_acid_sec_amide - 1] = true;
	      if (bond[get_bond (a_view, a_n) - 1].ring_count > 0)
		{
		  if (bond[get_bond (a_view, a_n) - 1].arom == true)
		    {
		      /*fg[fg_lactam_heteroarom]    := true else  */
		      fg[fg_oxohetarene - 1] = true;
		    }
		  else
		    fg[fg_lactam - 1] = true;
		}
	    }
	  /*if (is_dialkylamino(a_view,a_n)) or (is_alkylarylamino(a_view,a_n)) or */
	  /*   (is_diarylamino(a_view,a_n)) then  */
	  if (is_C_disubst_amino (a_view, a_n) &
	      (!is_subst_acylamino (a_view, a_n)))
	    {			/* v0.3j */
	      if (bond[get_bond (a_view, a_n) - 1].arom == false)
		fg[fg_carboxylic_acid_amide - 1] = true;
	      if (bond[get_bond (a_view, a_n) - 1].arom == false)
		fg[fg_carboxylic_acid_tert_amide - 1] = true;
	      if (bond[get_bond (a_view, a_n) - 1].ring_count > 0)
		{
		  if (bond[get_bond (a_view, a_n) - 1].arom == true)
		    {
		      /*fg[fg_lactam_heteroarom]    := true else  */
		      fg[fg_oxohetarene - 1] = true;
		    }
		  else
		    fg[fg_lactam - 1] = true;
		}
	    }
	  if (is_hydroxylamino (a_view, a_n))
	    fg[fg_hydroxamic_acid - 1] = true;
	  if (is_hydrazino (a_view, a_n))
	    fg[fg_carboxylic_acid_hydrazide - 1] = true;
	  if (is_azido (a_view, a_n))
	    fg[fg_carboxylic_acid_azide - 1] = true;
	}
      if (s_count == 1)
	{			/* anhydride is checked somewhere else */
	  if (bond[get_bond (a_view, a_s) - 1].arom == false)
	    fg[fg_thiocarboxylic_acid_deriv - 1] = true;
	  if (is_sulfanyl (a_view, a_s))
	    fg[fg_thiocarboxylic_acid - 1] = true;
	  if (is_alkylsulfanyl (a_view, a_s) || is_arylsulfanyl (a_view, a_s))
	    {
	      if (bond[get_bond (a_view, a_s) - 1].arom == false)
		fg[fg_thiocarboxylic_acid_ester - 1] = true;
	      if (bond[get_bond (a_view, a_s) - 1].ring_count > 0)
		{
		  if (bond[get_bond (a_view, a_s) - 1].arom == true)
		    {
		      /*fg[fg_thiolactone_heteroarom] := true else fg[fg_thiolactone] := true; */
		      fg[fg_oxohetarene - 1] = true;
		    }
		  else
		    fg[fg_thiolactone - 1] = true;
		}
	    }
	}
    }				/* end Oxo-C */
  if (is_thioxo_C (a_view))
    {
      /* fg[fg_thiocarboxylic_acid_deriv]  := true; */
      if (o_count == 1)
	{			/* anhydride is checked somewhere else */
	  if (bond[get_bond (a_view, a_o) - 1].arom == false)
	    fg[fg_thiocarboxylic_acid_deriv - 1] = true;
	  if (is_hydroxy (a_view, a_o))
	    fg[fg_thiocarboxylic_acid - 1] = true;	/* fixed in v0.3c */
	  if (is_alkoxy (a_view, a_o) || is_aryloxy (a_view, a_o))
	    {
	      if (bond[get_bond (a_view, a_s) - 1].arom == false)
		fg[fg_thiocarboxylic_acid_ester - 1] = true;
	      if (bond[get_bond (a_view, a_o) - 1].ring_count > 0)
		{
		  if (bond[get_bond (a_view, a_o) - 1].arom == true)
		    {
		      /*fg[fg_thiolactone_heteroarom] := true else fg[fg_thiolactone] := true; */
		      fg[fg_thioxohetarene - 1] = true;
		    }
		  else
		    fg[fg_thiolactone - 1] = true;
		}
	    }
	}
      if (n_count == 1)
	{
	  if (bond[get_bond (a_view, a_n) - 1].arom == false)
	    fg[fg_thiocarboxylic_acid_deriv - 1] = true;
	  else
	    {
	      /*fg[fg_thiolactam_heteroarom] := true;  (* catches also pyridazines, 1,2,3-triazines, etc. */
	      fg[fg_thioxohetarene - 1] = true;
	    }
	  /* catches also pyridazines, 1,2,3-triazines, etc. */
	  if (is_amino (a_view, a_n))
	    {
	      fg[fg_thiocarboxylic_acid_amide - 1] = true;
	      /* fg[fg_thiocarboxylic_acid_prim_amide] := true; */
	    }
	  /*if (is_alkylamino(a_view,a_n)) or (is_arylamino(a_view,a_n)) then  */
	  if (is_C_monosubst_amino (a_view, a_n) &
	      (!is_subst_acylamino (a_view, a_n)))
	    {			/* v0.3j */
	      if (bond[get_bond (a_view, a_n) - 1].arom == false)
		fg[fg_thiocarboxylic_acid_amide - 1] = true;
	      /*fg[fg_thiocarboxylic_acid_sec_amide]  := true; */
	      if (bond[get_bond (a_view, a_n) - 1].ring_count > 0)
		{
		  if (bond[get_bond (a_view, a_n) - 1].arom == true)
		    {
		      /*fg[fg_thiolactam_heteroarom] := true else fg[fg_thiolactam] := true; */
		      fg[fg_thioxohetarene - 1] = true;
		    }
		  else
		    fg[fg_thiolactam - 1] = true;
		}
	    }
	  /*if (is_dialkylamino(a_view,a_n)) or (is_alkylarylamino(a_view,a_n)) or */
	  /*   (is_diarylamino(a_view,a_n)) then  */
	  if (is_C_disubst_amino (a_view, a_n) &
	      (!is_subst_acylamino (a_view, a_n)))
	    {			/* v0.3j */
	      if (bond[get_bond (a_view, a_n) - 1].arom == false)
		fg[fg_thiocarboxylic_acid_amide - 1] = true;
	      /*fg[fg_thiocarboxylic_acid_tert_amide] := true; */
	      if (bond[get_bond (a_view, a_n) - 1].ring_count > 0)
		{
		  if (bond[get_bond (a_view, a_n) - 1].arom == true)
		    {
		      /*fg[fg_thiolactam_heteroarom] := true else fg[fg_thiolactam] := true; */
		      fg[fg_thioxohetarene - 1] = true;
		    }
		  else
		    fg[fg_thiolactam - 1] = true;
		}
	    }
	}
      if (s_count == 1)
	{			/* anhydride is checked somewhere else */
	  if (bond[get_bond (a_view, a_s) - 1].arom == false)
	    fg[fg_thiocarboxylic_acid_deriv - 1] = true;
	  if (is_sulfanyl (a_view, a_s))
	    fg[fg_thiocarboxylic_acid - 1] = true;
	  if (is_alkylsulfanyl (a_view, a_s) || is_arylsulfanyl (a_view, a_s))
	    {
	      if (bond[get_bond (a_view, a_s) - 1].arom == false)
		fg[fg_thiocarboxylic_acid_ester - 1] = true;
	      if (bond[get_bond (a_view, a_s) - 1].ring_count > 0)
		{
		  if (bond[get_bond (a_view, a_s) - 1].arom == true)
		    {
		      /*fg[fg_thiolactone_heteroarom] := true else fg[fg_thiolactone] := true; */
		      fg[fg_thioxohetarene - 1] = true;
		    }
		  else
		    fg[fg_thiolactone - 1] = true;
		}
	    }
	}
    }				/* end Thioxo-C */
  if (is_true_imino_C (a_view))
    {
      if (o_count == 1)
	{
	  if (bond[get_bond (a_view, a_o) - 1].arom == false)
	    fg[fg_carboxylic_acid_deriv - 1] = true;
	  if (is_alkoxy (a_view, a_o) || is_aryloxy (a_view, a_o))
	    {
	      if (bond[get_bond (a_view, a_o) - 1].arom == false)
		fg[fg_imido_ester - 1] = true;
	    }
	}
      if ((n_count == 1) && (bond[get_bond (a_view, a_n) - 1].arom == false))
	{
	  if (bond[get_bond (a_view, a_n) - 1].arom == false)
	    fg[fg_carboxylic_acid_deriv - 1] = true;
	  if (is_amino (a_view, a_n) || is_subst_amino (a_view, a_n))
	    {
	      if (bond[get_bond (a_view, a_n) - 1].arom == false)
		fg[fg_carboxylic_acid_deriv - 1] = true;
	      fg[fg_carboxylic_acid_amidine - 1] = true;
	    }
	  if (is_hydrazino (a_view, a_n))
	    {
	      if (bond[get_bond (a_view, a_n) - 1].arom == false)
		fg[fg_carboxylic_acid_amidrazone - 1] = true;
	    }
	}
      if ((n_count == 1) && (bond[get_bond (a_view, a_n) - 1].arom == true))
	/* catches also pyridazines, 1,2,3-triazines, etc. */
	fg[fg_iminohetarene - 1] = true;
      if (s_count == 1)
	{
	  if (bond[get_bond (a_view, a_s) - 1].arom == false)
	    fg[fg_carboxylic_acid_deriv - 1] = true;
	  if (is_alkylsulfanyl (a_view, a_s) || is_arylsulfanyl (a_view, a_s))
	    {
	      if (bond[get_bond (a_view, a_s) - 1].arom == false)
		fg[fg_imido_thioester - 1] = true;
	    }
	}
    }
  if (is_hydroximino_C (a_view))
    {
      if (bond[get_bond (a_view, a_n) - 1].arom == false)
	fg[fg_carboxylic_acid_deriv - 1] = true;
      if (o_count == 1)
	{
	  if (is_hydroxy (a_view, a_o))
	    fg[fg_hydroxamic_acid - 1] = true;
	}
    }
  if (!is_hydrazono_C (a_view))
    return;
  if (bond[get_bond (a_view, a_n) - 1].arom == false)
    fg[fg_carboxylic_acid_deriv - 1] = true;
  if (n_count == 1)
    {
      if (is_amino (a_view, a_n) || is_subst_amino (a_view, a_n))
	fg[fg_carboxylic_acid_amidrazone - 1] = true;
    }
}


static void
chk_co2_sp2 (a_view, a_ref)
     int a_view, a_ref;
{
  int i;
  neighbor_rec nb;
  str2 nb_el;
  int o_count = 0, or_count = 0, n_count = 0, nn_count = 0, nnx_count = 0,
    s_count = 0, sr_count = 0;
  int FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  get_neighbors (nb, a_view);
  FORLIM = atom[a_view - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      if (bond[get_bond (a_view, nb[i]) - 1].btype == 'S')
	{
	  strcpy (nb_el, atom[nb[i] - 1].element);
	  if (strcmp (nb_el, "C "))
	    {
	      if (!strcmp (nb_el, "O "))
		{
		  o_count++;
		  if (is_alkoxy (a_view, nb[i]) |
		      is_alkenyloxy (a_view, nb[i]) || is_aryloxy (a_view,
								   nb[i]))
		    /* v0.3j */
		    or_count++;
		}
	      if (!strcmp (nb_el, "N "))
		{
		  n_count++;
		  if (is_hydrazino (a_view, nb[i]))
		    nn_count++;
		  if (is_subst_hydrazino (a_view, nb[i]))	/* more general... */
		    nnx_count++;
		}
	      if (!strcmp (nb_el, "S "))
		{
		  s_count++;
		  if (is_alkylsulfanyl (a_view, nb[i]) |
		      is_arylsulfanyl (a_view, nb[i]))
		    sr_count++;
		}
	    }
	}
    }
  if (is_oxo_C (a_view))
    {
      if (o_count == 2)
	{
	  fg[fg_carbonic_acid_deriv - 1] = true;
	  if (or_count == 1)
	    fg[fg_carbonic_acid_monoester - 1] = true;
	  if (or_count == 2)
	    fg[fg_carbonic_acid_diester - 1] = true;
	}
      if (o_count == 1 && s_count == 1)
	{
	  fg[fg_thiocarbonic_acid_deriv - 1] = true;
	  if (or_count + sr_count == 1)
	    fg[fg_thiocarbonic_acid_monoester - 1] = true;
	  if (or_count + sr_count == 2)
	    fg[fg_thiocarbonic_acid_diester - 1] = true;
	}
      if (s_count == 2)
	{
	  fg[fg_thiocarbonic_acid_deriv - 1] = true;
	  if (sr_count == 1)
	    fg[fg_thiocarbonic_acid_monoester - 1] = true;
	  if (sr_count == 2)
	    fg[fg_thiocarbonic_acid_diester - 1] = true;
	}
      if (o_count == 1 && n_count == 1)
	{
	  fg[fg_carbamic_acid_deriv - 1] = true;
	  if (or_count == 0)
	    fg[fg_carbamic_acid - 1] = true;
	  if (or_count == 1)
	    fg[fg_carbamic_acid_ester - 1] = true;
	}
      if (s_count == 1 && n_count == 1)
	{
	  fg[fg_thiocarbamic_acid_deriv - 1] = true;
	  if (sr_count == 0)
	    fg[fg_thiocarbamic_acid - 1] = true;
	  if (sr_count == 1)
	    fg[fg_thiocarbamic_acid_ester - 1] = true;
	}
      if (n_count == 2)
	{
	  if (nn_count == 1)
	    fg[fg_semicarbazide - 1] = true;
	  else
	    {
	      if (nnx_count == 0)	/* excludes semicarbazones */
		fg[fg_urea - 1] = true;
	    }
	}
    }				/* end Oxo-C */
  if (is_thioxo_C (a_view))
    {
      if (o_count == 2)
	{
	  fg[fg_thiocarbonic_acid_deriv - 1] = true;
	  if (or_count == 1)
	    fg[fg_thiocarbonic_acid_monoester - 1] = true;
	  if (or_count == 2)
	    fg[fg_thiocarbonic_acid_diester - 1] = true;
	}
      if (o_count == 1 && s_count == 1)
	{
	  fg[fg_thiocarbonic_acid_deriv - 1] = true;
	  if (or_count + sr_count == 1)
	    fg[fg_thiocarbonic_acid_monoester - 1] = true;
	  if (or_count + sr_count == 2)
	    fg[fg_thiocarbonic_acid_diester - 1] = true;
	}
      if (s_count == 2)
	{
	  fg[fg_thiocarbonic_acid_deriv - 1] = true;
	  if (sr_count == 1)
	    fg[fg_thiocarbonic_acid_monoester - 1] = true;
	  if (sr_count == 2)
	    fg[fg_thiocarbonic_acid_diester - 1] = true;
	}
      if (o_count == 1 && n_count == 1)
	{
	  fg[fg_thiocarbamic_acid_deriv - 1] = true;
	  if (or_count == 0)
	    fg[fg_thiocarbamic_acid - 1] = true;
	  if (or_count == 1)
	    fg[fg_thiocarbamic_acid_ester - 1] = true;
	}
      if (s_count == 1 && n_count == 1)
	{
	  fg[fg_thiocarbamic_acid_deriv - 1] = true;
	  if (sr_count == 0)
	    fg[fg_thiocarbamic_acid - 1] = true;
	  if (sr_count == 1)
	    fg[fg_thiocarbamic_acid_ester - 1] = true;
	}
      if (n_count == 2)
	{
	  if (nn_count == 1)
	    fg[fg_thiosemicarbazide - 1] = true;
	  else
	    {
	      if (nnx_count == 0)	/* excludes thiosemicarbazones */
		fg[fg_thiourea - 1] = true;
	    }
	}
    }				/* end Thioxo-C */
  if (!
      (is_true_imino_C (a_view) &
       (bond[get_bond (a_view, a_ref) - 1].arom == false)))
    {
      return;
    }				/* end Imino-C */
  if (o_count == 1 && n_count == 1)
    fg[fg_isourea - 1] = true;
  if (s_count == 1 && n_count == 1)
    fg[fg_isothiourea - 1] = true;
  if (n_count == 2)
    fg[fg_guanidine - 1] = true;
}


static void
chk_co2_sp (a_view, a_ref)
     int a_view, a_ref;
{
  int i;
  neighbor_rec nb;
  str2 nb_el;
  int o_count = 0, n_count = 0, s_count = 0;
  int FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  get_neighbors (nb, a_view);
  FORLIM = atom[a_view - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      if (bond[get_bond (a_view, nb[i]) - 1].btype == 'D')
	{
	  strcpy (nb_el, atom[nb[i] - 1].element);
	  if (strcmp (nb_el, "C "))
	    {
	      if (!strcmp (nb_el, "O "))
		o_count++;
	      if (!strcmp (nb_el, "N "))
		n_count++;
	      if (!strcmp (nb_el, "S "))
		s_count++;
	    }
	}
    }
  if (o_count + s_count == 2)	/* new in v0.3b */
    fg[fg_co2_deriv - 1] = true;
  if (o_count == 1 && n_count == 1)
    fg[fg_isocyanate - 1] = true;
  if (s_count == 1 && n_count == 1)
    fg[fg_isothiocyanate - 1] = true;
  if (n_count == 2)
    fg[fg_carbodiimide - 1] = true;
}


static void
chk_triple (a1, a2)
     int a1, a2;
{
  str2 a1_el, a2_el;

  strcpy (a1_el, atom[a1 - 1].element);
  strcpy (a2_el, atom[a2 - 1].element);
  if ((!strcmp (a1_el, "C ") && !strcmp (a2_el, "C ")) &
      (bond[get_bond (a1, a2) - 1].arom == false))
    fg[fg_alkyne - 1] = true;
  if (is_nitrile (a1, a2))
    fg[fg_nitrile - 1] = true;
  if (is_isonitrile (a1, a2))
    fg[fg_isonitrile - 1] = true;
  if (is_cyanate (a1, a2))
    fg[fg_cyanate - 1] = true;
  if (is_thiocyanate (a1, a2))
    fg[fg_thiocyanate - 1] = true;
}


static void
chk_ccx (a_view, a_ref)
     int a_view, a_ref;
{
  int i;
  neighbor_rec nb;
  int oh_count = 0, or_count = 0, n_count = 0;
  int FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  get_neighbors (nb, a_ref);
  FORLIM = atom[a_ref - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      if (bond[get_bond (a_ref, nb[i]) - 1].btype == 'S')
	{
	  if (is_hydroxy (a_ref, nb[i]))
	    oh_count++;
	  if (is_alkoxy (a_ref, nb[i]) || is_aryloxy (a_ref, nb[i]) |
	      is_siloxy (a_ref, nb[i]))
	    or_count++;
	  if (!strcmp (atom[nb[i] - 1].atype, "N3 ") ||
	      !strcmp (atom[nb[i] - 1].atype, "NAM"))
	    n_count++;
	}
    }
  if (oh_count == 1)
    fg[fg_enol - 1] = true;
  if (or_count == 1)
    fg[fg_enolether - 1] = true;
  if (n_count == 1)
    fg[fg_enamine - 1] = true;
  /* new in v0.2f   (regard anything else as an alkene) */
  if (oh_count + or_count + n_count == 0)
    fg[fg_alkene - 1] = true;
}


static void
chk_xccx (a_view, a_ref)
     int a_view, a_ref;
{
  int i;
  neighbor_rec nb;
  int oh_count = 0, or_count = 0, n_count = 0;
  int FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  get_neighbors (nb, a_view);
  FORLIM = atom[a_view - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      if (bond[get_bond (a_view, nb[i]) - 1].btype == 'S')
	{
	  if (is_hydroxy (a_view, nb[i]))
	    oh_count++;
	  if (is_alkoxy (a_view, nb[i]) || is_aryloxy (a_view, nb[i]) |
	      is_siloxy (a_view, nb[i]))
	    or_count++;
	  if (!strcmp (atom[nb[i] - 1].atype, "N3 ") ||
	      !strcmp (atom[nb[i] - 1].atype, "NAM"))
	    n_count++;
	}
    }
  memset (nb, 0, sizeof (neighbor_rec));
  get_neighbors (nb, a_ref);
  FORLIM = atom[a_ref - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      if (bond[get_bond (a_ref, nb[i]) - 1].btype == 'S')
	{
	  if (is_hydroxy (a_ref, nb[i]))
	    oh_count++;
	  if (is_alkoxy (a_ref, nb[i]) || is_aryloxy (a_ref, nb[i]) |
	      is_siloxy (a_ref, nb[i]))
	    or_count++;
	  if (!strcmp (atom[nb[i] - 1].atype, "N3 ") ||
	      !strcmp (atom[nb[i] - 1].atype, "NAM"))
	    n_count++;
	}
    }
  if (oh_count == 2)
    fg[fg_enediol - 1] = true;
  /* new in v0.2f   (regard anything else as an alkene) */
  if (oh_count + or_count + n_count == 0)
    fg[fg_alkene - 1] = true;
}


static void
chk_n_o_dbl (a1, a2)
     int a1, a2;
{
  int i;
  neighbor_rec nb;
  str2 nb_el;
  int or_count = 0, n_count = 0, c_count = 0;
  int b;			/* v0.3j */
  int het_count = 0;		/* v0.3k */
  char bt;			/* v0.3k */
  float bo_sum = 0.0;		/* v0.3k */
  int FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  get_neighbors (nb, a1);
  FORLIM = atom[a1 - 1].neighbor_count;
  /* v0.3k */
  /* v0.3k */
  for (i = 0; i < FORLIM; i++)
    {
      if (nb[i] != a2)
	{
	  b = get_bond (a1, nb[i]);	/* v0.3j */
	  strcpy (nb_el, atom[nb[i] - 1].element);
	  bt = bond[b - 1].btype;	/* v0.3k */
	  if (strcmp (nb_el, "C ") && strcmp (nb_el, "H ")
	      /*&&  strcmp (nb_el, "D ") */  && strcmp (nb_el, "DU")
	      && strcmp (nb_el, "LP") && bond[b - 1].arom == false)
	    /* added 'D ' in v0.3n */
	    het_count++;
	  /* v0.3k: ignore hetero atoms */
	  /* in aromatic rings like isoxazole  */
	  if (bt == 'S')
	    bo_sum += 1.0;
	  if (bt == 'D')
	    bo_sum += 2.0;
	  if (bt == 'A')
	    bo_sum += 1.5;
	  if (!strcmp (nb_el, "O "))
	    or_count++;
	  if (!strcmp (nb_el, "N "))
	    n_count++;
	  if (!strcmp (nb_el, "C ") && bond[b - 1].btype == 'S')	/* v0.3k */
	    c_count++;
	  /* if (is_alkyl(a1,nb[i])) or (is_aryl(a1,nb[i])) then inc(c_count); */
	}
    }
  if (or_count + n_count + c_count == 1 && atom[a1 - 1].neighbor_count == 2)
    {				/* excludes nitro etc. */
      if (or_count == 1)
	fg[fg_nitrite - 1] = true;
      if (c_count == 1)
	fg[fg_nitroso_compound - 1] = true;
      if (n_count == 1)		/* instead of nitrosamine  v0.3j */
	fg[fg_nitroso_compound - 1] = true;
      /*if (n_count = 1) then fg[fg_nitrosamine]   := true;  (* still missing */
    }
  /*if ((c_count > 1) and (or_count = 0) and (n_count = 0)) then */
  /*  begin */
  /*    fg[fg_n_oxide] := true; */
  /*  end; */
  /* new approach in v0.3k */
  if (het_count == 0 && bo_sum > 2)	/* =O does not count! */
    fg[fg_n_oxide - 1] = true;
}


static void
chk_sulfoxide (a1, a2)
     int a1, a2;
{
  int i;
  neighbor_rec nb;
  str2 nb_el;
  int o_count = 0, c_count = 0;
  int FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  get_neighbors (nb, a1);
  FORLIM = atom[a1 - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      strcpy (nb_el, atom[nb[i] - 1].element);
      if (!strcmp (nb_el, "O "))
	o_count++;
      if (is_alkyl (a1, nb[i]) || is_aryl (a1, nb[i]))
	c_count++;
    }
  if (o_count == 1 && c_count == 2)
    fg[fg_sulfoxide - 1] = true;
}


static void
chk_double (a1, a2)
     int a1, a2;
{
  str2 a1_el, a2_el;

  strcpy (a1_el, atom[a1 - 1].element);
  strcpy (a2_el, atom[a2 - 1].element);
  if ((!strcmp (a1_el, "C ") && strcmp (a2_el, "C ")) &
      (bond[get_bond (a1, a2) - 1].arom == false))
    {
      if (hetbond_count (a1) == 2)
	chk_carbonyl_deriv (a1, a2);
      if (hetbond_count (a1) == 3)
	chk_carboxyl_deriv (a1, a2);
      if (hetbond_count (a1) == 4)
	{
	  if (!strcmp (atom[a1 - 1].atype, "C2 "))
	    chk_co2_sp2 (a1, a2);
	  if (!strcmp (atom[a1 - 1].atype, "C1 "))
	    chk_co2_sp (a1, a2);
	}
    }				/* end C=X */
  if ((!strcmp (atom[a1 - 1].atype, "C2 ")
       && !strcmp (atom[a2 - 1].atype,
		   "C2 ")) && (bond[get_bond (a1, a2) - 1].arom == false))
    {
      if ((hetbond_count (a1) == 0) && (hetbond_count (a2) == 2))
	fg[fg_ketene_acetal_deriv - 1] = true;
      if ((hetbond_count (a1) == 0) && (hetbond_count (a2) == 1))
	chk_ccx (a1, a2);
      if ((hetbond_count (a1) == 1) && (hetbond_count (a2) == 1))
	chk_xccx (a1, a2);
      if (((hetbond_count (a1) == 0) && (hetbond_count (a2) == 0)) &&
	  atom[a1 - 1].arom == false && atom[a2 - 1].arom == false)
	fg[fg_alkene - 1] = true;
    }
  if (((!strcmp (a1_el, "N ")
	&& !strcmp (a2_el,
		    "N ")) && (hetbond_count (a1) ==
			       2) && (hetbond_count (a2) ==
				      2) && (bond[get_bond (a1, a2) -
						  1].arom == false))
      && atom[a1 - 1].neighbor_count == 2 && atom[a2 - 1].neighbor_count == 2)
    fg[fg_azo_compound - 1] = true;
  if (!strcmp (a1_el, "N ") && !strcmp (a2_el, "O "))
    chk_n_o_dbl (a1, a2);
  if (!strcmp (a1_el, "S ") && !strcmp (a2_el, "O "))
    chk_sulfoxide (a1, a2);
}


static void
chk_c_hal (a1, a2)
     int a1, a2;
{
  str2 a2_el;

  strcpy (a2_el, atom[a2 - 1].element);
  fg[fg_halogen_deriv - 1] = true;
  if (atom[a1 - 1].arom)
    {
      fg[fg_aryl_halide - 1] = true;
      if (!strcmp (a2_el, "F "))
	fg[fg_aryl_fluoride - 1] = true;
      if (!strcmp (a2_el, "CL"))
	fg[fg_aryl_chloride - 1] = true;
      if (!strcmp (a2_el, "BR"))
	fg[fg_aryl_bromide - 1] = true;
      if (!strcmp (a2_el, "I "))
	fg[fg_aryl_iodide - 1] = true;
      return;
    }
  if ((strcmp (atom[a1 - 1].atype, "C3 ") == 0) && (hetbond_count (a1) <= 2))
    {				/* alkyl halides */
      fg[fg_alkyl_halide - 1] = true;
      if (!strcmp (a2_el, "F "))
	fg[fg_alkyl_fluoride - 1] = true;
      if (!strcmp (a2_el, "CL"))
	fg[fg_alkyl_chloride - 1] = true;
      if (!strcmp (a2_el, "BR"))
	fg[fg_alkyl_bromide - 1] = true;
      if (!strcmp (a2_el, "I "))
	fg[fg_alkyl_iodide - 1] = true;
    }
  if ((strcmp (atom[a1 - 1].atype, "C2 ") == 0) && (hetbond_count (a1) == 3))
    {				/* acyl halides and related compounds */
      if (is_oxo_C (a1))
	{
	  fg[fg_acyl_halide - 1] = true;
	  if (!strcmp (a2_el, "F "))
	    fg[fg_acyl_fluoride - 1] = true;
	  if (!strcmp (a2_el, "CL"))
	    fg[fg_acyl_chloride - 1] = true;
	  if (!strcmp (a2_el, "BR"))
	    fg[fg_acyl_bromide - 1] = true;
	  if (!strcmp (a2_el, "I "))
	    fg[fg_acyl_iodide - 1] = true;
	}
      if (is_thioxo_C (a1))
	fg[fg_thiocarboxylic_acid_deriv - 1] = true;
      if (is_imino_C (a1))
	fg[fg_imidoyl_halide - 1] = true;
    }
  if (!
      ((strcmp (atom[a1 - 1].atype, "C2 ") == 0)
       && (hetbond_count (a1) == 4)))
    /* chloroformates etc. */
    return;
  /* still missing: polyhalogen compounds (-CX2H, -CX3) */
  fg[fg_co2_deriv - 1] = true;
  if (is_oxo_C (a1))
    {
      fg[fg_carbonic_acid_deriv - 1] = true;
      if (is_alkoxycarbonyl (a2, a1) || is_aryloxycarbonyl (a2, a1))
	fg[fg_carbonic_acid_ester_halide - 1] = true;
      if (is_carbamoyl (a2, a1))
	{
	  fg[fg_carbamic_acid_deriv - 1] = true;
	  fg[fg_carbamic_acid_halide - 1] = true;
	}
    }
  if (!is_thioxo_C (a1))
    return;
  fg[fg_thiocarbonic_acid_deriv - 1] = true;
  if (is_alkoxythiocarbonyl (a2, a1) || is_aryloxythiocarbonyl (a2, a1))
    fg[fg_thiocarbonic_acid_ester_halide - 1] = true;
  if (is_thiocarbamoyl (a2, a1))
    {
      fg[fg_thiocarbamic_acid_deriv - 1] = true;
      fg[fg_thiocarbamic_acid_halide - 1] = true;
      /* end of non-aromatic halogen compounds */
    }
}


static void
chk_c_o (a1, a2)
     int a1, a2;
{
  /* ignore heteroaromatic rings (like furan, thiophene, etc.) */
  if (bond[get_bond (a1, a2) - 1].arom == true)
    return;
  if (is_true_alkyl (a2, a1) && is_hydroxy (a1, a2))
    {
      fg[fg_hydroxy - 1] = true;
      fg[fg_alcohol - 1] = true;
      if (atom[a1 - 1].neighbor_count <= 2)
	fg[fg_prim_alcohol - 1] = true;
      if (atom[a1 - 1].neighbor_count == 3)
	fg[fg_sec_alcohol - 1] = true;
      if (atom[a1 - 1].neighbor_count == 4)
	fg[fg_tert_alcohol - 1] = true;
    }
  if (is_aryl (a2, a1) && is_hydroxy (a1, a2))
    {
      fg[fg_hydroxy - 1] = true;
      fg[fg_phenol - 1] = true;
    }
  if (is_true_alkyl (a2, a1) && is_true_alkoxy (a1, a2))
    {
      fg[fg_ether - 1] = true;
      fg[fg_dialkylether - 1] = true;
    }
  if ((is_true_alkyl (a2, a1) && is_aryloxy (a1, a2)) |
      (is_aryl (a2, a1) && is_true_alkoxy (a1, a2)))
    {
      fg[fg_ether - 1] = true;
      fg[fg_alkylarylether - 1] = true;
    }
  if (is_aryl (a2, a1) && is_aryloxy (a1, a2))
    {
      fg[fg_ether - 1] = true;
      fg[fg_diarylether - 1] = true;
    }
  if ((is_true_alkyl (a2, a1) || is_aryl (a2, a1)) && is_alkynyloxy (a1, a2))
    {
      fg[fg_ether - 1] = true;
      ether_generic = true;
    }
  if (is_alkynyl (a2, a1) && is_hydroxy (a1, a2))
    {
      fg[fg_hydroxy - 1] = true;
      hydroxy_generic = true;
    }

}


static void
chk_c_s (a1, a2)
     int a1, a2;
{
  int i;
  neighbor_rec nb;
  str2 nb_el;
  int o_count = 0, oh_count = 0, or_count = 0, n_count = 0, c_count = 0,
    hal_count = 0;
  int FORLIM;

  /* ignore heteroaromatic rings (like furan, thiophene, etc.) */
  if (bond[get_bond (a1, a2) - 1].arom == true)
    return;
  if (is_alkyl (a2, a1) && is_sulfanyl (a1, a2))
    {
      fg[fg_thiol - 1] = true;
      fg[fg_alkylthiol - 1] = true;
    }
  if (is_aryl (a2, a1) && is_sulfanyl (a1, a2))
    {
      fg[fg_thiol - 1] = true;
      fg[fg_arylthiol - 1] = true;
    }
  if (is_true_alkyl (a2, a1) && is_true_alkylsulfanyl (a1, a2))
    fg[fg_thioether - 1] = true;
  if ((is_true_alkyl (a2, a1) && is_arylsulfanyl (a1, a2)) |
      (is_aryl (a2, a1) && is_true_alkylsulfanyl (a1, a2)))
    fg[fg_thioether - 1] = true;
  if (is_aryl (a2, a1) && is_arylsulfanyl (a1, a2))
    fg[fg_thioether - 1] = true;
  /* check for sulfinic/sulfenic acid derivatives */
  memset (nb, 0, sizeof (neighbor_rec));
  get_neighbors (nb, a2);
  FORLIM = atom[a2 - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      strcpy (nb_el, atom[nb[i] - 1].element);
      if (is_alkyl (a2, nb[i]) || is_aryl (a2, nb[i]))
	c_count++;
      if (is_hydroxy (a2, nb[i]))
	oh_count++;
      if (is_alkoxy (a2, nb[i]) || is_aryloxy (a2, nb[i]))
	or_count++;
      if (is_amino (a2, nb[i]) || is_subst_amino (a2, nb[i]))
	n_count++;
      if (!strcmp (nb_el, "F ") || !strcmp (nb_el, "CL") ||
	  !strcmp (nb_el, "BR") || !strcmp (nb_el, "I "))
	hal_count++;
      if (!strcmp (nb_el, "O "))
	o_count++;
    }
  if (c_count != 1)
    return;
  if (atom[a2 - 1].neighbor_count == 3 && o_count - oh_count - or_count == 1)
    {				/* sulfinic acid && derivs */
      fg[fg_sulfinic_acid_deriv - 1] = true;
      if (oh_count == 1)
	fg[fg_sulfinic_acid - 1] = true;
      if (or_count == 1)
	fg[fg_sulfinic_acid_ester - 1] = true;
      if (hal_count == 1)
	fg[fg_sulfinic_acid_halide - 1] = true;
      if (n_count == 1)
	fg[fg_sulfinic_acid_amide - 1] = true;
    }
  if (atom[a2 - 1].neighbor_count != 2 || o_count - oh_count - or_count != 0)
    /* sulfenic acid && derivs */
    return;

  fg[fg_sulfenic_acid_deriv - 1] = true;
  if (oh_count == 1)
    fg[fg_sulfenic_acid - 1] = true;
  if (or_count == 1)
    fg[fg_sulfenic_acid_ester - 1] = true;
  if (hal_count == 1)
    fg[fg_sulfenic_acid_halide - 1] = true;
  if (n_count == 1)
    fg[fg_sulfenic_acid_amide - 1] = true;
}


static void
chk_c_n (a1, a2)
     int a1, a2;
{
  /* ignore heteroaromatic rings (like furan, thiophene, pyrrol, etc.) */
  if (atom[a2 - 1].arom == true)
    return;
  if (is_true_alkyl (a2, a1) && is_amino (a1, a2))
    {
      fg[fg_amine - 1] = true;
      fg[fg_prim_amine - 1] = true;
      fg[fg_prim_aliph_amine - 1] = true;
    }
  if (is_aryl (a2, a1) && is_amino (a1, a2))
    {
      fg[fg_amine - 1] = true;
      fg[fg_prim_amine - 1] = true;
      fg[fg_prim_arom_amine - 1] = true;
    }
  if (is_true_alkyl (a2, a1) && is_true_alkylamino (a1, a2))
    {
      fg[fg_amine - 1] = true;
      fg[fg_sec_amine - 1] = true;
      fg[fg_sec_aliph_amine - 1] = true;
    }
  if (is_aryl (a2, a1) && is_true_alkylamino (a1, a2))
    {
      fg[fg_amine - 1] = true;
      fg[fg_sec_amine - 1] = true;
      fg[fg_sec_mixed_amine - 1] = true;
    }
  if (is_aryl (a2, a1) && is_arylamino (a1, a2))
    {
      fg[fg_amine - 1] = true;
      fg[fg_sec_amine - 1] = true;
      fg[fg_sec_arom_amine - 1] = true;
    }
  if (is_true_alkyl (a2, a1) && is_true_dialkylamino (a1, a2))
    {
      fg[fg_amine - 1] = true;
      fg[fg_tert_amine - 1] = true;
      fg[fg_tert_aliph_amine - 1] = true;
    }
  if ((is_true_alkyl (a2, a1) && is_diarylamino (a1, a2)) |
      (is_aryl (a2, a1) && is_true_dialkylamino (a1, a2)))
    {
      fg[fg_amine - 1] = true;
      fg[fg_tert_amine - 1] = true;
      fg[fg_tert_mixed_amine - 1] = true;
    }
  if (is_aryl (a2, a1) && is_diarylamino (a1, a2))
    {
      fg[fg_amine - 1] = true;
      fg[fg_tert_amine - 1] = true;
      fg[fg_tert_arom_amine - 1] = true;
    }
  if ((is_alkyl (a2, a1) || is_aryl (a2, a1) || is_alkenyl (a2, a1) |
       is_alkynyl (a2, a1)) && is_hydroxylamino (a1, a2) && (is_acyl_gen (a2,
									  a1)
							     == false))
    /* v0.3k */
    fg[fg_hydroxylamine - 1] = true;
  /* v0.3k */
  /* v0.3k  */
  if ((is_alkyl (a2, a1) || is_aryl (a2, a1) || is_acyl (a2, a1) |
       is_alkenyl (a2, a1) || is_alkynyl (a2, a1)) && is_hydrazino (a1, a2))
    fg[fg_hydrazine - 1] = true;
  if ((is_alkyl (a2, a1) || is_aryl (a2, a1) || is_alkenyl (a2, a1) |
       is_alkynyl (a2, a1)) && is_azido (a1, a2))
    /* v0.3k */
    fg[fg_azide - 1] = true;
  if ((is_alkyl (a2, a1) || is_aryl (a2, a1) || is_alkenyl (a2, a1) |
       is_alkynyl (a2, a1)) && is_diazonium (a1, a2))
    /* v0.3k */
    fg[fg_diazonium_salt - 1] = true;
  if ((is_alkyl (a2, a1) || is_aryl (a2, a1) || is_alkenyl (a2, a1) |
       is_alkynyl (a2, a1)) && is_nitro (a1, a2))
    /* v0.3k */
    fg[fg_nitro_compound - 1] = true;
  if (is_alkynyl (a2, a1) &
      (is_amino (a1, a2) || is_C_monosubst_amino (a1, a2) |
       (is_C_disubst_amino (a1, a2) && (!is_acylamino (a1, a2)))))
    {
      fg[fg_amine - 1] = true;
      amine_generic = true;
    }
}


static void
chk_c_c (a1, a2)
     int a1, a2;
{
  int i;
  neighbor_rec nb;
  int oh_count, nhr_count, FORLIM;

  /* ignore aromatic rings */
  if (atom[a2 - 1].arom == true)
    return;
  /*check for 1,2-diols and 1,2-aminoalcoholes */
  if (!strcmp (atom[a1 - 1].atype, "C3 ")
      && !strcmp (atom[a2 - 1].atype, "C3 "))
    {
      if ((hetbond_count (a1) == 1) && (hetbond_count (a2) == 1))
	{
	  oh_count = 0;
	  nhr_count = 0;
	  memset (nb, 0, sizeof (neighbor_rec));
	  get_neighbors (nb, a1);
	  FORLIM = atom[a1 - 1].neighbor_count;
	  for (i = 0; i < FORLIM; i++)
	    {
	      if (nb[i] != a2)
		{
		  if (is_hydroxy (a1, nb[i]))
		    oh_count++;
		  if (is_amino (a1, nb[i]) || is_alkylamino (a1, nb[i]) |
		      is_arylamino (a1, nb[i]))
		    nhr_count++;
		}
	    }
	  memset (nb, 0, sizeof (neighbor_rec));
	  get_neighbors (nb, a2);
	  FORLIM = atom[a2 - 1].neighbor_count;
	  for (i = 0; i < FORLIM; i++)
	    {
	      if (nb[i] != a1)
		{
		  if (is_hydroxy (a2, nb[i]))
		    oh_count++;
		  if (is_amino (a2, nb[i]) || is_alkylamino (a2, nb[i]) |
		      is_arylamino (a2, nb[i]))
		    nhr_count++;
		}
	    }
	  if (oh_count == 2)
	    fg[fg_1_2_diol - 1] = true;
	  if (oh_count == 1 && nhr_count == 1)
	    fg[fg_1_2_aminoalcohol - 1] = true;
	}
    }
  /* check for alpha-aminoacids and alpha-hydroxyacids */
  if (strcmp (atom[a1 - 1].atype, "C3 ")
      || strcmp (atom[a2 - 1].atype, "C2 "))
    return;
  if (!((hetbond_count (a1) == 1) && (hetbond_count (a2) == 3)))
    return;
  oh_count = 0;
  nhr_count = 0;
  memset (nb, 0, sizeof (neighbor_rec));
  get_neighbors (nb, a1);
  FORLIM = atom[a1 - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      if (nb[i] != a2)
	{
	  if (is_hydroxy (a1, nb[i]))
	    oh_count++;
	  if (is_amino (a1, nb[i]) || is_alkylamino (a1, nb[i]) |
	      is_arylamino (a1, nb[i]))
	    nhr_count++;
	}
    }
  memset (nb, 0, sizeof (neighbor_rec));
  get_neighbors (nb, a2);
  FORLIM = atom[a2 - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      if (nb[i] != a1)
	{
	  if (is_hydroxy (a2, nb[i]))
	    oh_count++;
	}
    }
  if ((oh_count == 2) && is_oxo_C (a2))
    fg[fg_alpha_hydroxyacid - 1] = true;
  if ((oh_count == 1 && nhr_count == 1) && is_oxo_C (a2))
    fg[fg_alpha_aminoacid - 1] = true;
}


static void
chk_x_y_single (a_view, a_ref)
     int a_view, a_ref;
{
  if (!strcmp (atom[a_view - 1].atype, "O3 ") &&
      !strcmp (atom[a_ref - 1].atype, "O3 "))
    {
      if (is_hydroxy (a_ref, a_view) || is_hydroxy (a_view, a_ref))
	fg[fg_hydroperoxide - 1] = true;
      if ((is_alkoxy (a_ref, a_view) || is_aryloxy (a_ref, a_view) |
	   is_siloxy (a_ref, a_view)) && (is_alkoxy (a_view,
						     a_ref) |
					  is_aryloxy (a_view,
						      a_ref) |
					  is_siloxy (a_view, a_ref)))
	fg[fg_peroxide - 1] = true;
    }				/* still missing: peracid */
  if (!strcmp (atom[a_view - 1].atype, "S3 ") &&
      !strcmp (atom[a_ref - 1].atype, "S3 "))
    {
      if (atom[a_view - 1].neighbor_count == 2 &&
	  atom[a_ref - 1].neighbor_count == 2)
	fg[fg_disulfide - 1] = true;
    }
  if ((!strcmp (atom[a_view - 1].element, "N ") &&
       !strcmp (atom[a_ref - 1].element,
		"N ")) && (hetbond_count (a_view) ==
			   1) && (hetbond_count (a_ref) == 1))
    {
      /*if ((is_amino(a_ref,a_view)) or  */
      /*    (is_subst_amino(a_ref,a_view)) or */
      /*    (is_acylamino(a_ref,a_view))) and */
      /*   ((is_amino(a_view,a_ref)) or  */
      /*    (is_subst_amino(a_view,a_ref)) or */
      /*    (is_acylamino(a_ref,a_view))) then  */
      if (bond[get_bond (a_view, a_ref) - 1].arom == false)
	fg[fg_hydrazine - 1] = true;
    }
  if (!strcmp (atom[a_view - 1].element, "N ") &&
      !strcmp (atom[a_ref - 1].atype, "O3 "))
    {				/* bond is in "opposite" direction */
      if ((is_alkoxy (a_view, a_ref) || is_aryloxy (a_view, a_ref)) &
	  is_nitro (a_ref, a_view))
	fg[fg_nitrate - 1] = true;
      if ((is_nitro (a_ref, a_view) == false
	   && atom[a_view - 1].arom == false) && (is_amino (a_ref,
							    a_view) |
						  is_subst_amino (a_ref,
								  a_view)) &
	  (is_acylamino (a_ref, a_view) == false))
	fg[fg_hydroxylamine - 1] = true;	/* new in v0.3c */
    }
  if (!strcmp (atom[a_view - 1].element, "S ") &&
      !strcmp (atom[a_ref - 1].element, "O "))
    chk_sulfoxide (a_view, a_ref);
}


static void
chk_single (a1, a2)
     int a1, a2;
{
  str2 a1_el, a2_el;

  strcpy (a1_el, atom[a1 - 1].element);
  strcpy (a2_el, atom[a2 - 1].element);
  if (!strcmp (a1_el, "C ") &&
      (!strcmp (a2_el, "F ") || !strcmp (a2_el, "CL")
       || !strcmp (a2_el, "BR") || !strcmp (a2_el, "I ")))
    chk_c_hal (a1, a2);
  if (!strcmp (a1_el, "C ") && !strcmp (a2_el, "O "))
    chk_c_o (a1, a2);
  if (!strcmp (a1_el, "C ") && !strcmp (a2_el, "S "))
    chk_c_s (a1, a2);
  if (!strcmp (a1_el, "C ") && !strcmp (a2_el, "N "))
    chk_c_n (a1, a2);
  if ((strcmp (a1_el, "C ") == 0) && atom[a2 - 1].metal && (is_cyano_c (a1) ==
							    false))
    {
      fg[fg_organometallic - 1] = true;
      if (!strcmp (a2_el, "LI"))
	fg[fg_organolithium - 1] = true;
      if (!strcmp (a2_el, "MG"))
	fg[fg_organomagnesium - 1] = true;
    }
  if (!strcmp (a1_el, "C ") && !strcmp (a2_el, "C "))
    chk_c_c (a1, a2);
  if (strcmp (a1_el, "C ") && strcmp (a2_el, "C "))
    chk_x_y_single (a1, a2);
}


static void
chk_carbonyl_deriv_sp3 (a_ref)
     int a_ref;
{
  int i;
  neighbor_rec nb;
  int oh_count = 0, or_count = 0, n_count = 0, sh_count = 0, sr_count = 0;
  int FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  get_neighbors (nb, a_ref);
  FORLIM = atom[a_ref - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      if (is_hydroxy (a_ref, nb[i]))
	oh_count++;
      if (is_alkoxy (a_ref, nb[i]) || is_aryloxy (a_ref, nb[i]) |
	  is_alkenyloxy (a_ref, nb[i]) || is_alkynyloxy (a_ref, nb[i]))
	or_count++;
      if (is_sulfanyl (a_ref, nb[i]))
	sh_count++;
      if (is_alkylsulfanyl (a_ref, nb[i]) || is_arylsulfanyl (a_ref, nb[i]) |
	  is_alkenylsulfanyl (a_ref, nb[i]) || is_alkynylsulfanyl (a_ref,
								   nb[i]))
	sr_count++;
      if (!strcmp (atom[nb[i] - 1].atype, "N3 ") ||
	  !strcmp (atom[nb[i] - 1].atype, "NAM"))
	n_count++;
    }
  if (oh_count == 2)
    fg[fg_carbonyl_hydrate - 1] = true;
  if (oh_count == 1 && or_count == 1)
    fg[fg_hemiacetal - 1] = true;
  if (or_count == 2)
    fg[fg_acetal - 1] = true;
  if ((oh_count == 1 || or_count == 1) && n_count == 1)
    fg[fg_hemiaminal - 1] = true;
  if (n_count == 2)
    fg[fg_aminal - 1] = true;
  if ((sh_count == 1 || sr_count == 1) && n_count == 1)
    fg[fg_thiohemiaminal - 1] = true;
  if (sr_count == 2 || or_count == 1 && sr_count == 1)
    fg[fg_thioacetal - 1] = true;
}


static void
chk_carboxyl_deriv_sp3 (a_ref)
     int a_ref;
{
  int i;
  neighbor_rec nb;
  int or_count = 0, oh_count = 0, n_count = 0;	/* oh_count new in v0.3c */
  int electroneg_count = 0;	/* new in v0.3j */
  int hal_count = 0;
  str2 nb_el;
  int FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  get_neighbors (nb, a_ref);
  FORLIM = atom[a_ref - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      strcpy (nb_el, atom[nb[i] - 1].element);	/* v0.3j */
      if (is_electroneg (nb_el))
	electroneg_count++;
      if (!strcmp (nb_el, "F ") || !strcmp (nb_el, "CL") ||
	  !strcmp (nb_el, "BR") || !strcmp (nb_el, "I ")
	  || !strcmp (nb_el, "AT"))
	hal_count++;
      if (is_alkoxy (a_ref, nb[i]) || is_aryloxy (a_ref, nb[i]) |
	  is_siloxy (a_ref, nb[i]))
	or_count++;
      if (is_hydroxy (a_ref, nb[i]))	/* new in v0.3c    */
	oh_count++;
      if (!strcmp (atom[nb[i] - 1].atype, "N3 ") ||
	  !strcmp (atom[nb[i] - 1].atype, "NAM"))
	n_count++;
    }
  /*if (or_count + n_count > 1) then fg[fg_orthocarboxylic_acid_deriv] := true;  (* until v0.3i */
  if (electroneg_count == 3 && hal_count < 3)	/* v0.3j */
    fg[fg_orthocarboxylic_acid_deriv - 1] = true;
  if (or_count == 3)
    fg[fg_carboxylic_acid_orthoester - 1] = true;
  if (or_count == 2 && n_count == 1)
    fg[fg_carboxylic_acid_amide_acetal - 1] = true;
  if (oh_count > 0 && oh_count + or_count + n_count == 3)	/* new in v0.3c */
    fg[fg_orthocarboxylic_acid_deriv - 1] = true;
}


static void
chk_anhydride (a_ref)
     int a_ref;
{
  int i;
  neighbor_rec nb;
  int acyl_count = 0;
  int FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  get_neighbors (nb, a_ref);
  FORLIM = atom[a_ref - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      if (is_acyl (a_ref, nb[i]) || is_carbamoyl (a_ref, nb[i]))
	acyl_count++;
    }
  if (acyl_count == 2 && !strcmp (atom[a_ref - 1].atype, "O3 "))
    {
      fg[fg_carboxylic_acid_deriv - 1] = true;
      fg[fg_carboxylic_acid_anhydride - 1] = true;
    }
}


static void
chk_imide (a_ref)
     int a_ref;
{
  int i;
  neighbor_rec nb;
  int acyl_count = 0;
  int FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  get_neighbors (nb, a_ref);
  FORLIM = atom[a_ref - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      if (is_acyl_gen (a_ref, nb[i]) || is_carbamoyl (a_ref, nb[i]))	/* v0.3j */
	acyl_count++;
    }
  if (acyl_count < 2 || strcmp (atom[a_ref - 1].element, "N "))
    /* v0.3j: accept also N-acyl-imides */
    return;
  fg[fg_carboxylic_acid_deriv - 1] = true;
  fg[fg_carboxylic_acid_imide - 1] = true;
  if (atom[a_ref - 1].neighbor_count == 2)
    fg[fg_carboxylic_acid_unsubst_imide - 1] = true;
  if (atom[a_ref - 1].neighbor_count == 3)
    fg[fg_carboxylic_acid_subst_imide - 1] = true;
}


static void
chk_12diphenol (a_view, a_ref)
     int a_view, a_ref;
{
  int i;
  neighbor_rec nb;
  int oh_count = 0;
  int FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  get_neighbors (nb, a_view);
  FORLIM = atom[a_view - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      if (bond[get_bond (a_view, nb[i]) - 1].btype == 'S')
	{
	  if (is_hydroxy (a_view, nb[i]))
	    oh_count++;
	}
    }
  memset (nb, 0, sizeof (neighbor_rec));
  get_neighbors (nb, a_ref);
  FORLIM = atom[a_ref - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    {
      if (bond[get_bond (a_ref, nb[i]) - 1].btype == 'S')
	{
	  if (is_hydroxy (a_ref, nb[i]))
	    oh_count++;
	}
    }
  if (oh_count == 2)
    fg[fg_1_2_diphenol - 1] = true;
}


static void
chk_arom_fg (a1, a2)
     int a1, a2;
{
  if ((hetbond_count (a1) == 1) && (hetbond_count (a2) == 1))
    chk_12diphenol (a1, a2);
}


static boolean
is_arene (r_id)
     int r_id;
{
  int i, j;
  boolean r = true;
  ringpath_type testring;
  int ring_size, a_prev, a_ref;

/* p2c: checkmol.pas: Note: Eliminated unused assignment statement [338] */
  if (r_id < 1 || r_id > n_rings)
    return false;
  memset (testring, 0, sizeof (ringpath_type));
  ring_size = ringprop[r_id - 1].size;	/* v0.3j */
  /*for j := 1 to max_ringsize do if ring^[r_id,j] > 0 then testring[j] := ring^[r_id,j]; */
  for (j = 0; j < ring_size; j++)	/* v0.3j */
    testring[j] = ring[r_id - 1][j];
  /*ring_size := path_length(testring); */
  if (ring_size <= 2)
    return false;
  a_prev = testring[ring_size - 1];
  for (i = 0; i < ring_size; i++)
    {
      a_ref = testring[i];
      if (bond[get_bond (a_prev, a_ref) - 1].arom == false)
	r = false;
      a_prev = a_ref;
    }
  return r;
}


static boolean
is_heterocycle (r_id)
     int r_id;
{
  int i, j;
  boolean r = false;
  ringpath_type testring;
  int ring_size, a_ref;

  if (r_id < 1 || r_id > n_rings)
    return false;
  memset (testring, 0, sizeof (ringpath_type));
  ring_size = ringprop[r_id - 1].size;	/* v0.3j */
  /*for j := 1 to max_ringsize do if ring^[r_id,j] > 0 then testring[j] := ring^[r_id,j]; */
  for (j = 0; j < ring_size; j++)	/* v0.3j */
    testring[j] = ring[r_id - 1][j];
  /*ring_size := path_length(testring); */
  if (ring_size <= 2)
    return false;
  for (i = 0; i < ring_size; i++)
    {
      a_ref = testring[i];
      if (strcmp (atom[a_ref - 1].element, "C "))
	r = true;
    }
  return r;
}


static void
chk_oxo_thioxo_imino_hetarene (r_id)
     int r_id;
{
  int i, j;
  ringpath_type testring;
  int ring_size, a_ref;

  if (r_id < 1 || r_id > n_rings)
    return;
  memset (testring, 0, sizeof (ringpath_type));
  ring_size = ringprop[r_id - 1].size;	/* v0.3j */
  /*for j := 1 to max_ringsize do if ring^[r_id,j] > 0 then testring[j] := ring^[r_id,j]; */
  for (j = 0; j < ring_size; j++)	/* v0.3j */
    testring[j] = ring[r_id - 1][j];
  /*ring_size := path_length(testring); */
  /*if (is_arene(r_id)) and (odd(ring_size) = false) then */
  if (!is_arene (r_id))		/* v0.3j */
    return;
  for (i = 0; i < ring_size; i++)
    {
      a_ref = testring[i];
      if (is_oxo_C (a_ref))
	fg[fg_oxohetarene - 1] = true;
      if (is_thioxo_C (a_ref))
	fg[fg_thioxohetarene - 1] = true;
      if (is_true_exocyclic_imino_C (a_ref, r_id))	/* v0.3j */
	fg[fg_iminohetarene - 1] = true;
    }
}


static void
chk_ion (a_ref)
     int a_ref;
{
  int i;
  neighbor_rec nb;
  int charge, FORLIM;

  memset (nb, 0, sizeof (neighbor_rec));
  get_neighbors (nb, a_ref);
  charge = atom[a_ref - 1].formal_charge;
  if (charge == 0)
    /* check if charge is neutralized by an adjacent opposite charge */
    return;
  FORLIM = atom[a_ref - 1].neighbor_count;
  for (i = 0; i < FORLIM; i++)
    charge += atom[nb[i] - 1].formal_charge;
  if (charge > 0)
    fg[fg_cation - 1] = true;
  if (charge < 0)
    fg[fg_anion - 1] = true;
}


static void
chk_functionalgroups ()
{
  int i, a1, a2;
  char bt;
  int pos_chg = 0, neg_chg = 0;
  int FORLIM;

  if (n_atoms < 1 || n_bonds < 1)
    return;
  FORLIM = n_atoms;
  for (i = 1; i <= FORLIM; i++)
    {				/* a few groups are best discovered in the atom list */
      if (!strcmp (atom[i - 1].atype, "SO2"))
	chk_so2_deriv (i);
      /*if (atom^[i].atype = 'SO ') then fg[fg_sulfoxide] := true;  (* do another check in the bond list!! */
      if (!strcmp (atom[i - 1].element, "P "))
	chk_p_deriv (i);
      if (!strcmp (atom[i - 1].element, "B "))
	chk_b_deriv (i);
      if (!strcmp (atom[i - 1].atype, "N3+") || atom[i - 1].formal_charge > 0)
	chk_ammon (i);
      if ((strcmp (atom[i - 1].atype, "C3 ") == 0)
	  && (hetbond_count (i) == 2))
	chk_carbonyl_deriv_sp3 (i);
      if ((strcmp (atom[i - 1].atype, "C3 ") == 0)
	  && (hetbond_count (i) == 3))
	chk_carboxyl_deriv_sp3 (i);
      if (!strcmp (atom[i - 1].atype, "O3 ")
	  && atom[i - 1].neighbor_count == 2)
	chk_anhydride (i);
      if ((!strcmp (atom[i - 1].atype, "N3 ")
	   || !strcmp (atom[i - 1].atype, "NAM"))
	  && atom[i - 1].neighbor_count >= 2)
	chk_imide (i);
      if (atom[i - 1].formal_charge > 0)
	pos_chg += atom[i - 1].formal_charge;
      if (atom[i - 1].formal_charge < 0)
	neg_chg += atom[i - 1].formal_charge;
      chk_ion (i);
    }
  FORLIM = n_bonds;
  for (i = 0; i < FORLIM; i++)
    {				/* most groups are best discovered in the bond list */
      a1 = bond[i].a1;
      a2 = bond[i].a2;
      bt = bond[i].btype;
      if (atom[a1 - 1].heavy && atom[a2 - 1].heavy)
	{
	  orient_bond (&a1, &a2);
	  if (bt == 'T')
	    chk_triple (a1, a2);
	  if (bt == 'D')
	    chk_double (a1, a2);
	  if (bt == 'S')
	    chk_single (a1, a2);
	  if (bond[i].arom)
	    chk_arom_fg (a1, a2);
	}
    }
  if (n_rings > 0)
    {
      FORLIM = n_rings;
      for (i = 1; i <= FORLIM; i++)
	{
	  chk_oxo_thioxo_imino_hetarene (i);
	  if (is_arene (i))
	    fg[fg_aromatic - 1] = true;
	  if (is_heterocycle (i))
	    fg[fg_heterocycle - 1] = true;
	}
    }
  if (pos_chg + neg_chg > 0)
    fg[fg_cation - 1] = true;
  if (pos_chg + neg_chg < 0)
    fg[fg_anion - 1] = true;
}


static void
write_fg_text ()
{
  if (fg[fg_cation - 1])
    printf ("cation\n");
  if (fg[fg_anion - 1])
    printf ("anion\n");
  /*  if fg[fg_carbonyl]                       then writeln('carbonyl compound'); */
  if (fg[fg_aldehyde - 1])
    printf ("aldehyde\n");
  if (fg[fg_ketone - 1])
    printf ("ketone\n");
  /*  if fg[fg_thiocarbonyl]                   then writeln('thiocarbonyl compound'); */
  if (fg[fg_thioaldehyde - 1])
    printf ("thioaldehyde\n");
  if (fg[fg_thioketone - 1])
    printf ("thioketone\n");
  if (fg[fg_imine - 1])
    printf ("imine\n");
  if (fg[fg_hydrazone - 1])
    printf ("hydrazone\n");
  if (fg[fg_semicarbazone - 1])
    printf ("semicarbazone\n");
  if (fg[fg_thiosemicarbazone - 1])
    printf ("thiosemicarbazone\n");
  if (fg[fg_oxime - 1])
    printf ("oxime\n");
  if (fg[fg_oxime_ether - 1])
    printf ("oxime ether\n");
  if (fg[fg_ketene - 1])
    printf ("ketene\n");
  if (fg[fg_ketene_acetal_deriv - 1])
    printf ("ketene acetal or derivative\n");
  if (fg[fg_carbonyl_hydrate - 1])
    printf ("carbonyl hydrate\n");
  if (fg[fg_hemiacetal - 1])
    printf ("hemiacetal\n");
  if (fg[fg_acetal - 1])
    printf ("acetal\n");
  if (fg[fg_hemiaminal - 1])
    printf ("hemiaminal\n");
  if (fg[fg_aminal - 1])
    printf ("aminal\n");
  if (fg[fg_thiohemiaminal - 1])
    printf ("hemithioaminal\n");
  if (fg[fg_thioacetal - 1])
    printf ("thioacetal\n");
  if (fg[fg_enamine - 1])
    printf ("enamine\n");
  if (fg[fg_enol - 1])
    printf ("enol\n");
  if (fg[fg_enolether - 1])
    printf ("enol ether\n");
  if (fg[fg_hydroxy - 1] && hydroxy_generic)
    printf ("hydroxy compound\n");
  /*  if fg[fg_alcohol]                        then writeln('alcohol'); */
  if (fg[fg_prim_alcohol - 1])
    printf ("primary alcohol\n");
  if (fg[fg_sec_alcohol - 1])
    printf ("secondary alcohol\n");
  if (fg[fg_tert_alcohol - 1])
    printf ("tertiary alcohol\n");
  if (fg[fg_1_2_diol - 1])
    printf ("1,2-diol\n");
  if (fg[fg_1_2_aminoalcohol - 1])
    printf ("1,2-aminoalcohol\n");
  if (fg[fg_phenol - 1])
    printf ("phenol or hydroxyhetarene\n");
  if (fg[fg_1_2_diphenol - 1])
    printf ("1,2-diphenol\n");
  if (fg[fg_enediol - 1])
    printf ("enediol\n");
  if (fg[fg_ether - 1] && ether_generic)
    printf ("ether\n");
  if (fg[fg_dialkylether - 1])
    printf ("dialkyl ether\n");
  if (fg[fg_alkylarylether - 1])
    printf ("alkyl aryl ether \n");
  if (fg[fg_diarylether - 1])
    printf ("diaryl ether\n");
  if (fg[fg_thioether - 1])
    printf ("thioether\n");
  if (fg[fg_disulfide - 1])
    printf ("disulfide\n");
  if (fg[fg_peroxide - 1])
    printf ("peroxide\n");
  if (fg[fg_hydroperoxide - 1])
    printf ("hydroperoxide \n");
  if (fg[fg_hydrazine - 1])
    printf ("hydrazine derivative\n");
  if (fg[fg_hydroxylamine - 1])
    printf ("hydroxylamine\n");
  if (fg[fg_amine - 1] && amine_generic)
    printf ("amine\n");
  if (fg[fg_prim_amine - 1])
    printf ("primary amine\n");
  if (fg[fg_prim_aliph_amine - 1])
    printf ("primary aliphatic amine (alkylamine)\n");
  if (fg[fg_prim_arom_amine - 1])
    printf ("primary aromatic amine\n");
  if (fg[fg_sec_amine - 1])
    printf ("secondary amine\n");
  if (fg[fg_sec_aliph_amine - 1])
    printf ("secondary aliphatic amine (dialkylamine)\n");
  if (fg[fg_sec_mixed_amine - 1])
    printf ("secondary aliphatic/aromatic amine (alkylarylamine)\n");
  if (fg[fg_sec_arom_amine - 1])
    printf ("secondary aromatic amine (diarylamine)\n");
  if (fg[fg_tert_amine - 1])
    printf ("tertiary amine\n");
  if (fg[fg_tert_aliph_amine - 1])
    printf ("tertiary aliphatic amine (trialkylamine)\n");
  if (fg[fg_tert_mixed_amine - 1])
    printf ("tertiary aliphatic/aromatic amine (alkylarylamine)\n");
  if (fg[fg_tert_arom_amine - 1])
    printf ("tertiary aromatic amine (triarylamine)\n");
  if (fg[fg_quart_ammonium - 1])
    printf ("quaternary ammonium salt\n");
  if (fg[fg_n_oxide - 1])
    printf ("N-oxide\n");
  /* new in v0.2f */
  if (fg[fg_halogen_deriv - 1])
    {
      if (!fg[fg_alkyl_halide - 1] && !fg[fg_aryl_halide - 1] &&
	  !fg[fg_acyl_halide - 1])
	printf ("halogen derivative\n");
    }
  /*  if fg[fg_alkyl_halide]                   then writeln('alkyl halide'); */
  if (fg[fg_alkyl_fluoride - 1])
    printf ("alkyl fluoride\n");
  if (fg[fg_alkyl_chloride - 1])
    printf ("alkyl chloride\n");
  if (fg[fg_alkyl_bromide - 1])
    printf ("alkyl bromide\n");
  if (fg[fg_alkyl_iodide - 1])
    printf ("alkyl iodide\n");
  /*  if fg[fg_aryl_halide]                    then writeln('aryl halide'); */
  if (fg[fg_aryl_fluoride - 1])
    printf ("aryl fluoride\n");
  if (fg[fg_aryl_chloride - 1])
    printf ("aryl chloride\n");
  if (fg[fg_aryl_bromide - 1])
    printf ("aryl bromide\n");
  if (fg[fg_aryl_iodide - 1])
    printf ("aryl iodide\n");
  if (fg[fg_organometallic - 1])
    printf ("organometallic compound\n");
  if (fg[fg_organolithium - 1])
    printf ("organolithium compound\n");
  if (fg[fg_organomagnesium - 1])
    printf ("organomagnesium compound\n");
  /*  if fg[fg_carboxylic_acid_deriv]          then writeln('carboxylic acid derivative'); */
  if (fg[fg_carboxylic_acid - 1])
    printf ("carboxylic acid\n");
  if (fg[fg_carboxylic_acid_salt - 1])
    printf ("carboxylic acid salt\n");
  if (fg[fg_carboxylic_acid_ester - 1])
    printf ("carboxylic acid ester\n");
  if (fg[fg_lactone - 1])
    printf ("lactone\n");
  /*  if fg[fg_carboxylic_acid_amide]          then writeln('carboxylic acid amide'); */
  if (fg[fg_carboxylic_acid_prim_amide - 1])
    printf ("primary carboxylic acid amide\n");
  if (fg[fg_carboxylic_acid_sec_amide - 1])
    printf ("secondary carboxylic acid amide\n");
  if (fg[fg_carboxylic_acid_tert_amide - 1])
    printf ("tertiary carboxylic acid amide\n");
  if (fg[fg_lactam - 1])
    printf ("lactam\n");
  if (fg[fg_carboxylic_acid_hydrazide - 1])
    printf ("carboxylic acid hydrazide\n");
  if (fg[fg_carboxylic_acid_azide - 1])
    printf ("carboxylic acid azide\n");
  if (fg[fg_hydroxamic_acid - 1])
    printf ("hydroxamic acid\n");
  if (fg[fg_carboxylic_acid_amidine - 1])
    printf ("carboxylic acid amidine\n");
  if (fg[fg_carboxylic_acid_amidrazone - 1])
    printf ("carboxylic acid amidrazone\n");
  if (fg[fg_nitrile - 1])
    printf ("carbonitrile\n");
  /*  if fg[fg_acyl_halide]                    then writeln('acyl halide'); */
  if (fg[fg_acyl_fluoride - 1])
    printf ("acyl fluoride\n");
  if (fg[fg_acyl_chloride - 1])
    printf ("acyl chloride\n");
  if (fg[fg_acyl_bromide - 1])
    printf ("acyl bromide\n");
  if (fg[fg_acyl_iodide - 1])
    printf ("acyl iodide\n");
  if (fg[fg_acyl_cyanide - 1])
    printf ("acyl cyanide\n");
  if (fg[fg_imido_ester - 1])
    printf ("imido ester\n");
  if (fg[fg_imidoyl_halide - 1])
    printf ("imidoyl halide\n");
  /*  if fg[fg_thiocarboxylic_acid_deriv]      then writeln('thiocarboxylic acid derivative'); */
  if (fg[fg_thiocarboxylic_acid - 1])
    printf ("thiocarboxylic acid\n");
  if (fg[fg_thiocarboxylic_acid_ester - 1])
    printf ("thiocarboxylic acid ester\n");
  if (fg[fg_thiolactone - 1])
    printf ("thiolactone\n");
  if (fg[fg_thiocarboxylic_acid_amide - 1])
    printf ("thiocarboxylic acid amide\n");
  if (fg[fg_thiolactam - 1])
    printf ("thiolactam\n");
  if (fg[fg_imido_thioester - 1])
    printf ("imidothioester\n");
  if (fg[fg_oxohetarene - 1])
    printf ("oxo(het)arene\n");
  if (fg[fg_thioxohetarene - 1])
    printf ("thioxo(het)arene\n");
  if (fg[fg_iminohetarene - 1])
    printf ("imino(het)arene\n");
  if (fg[fg_orthocarboxylic_acid_deriv - 1])
    printf ("orthocarboxylic acid derivative\n");
  if (fg[fg_carboxylic_acid_orthoester - 1])
    printf ("orthoester\n");
  if (fg[fg_carboxylic_acid_amide_acetal - 1])
    printf ("amide acetal\n");
  if (fg[fg_carboxylic_acid_anhydride - 1])
    printf ("carboxylic acid anhydride\n");
  /*  if fg[fg_carboxylic_acid_imide]          then writeln('carboxylic acid imide'); */
  if (fg[fg_carboxylic_acid_unsubst_imide - 1])
    printf ("carboxylic acid imide, N-unsubstituted\n");
  if (fg[fg_carboxylic_acid_subst_imide - 1])
    printf ("carboxylic acid imide, N-substituted\n");
  if (fg[fg_co2_deriv - 1])
    printf ("CO2 derivative (general)\n");
  if (fg[fg_carbonic_acid_deriv - 1] &&
      !(fg[fg_carbonic_acid_monoester - 1]
	|| fg[fg_carbonic_acid_diester - 1]
	|| fg[fg_carbonic_acid_ester_halide - 1]))
    /* changed in v0.3c */
    printf ("carbonic acid derivative\n");
  if (fg[fg_carbonic_acid_monoester - 1])
    printf ("carbonic acid monoester\n");
  if (fg[fg_carbonic_acid_diester - 1])
    printf ("carbonic acid diester\n");
  if (fg[fg_carbonic_acid_ester_halide - 1])
    printf ("carbonic acid ester halide (alkyl/aryl haloformate)\n");
  if (fg[fg_thiocarbonic_acid_deriv - 1])
    printf ("thiocarbonic acid derivative\n");
  if (fg[fg_thiocarbonic_acid_monoester - 1])
    printf ("thiocarbonic acid monoester\n");
  if (fg[fg_thiocarbonic_acid_diester - 1])
    printf ("thiocarbonic acid diester\n");
  if (fg[fg_thiocarbonic_acid_ester_halide - 1])
    printf ("thiocarbonic acid ester halide (alkyl/aryl halothioformate)\n");
  if (fg[fg_carbamic_acid_deriv - 1]
      && !(fg[fg_carbamic_acid - 1] || fg[fg_carbamic_acid_ester - 1]
	   || fg[fg_carbamic_acid_halide - 1]))
    /* changed in v0.3c */
    printf ("carbamic acid derivative\n");
  if (fg[fg_carbamic_acid - 1])
    printf ("carbamic acid\n");
  if (fg[fg_carbamic_acid_ester - 1])
    printf ("carbamic acid ester (urethane)\n");
  if (fg[fg_carbamic_acid_halide - 1])
    printf ("carbamic acid halide (haloformic acid amide)\n");
  if (fg[fg_thiocarbamic_acid_deriv - 1] &&
      !(fg[fg_thiocarbamic_acid - 1] || fg[fg_thiocarbamic_acid_ester - 1]
	|| fg[fg_thiocarbamic_acid_halide - 1]))
    /* changed in v0.3c */
    printf ("thiocarbamic acid derivative\n");
  if (fg[fg_thiocarbamic_acid - 1])
    printf ("thiocarbamic acid\n");
  if (fg[fg_thiocarbamic_acid_ester - 1])
    printf ("thiocarbamic acid ester\n");
  if (fg[fg_thiocarbamic_acid_halide - 1])
    printf ("thiocarbamic acid halide (halothioformic acid amide)\n");
  if (fg[fg_urea - 1])
    printf ("urea\n");
  if (fg[fg_isourea - 1])
    printf ("isourea\n");
  if (fg[fg_thiourea - 1])
    printf ("thiourea\n");
  if (fg[fg_isothiourea - 1])
    printf ("isothiourea\n");
  if (fg[fg_guanidine - 1])
    printf ("guanidine\n");
  if (fg[fg_semicarbazide - 1])
    printf ("semicarbazide\n");
  if (fg[fg_thiosemicarbazide - 1])
    printf ("thiosemicarbazide\n");
  if (fg[fg_azide - 1])
    printf ("azide\n");
  if (fg[fg_azo_compound - 1])
    printf ("azo compound\n");
  if (fg[fg_diazonium_salt - 1])
    printf ("diazonium salt\n");
  if (fg[fg_isonitrile - 1])
    printf ("isonitrile\n");
  if (fg[fg_cyanate - 1])
    printf ("cyanate\n");
  if (fg[fg_isocyanate - 1])
    printf ("isocyanate\n");
  if (fg[fg_thiocyanate - 1])
    printf ("thiocyanate\n");
  if (fg[fg_isothiocyanate - 1])
    printf ("isothiocyanate\n");
  if (fg[fg_carbodiimide - 1])
    printf ("carbodiimide\n");
  if (fg[fg_nitroso_compound - 1])
    printf ("nitroso compound\n");
  if (fg[fg_nitro_compound - 1])
    printf ("nitro compound\n");
  if (fg[fg_nitrite - 1])
    printf ("nitrite\n");
  if (fg[fg_nitrate - 1])
    printf ("nitrate\n");
  /*  if fg[fg_sulfuric_acid_deriv]            then writeln('sulfuric acid derivative'); */
  if (fg[fg_sulfuric_acid - 1])
    printf ("sulfuric acid\n");
  if (fg[fg_sulfuric_acid_monoester - 1])
    printf ("sulfuric acid monoester\n");
  if (fg[fg_sulfuric_acid_diester - 1])
    printf ("sulfuric acid diester\n");
  if (fg[fg_sulfuric_acid_amide_ester - 1])
    printf ("sulfuric acid amide ester\n");
  if (fg[fg_sulfuric_acid_amide - 1])
    printf ("sulfuric acid amide\n");
  if (fg[fg_sulfuric_acid_diamide - 1])
    printf ("sulfuric acid diamide\n");
  if (fg[fg_sulfuryl_halide - 1])
    printf ("sulfuryl halide\n");
  /*  if fg[fg_sulfonic_acid_deriv]            then writeln('sulfonic acid derivative '); */
  if (fg[fg_sulfonic_acid - 1])
    printf ("sulfonic acid\n");
  if (fg[fg_sulfonic_acid_ester - 1])
    printf ("sulfonic acid ester\n");
  if (fg[fg_sulfonamide - 1])
    printf ("sulfonamide\n");
  if (fg[fg_sulfonyl_halide - 1])
    printf ("sulfonyl halide\n");
  if (fg[fg_sulfone - 1])
    printf ("sulfone\n");
  if (fg[fg_sulfoxide - 1])
    printf ("sulfoxide\n");
  /*  if fg[fg_sulfinic_acid_deriv]            then writeln('sulfinic acid derivative'); */
  if (fg[fg_sulfinic_acid - 1])
    printf ("sulfinic acid\n");
  if (fg[fg_sulfinic_acid_ester - 1])
    printf ("sulfinic acid ester\n");
  if (fg[fg_sulfinic_acid_halide - 1])
    printf ("sulfinic acid halide\n");
  if (fg[fg_sulfinic_acid_amide - 1])
    printf ("sulfinic acid amide\n");
  /*  if fg[fg_sulfenic_acid_deriv]            then writeln('sulfenic acid derivative'); */
  if (fg[fg_sulfenic_acid - 1])
    printf ("sulfenic acid\n");
  if (fg[fg_sulfenic_acid_ester - 1])
    printf ("sulfenic acid ester\n");
  if (fg[fg_sulfenic_acid_halide - 1])
    printf ("sulfenic acid halide\n");
  if (fg[fg_sulfenic_acid_amide - 1])
    printf ("sulfenic acid amide\n");
  if (fg[fg_thiol - 1])
    printf ("thiol (sulfanyl compound)\n");
  if (fg[fg_alkylthiol - 1])
    printf ("alkylthiol\n");
  if (fg[fg_arylthiol - 1])
    printf ("arylthiol\n");
  /*  if fg[fg_phosphoric_acid_deriv]          then writeln('phosphoric acid derivative'); */
  if (fg[fg_phosphoric_acid - 1])
    printf ("phosphoric acid\n");
  if (fg[fg_phosphoric_acid_ester - 1])
    printf ("phosphoric acid ester\n");
  if (fg[fg_phosphoric_acid_halide - 1])
    printf ("phosphoric acid halide\n");
  if (fg[fg_phosphoric_acid_amide - 1])
    printf ("phosphoric acid amide\n");
  /*  if fg[fg_thiophosphoric_acid_deriv]      then writeln('thiophosphoric acid derivative'); */
  if (fg[fg_thiophosphoric_acid - 1])
    printf ("thiophosphoric acid\n");
  if (fg[fg_thiophosphoric_acid_ester - 1])
    printf ("thiophosphoric acid ester\n");
  if (fg[fg_thiophosphoric_acid_halide - 1])
    printf ("thiophosphoric acid halide\n");
  if (fg[fg_thiophosphoric_acid_amide - 1])
    printf ("thiophosphoric acid amide\n");
  if (fg[fg_phosphonic_acid_deriv - 1])
    printf ("phosphonic acid derivative \n");
  if (fg[fg_phosphonic_acid - 1])
    printf ("phosphonic acid\n");
  if (fg[fg_phosphonic_acid_ester - 1])
    printf ("phosphonic acid ester\n");
  if (fg[fg_phosphine - 1])
    printf ("phosphine\n");
  if (fg[fg_phosphinoxide - 1])
    printf ("phosphine oxide\n");
  if (fg[fg_boronic_acid_deriv - 1])
    printf ("boronic acid derivative\n");
  if (fg[fg_boronic_acid - 1])
    printf ("boronic acid\n");
  if (fg[fg_boronic_acid_ester - 1])
    printf ("boronic acid ester\n");
  if (fg[fg_alkene - 1])
    printf ("alkene\n");
  if (fg[fg_alkyne - 1])
    printf ("alkyne\n");
  if (fg[fg_aromatic - 1])
    printf ("aromatic compound\n");
  if (fg[fg_heterocycle - 1])
    printf ("heterocyclic compound\n");
  if (fg[fg_alpha_aminoacid - 1])
    printf ("alpha-aminoacid\n");
  if (fg[fg_alpha_hydroxyacid - 1])
    printf ("alpha-hydroxyacid\n");
}


static void
write_fg_text_de ()
{
  if (fg[fg_cation - 1])
    printf ("Kation\n");
  if (fg[fg_anion - 1])
    printf ("Anion\n");
  /*  if fg[fg_carbonyl]                       then writeln('Carbonylverbindung'); */
  if (fg[fg_aldehyde - 1])
    printf ("Aldehyd\n");
  if (fg[fg_ketone - 1])
    printf ("Keton\n");
  /*  if fg[fg_thiocarbonyl]                   then writeln('Thiocarbonylverbindung'); */
  if (fg[fg_thioaldehyde - 1])
    printf ("Thioaldehyd\n");
  if (fg[fg_thioketone - 1])
    printf ("Thioketon\n");
  if (fg[fg_imine - 1])
    printf ("Imin\n");
  if (fg[fg_hydrazone - 1])
    printf ("Hydrazon\n");
  if (fg[fg_semicarbazone - 1])
    printf ("Semicarbazon\n");
  if (fg[fg_thiosemicarbazone - 1])
    printf ("Thiosemicarbazon\n");
  if (fg[fg_oxime - 1])
    printf ("Oxim\n");
  if (fg[fg_oxime_ether - 1])
    printf ("Oximether\n");
  if (fg[fg_ketene - 1])
    printf ("Keten\n");
  if (fg[fg_ketene_acetal_deriv - 1])
    printf ("Keten-Acetal oder Derivat\n");
  if (fg[fg_carbonyl_hydrate - 1])
    printf ("Carbonyl-Hydrat\n");
  if (fg[fg_hemiacetal - 1])
    printf ("Halbacetal\n");
  if (fg[fg_acetal - 1])
    printf ("Acetal\n");
  if (fg[fg_hemiaminal - 1])
    printf ("Halbaminal\n");
  if (fg[fg_aminal - 1])
    printf ("Aminal\n");
  if (fg[fg_thiohemiaminal - 1])
    printf ("Thiohalbaminal\n");
  if (fg[fg_thioacetal - 1])
    printf ("Thioacetal\n");
  if (fg[fg_enamine - 1])
    printf ("Enamin\n");
  if (fg[fg_enol - 1])
    printf ("Enol\n");
  if (fg[fg_enolether - 1])
    printf ("Enolether\n");
  if (fg[fg_hydroxy - 1] && hydroxy_generic)
    printf ("Hydroxy-Verbindung\n");
  /*  if fg[fg_alcohol]                        then writeln('Alkohol'); */
  if (fg[fg_prim_alcohol - 1])
    printf ("prim\344rer Alkohol\n");
/* p2c: checkmol.pas, line 7283: Note: character >= 128 encountered [281] */
  if (fg[fg_sec_alcohol - 1])
    printf ("sekund\344rer Alkohol\n");
/* p2c: checkmol.pas, line 7284: Note: character >= 128 encountered [281] */
  if (fg[fg_tert_alcohol - 1])
    printf ("terti\344rer Alkohol\n");
/* p2c: checkmol.pas, line 7285: Note: character >= 128 encountered [281] */
  if (fg[fg_1_2_diol - 1])
    printf ("1,2-Diol\n");
  if (fg[fg_1_2_aminoalcohol - 1])
    printf ("1,2-Aminoalkohol\n");
  if (fg[fg_phenol - 1])
    printf ("Phenol oder Hydroxyhetaren\n");
  if (fg[fg_1_2_diphenol - 1])
    printf ("1,2-Diphenol\n");
  if (fg[fg_enediol - 1])
    printf ("Endiol\n");
  if (fg[fg_ether - 1] && ether_generic)
    printf ("Ether\n");
  if (fg[fg_dialkylether - 1])
    printf ("Dialkylether\n");
  if (fg[fg_alkylarylether - 1])
    printf ("Alkylarylether \n");
  if (fg[fg_diarylether - 1])
    printf ("Diarylether\n");
  if (fg[fg_thioether - 1])
    printf ("Thioether\n");
  if (fg[fg_disulfide - 1])
    printf ("Disulfid\n");
  if (fg[fg_peroxide - 1])
    printf ("Peroxid\n");
  if (fg[fg_hydroperoxide - 1])
    printf ("Hydroperoxid\n");
  if (fg[fg_hydrazine - 1])
    printf ("Hydrazin-Derivat\n");
  if (fg[fg_hydroxylamine - 1])
    printf ("Hydroxylamin\n");
  if (fg[fg_amine - 1] && amine_generic)
    printf ("Amin\n");
  if (fg[fg_prim_amine - 1])
    printf ("prim\344res Amin\n");
/* p2c: checkmol.pas, line 7302: Note: character >= 128 encountered [281] */
  if (fg[fg_prim_aliph_amine - 1])
    printf ("prim\344res aliphatisches Amin (Alkylamin)\n");
/* p2c: checkmol.pas, line 7303: Note: character >= 128 encountered [281] */
  if (fg[fg_prim_arom_amine - 1])
    printf ("prim\344res aromatisches Amin\n");
/* p2c: checkmol.pas, line 7304: Note: character >= 128 encountered [281] */
  if (fg[fg_sec_amine - 1])
    printf ("sekund\344res Amin\n");
/* p2c: checkmol.pas, line 7305: Note: character >= 128 encountered [281] */
  if (fg[fg_sec_aliph_amine - 1])
    printf ("sekund\344res aliphatisches Amin (Dialkylamin)\n");
/* p2c: checkmol.pas, line 7306: Note: character >= 128 encountered [281] */
  if (fg[fg_sec_mixed_amine - 1])
    printf
      ("sekund\344res aliphatisches/aromatisches Amin (Alkylarylamin)\n");
/* p2c: checkmol.pas, line 7307: Note: character >= 128 encountered [281] */
  if (fg[fg_sec_arom_amine - 1])
    printf ("sekund\344res aromatisches Amin (Diarylamin)\n");
/* p2c: checkmol.pas, line 7308: Note: character >= 128 encountered [281] */
  if (fg[fg_tert_amine - 1])
    printf ("terti\344res Amin\n");
/* p2c: checkmol.pas, line 7309: Note: character >= 128 encountered [281] */
  if (fg[fg_tert_aliph_amine - 1])
    printf ("terti\344res aliphatisches Amin (Trialkylamin)\n");
/* p2c: checkmol.pas, line 7310: Note: character >= 128 encountered [281] */
  if (fg[fg_tert_mixed_amine - 1])
    printf ("terti\344res aliphatisches/aromatisches Amin (Alkylarylamin)\n");
/* p2c: checkmol.pas, line 7311: Note: character >= 128 encountered [281] */
  if (fg[fg_tert_arom_amine - 1])
    printf ("terti\344res aromatisches Amin (Triarylamin)\n");
/* p2c: checkmol.pas, line 7312: Note: character >= 128 encountered [281] */
  if (fg[fg_quart_ammonium - 1])
    printf ("quart\344res Ammoniumsalz\n");
/* p2c: checkmol.pas, line 7313: Note: character >= 128 encountered [281] */
  if (fg[fg_n_oxide - 1])
    printf ("N-Oxid\n");
  /* new in v0.2f */
  if (fg[fg_halogen_deriv - 1])
    {
      if (!fg[fg_alkyl_halide - 1] && !fg[fg_aryl_halide - 1] &&
	  !fg[fg_acyl_halide - 1])
	printf ("Halogenverbindung\n");
    }
  /*  if fg[fg_alkyl_halide]                   then writeln('Alkylhalogenid'); */
  if (fg[fg_alkyl_fluoride - 1])
    printf ("Alkylfluorid\n");
  if (fg[fg_alkyl_chloride - 1])
    printf ("Alkylchlorid\n");
  if (fg[fg_alkyl_bromide - 1])
    printf ("Alkylbromid\n");
  if (fg[fg_alkyl_iodide - 1])
    printf ("Alkyliodid\n");
  /*  if fg[fg_aryl_halide]                    then writeln('Arylhalogenid'); */
  if (fg[fg_aryl_fluoride - 1])
    printf ("Arylfluorid\n");
  if (fg[fg_aryl_chloride - 1])
    printf ("Arylchlorid\n");
  if (fg[fg_aryl_bromide - 1])
    printf ("Arylbromid\n");
  if (fg[fg_aryl_iodide - 1])
    printf ("Aryliodid\n");
  if (fg[fg_organometallic - 1])
    printf ("Organometall-Verbindung\n");
  if (fg[fg_organolithium - 1])
    printf ("Organolithium-Verbindung\n");
  if (fg[fg_organomagnesium - 1])
    printf ("Organomagnesium-Verbindung\n");
  /*  if fg[fg_carboxylic_acid_deriv]          then writeln('Carbonsäure-Derivat'); */
  if (fg[fg_carboxylic_acid - 1])
    printf ("Carbons\344ure\n");
/* p2c: checkmol.pas, line 7335: Note: character >= 128 encountered [281] */
  if (fg[fg_carboxylic_acid_salt - 1])
    printf ("Carbons\344uresalz\n");
/* p2c: checkmol.pas, line 7336: Note: character >= 128 encountered [281] */
  if (fg[fg_carboxylic_acid_ester - 1])
    printf ("Carbons\344ureester\n");
/* p2c: checkmol.pas, line 7337: Note: character >= 128 encountered [281] */
  if (fg[fg_lactone - 1])
    printf ("Lacton\n");
  /*  if fg[fg_carboxylic_acid_amide]          then writeln('Carbonsäureamid'); */
  if (fg[fg_carboxylic_acid_prim_amide - 1])
    printf ("prim\344res Carbons\344ureamid\n");
/* p2c: checkmol.pas, line 7340:
 * Note: characters >= 128 encountered [281] */
  if (fg[fg_carboxylic_acid_sec_amide - 1])
    printf ("sekund\344res Carbons\344ureamid\n");
/* p2c: checkmol.pas, line 7341:
 * Note: characters >= 128 encountered [281] */
  if (fg[fg_carboxylic_acid_tert_amide - 1])
    printf ("terti\344res Carbons\344ureamid\n");
/* p2c: checkmol.pas, line 7342:
 * Note: characters >= 128 encountered [281] */
  if (fg[fg_lactam - 1])
    printf ("Lactam\n");
  if (fg[fg_carboxylic_acid_hydrazide - 1])
    printf ("Carbons\344urehydrazid\n");
/* p2c: checkmol.pas, line 7344: Note: character >= 128 encountered [281] */
  if (fg[fg_carboxylic_acid_azide - 1])
    printf ("Carbons\344ureazid\n");
/* p2c: checkmol.pas, line 7345: Note: character >= 128 encountered [281] */
  if (fg[fg_hydroxamic_acid - 1])
    printf ("Hydroxams\344ure\n");
/* p2c: checkmol.pas, line 7346: Note: character >= 128 encountered [281] */
  if (fg[fg_carboxylic_acid_amidine - 1])
    printf ("Carbons\344ureamidin\n");
/* p2c: checkmol.pas, line 7347: Note: character >= 128 encountered [281] */
  if (fg[fg_carboxylic_acid_amidrazone - 1])
    printf ("Carbons\344ureamidrazon\n");
/* p2c: checkmol.pas, line 7348: Note: character >= 128 encountered [281] */
  if (fg[fg_nitrile - 1])
    printf ("Carbonitril\n");
  /*  if fg[fg_acyl_halide]                    then writeln('Acylhalogenid'); */
  if (fg[fg_acyl_fluoride - 1])
    printf ("Acylfluorid\n");
  if (fg[fg_acyl_chloride - 1])
    printf ("Acylchlorid\n");
  if (fg[fg_acyl_bromide - 1])
    printf ("Acylbromid\n");
  if (fg[fg_acyl_iodide - 1])
    printf ("Acyliodid\n");
  if (fg[fg_acyl_cyanide - 1])
    printf ("Acylcyanid\n");
  if (fg[fg_imido_ester - 1])
    printf ("Imidoester\n");
  if (fg[fg_imidoyl_halide - 1])
    printf ("Imidoylhalogenid\n");
  /*  if fg[fg_thiocarboxylic_acid_deriv]      then writeln('Thiocarbonsäure-Derivat'); */
  if (fg[fg_thiocarboxylic_acid - 1])
    printf ("Thiocarbons\344ure\n");
/* p2c: checkmol.pas, line 7359: Note: character >= 128 encountered [281] */
  if (fg[fg_thiocarboxylic_acid_ester - 1])
    printf ("Thiocarbons\344ureester\n");
/* p2c: checkmol.pas, line 7360: Note: character >= 128 encountered [281] */
  if (fg[fg_thiolactone - 1])
    printf ("Thiolacton\n");
  if (fg[fg_thiocarboxylic_acid_amide - 1])
    printf ("Thiocarbons\344ureamid\n");
/* p2c: checkmol.pas, line 7362: Note: character >= 128 encountered [281] */
  if (fg[fg_thiolactam - 1])
    printf ("Thiolactam\n");
  if (fg[fg_imido_thioester - 1])
    printf ("Imidothioester\n");
  if (fg[fg_oxohetarene - 1])
    printf ("Oxo(het)aren\n");
  if (fg[fg_thioxohetarene - 1])
    printf ("Thioxo(het)aren\n");
  if (fg[fg_iminohetarene - 1])
    printf ("Imino(het)aren\n");
  if (fg[fg_orthocarboxylic_acid_deriv - 1])
    printf ("Orthocarbons\344ure-Derivat\n");
/* p2c: checkmol.pas, line 7368: Note: character >= 128 encountered [281] */
  if (fg[fg_carboxylic_acid_orthoester - 1])
    printf ("Orthoester\n");
  if (fg[fg_carboxylic_acid_amide_acetal - 1])
    printf ("Amidacetal\n");
  if (fg[fg_carboxylic_acid_anhydride - 1])
    printf ("Carbons\344ureanhydrid\n");
/* p2c: checkmol.pas, line 7371: Note: character >= 128 encountered [281] */
  /*  if fg[fg_carboxylic_acid_imide]          then writeln('Carbonsäureimid'); */
  if (fg[fg_carboxylic_acid_unsubst_imide - 1])
    printf ("Carbons\344ureimid, N-unsubstituiert\n");
/* p2c: checkmol.pas, line 7373: Note: character >= 128 encountered [281] */
  if (fg[fg_carboxylic_acid_subst_imide - 1])
    printf ("Carbons\344ureimid, N-substituiert\n");
/* p2c: checkmol.pas, line 7374: Note: character >= 128 encountered [281] */
  if (fg[fg_co2_deriv - 1])
    printf ("CO2-Derivat (allgemein)\n");
  if (fg[fg_carbonic_acid_deriv - 1] &&
      !(fg[fg_carbonic_acid_monoester - 1]
	|| fg[fg_carbonic_acid_diester - 1]
	|| fg[fg_carbonic_acid_ester_halide - 1]))
    /* changed in v0.3c */
    printf ("Kohlens\344ure-Derivat\n");
/* p2c: checkmol.pas, line 7379: Note: character >= 128 encountered [281] */
  if (fg[fg_carbonic_acid_monoester - 1])
    printf ("Kohlens\344uremonoester\n");
/* p2c: checkmol.pas, line 7380: Note: character >= 128 encountered [281] */
  if (fg[fg_carbonic_acid_diester - 1])
    printf ("Kohlens\344urediester\n");
/* p2c: checkmol.pas, line 7381: Note: character >= 128 encountered [281] */
  if (fg[fg_carbonic_acid_ester_halide - 1])
    printf ("Kohlens\344ureesterhalogenid (Alkyl/Aryl-Halogenformiat)\n");
/* p2c: checkmol.pas, line 7382: Note: character >= 128 encountered [281] */
  if (fg[fg_thiocarbonic_acid_deriv - 1])
    printf ("Thiokohlens\344ure-Derivat\n");
/* p2c: checkmol.pas, line 7383: Note: character >= 128 encountered [281] */
  if (fg[fg_thiocarbonic_acid_monoester - 1])
    printf ("Thiokohlens\344uremonoester\n");
/* p2c: checkmol.pas, line 7384: Note: character >= 128 encountered [281] */
  if (fg[fg_thiocarbonic_acid_diester - 1])
    printf ("Thiokohlens\344urediester\n");
/* p2c: checkmol.pas, line 7385: Note: character >= 128 encountered [281] */
  if (fg[fg_thiocarbonic_acid_ester_halide - 1])
    printf
      ("Thiokohlens\344ureesterhalogenid (Alkyl/Aryl-Halogenthioformiat)\n");
/* p2c: checkmol.pas, line 7386: Note: character >= 128 encountered [281] */
  if (fg[fg_carbamic_acid_deriv - 1] &&
      !(fg[fg_carbamic_acid - 1] || fg[fg_carbamic_acid_ester - 1] ||
	fg[fg_carbamic_acid_halide - 1]))
    /* changed in v0.3c */
    printf ("Carbamins\344ure-Derivat\n");
/* p2c: checkmol.pas, line 7390: Note: character >= 128 encountered [281] */
  if (fg[fg_carbamic_acid - 1])
    printf ("Carbamins\344ure\n");
/* p2c: checkmol.pas, line 7391: Note: character >= 128 encountered [281] */
  if (fg[fg_carbamic_acid_ester - 1])
    printf ("Carbamins\344ureester (Urethan)\n");
/* p2c: checkmol.pas, line 7392: Note: character >= 128 encountered [281] */
  if (fg[fg_carbamic_acid_halide - 1])
    printf ("Carbamins\344urehalogenid (Halogenformamid)\n");
/* p2c: checkmol.pas, line 7393: Note: character >= 128 encountered [281] */
  if (fg[fg_thiocarbamic_acid_deriv - 1] &&
      !(fg[fg_thiocarbamic_acid - 1] || fg[fg_thiocarbamic_acid_ester - 1]
	|| fg[fg_thiocarbamic_acid_halide - 1]))
    /* changed in v0.3c */
    printf ("Thiocarbamins\344ure-Derivat\n");
/* p2c: checkmol.pas, line 7397: Note: character >= 128 encountered [281] */
  if (fg[fg_thiocarbamic_acid - 1])
    printf ("Thiocarbamins\344ure\n");
/* p2c: checkmol.pas, line 7398: Note: character >= 128 encountered [281] */
  if (fg[fg_thiocarbamic_acid_ester - 1])
    printf ("Thiocarbamins\344ureester\n");
/* p2c: checkmol.pas, line 7399: Note: character >= 128 encountered [281] */
  if (fg[fg_thiocarbamic_acid_halide - 1])
    printf ("Thiocarbamins\344urehalogenid (Halogenthioformamid)\n");
/* p2c: checkmol.pas, line 7400: Note: character >= 128 encountered [281] */
  if (fg[fg_urea - 1])
    printf ("Harnstoff\n");
  if (fg[fg_isourea - 1])
    printf ("Isoharnstoff\n");
  if (fg[fg_thiourea - 1])
    printf ("Thioharnstoff\n");
  if (fg[fg_isothiourea - 1])
    printf ("Isothioharnstoff\n");
  if (fg[fg_guanidine - 1])
    printf ("Guanidin\n");
  if (fg[fg_semicarbazide - 1])
    printf ("Semicarbazid\n");
  if (fg[fg_thiosemicarbazide - 1])
    printf ("Thiosemicarbazid\n");
  if (fg[fg_azide - 1])
    printf ("Azid\n");
  if (fg[fg_azo_compound - 1])
    printf ("Azoverbindung\n");
  if (fg[fg_diazonium_salt - 1])
    printf ("Diazoniumsalz\n");
  if (fg[fg_isonitrile - 1])
    printf ("Isonitril\n");
  if (fg[fg_cyanate - 1])
    printf ("Cyanat\n");
  if (fg[fg_isocyanate - 1])
    printf ("Isocyanat\n");
  if (fg[fg_thiocyanate - 1])
    printf ("Thiocyanat\n");
  if (fg[fg_isothiocyanate - 1])
    printf ("Isothiocyanat\n");
  if (fg[fg_carbodiimide - 1])
    printf ("Carbodiimid\n");
  if (fg[fg_nitroso_compound - 1])
    printf ("Nitroso-Verbindung\n");
  if (fg[fg_nitro_compound - 1])
    printf ("Nitro-Verbindung\n");
  if (fg[fg_nitrite - 1])
    printf ("Nitrit\n");
  if (fg[fg_nitrate - 1])
    printf ("Nitrat\n");
  /*  if fg[fg_sulfuric_acid_deriv]            then writeln('Schwefelsäure-Derivat'); */
  if (fg[fg_sulfuric_acid - 1])
    printf ("Schwefels\344ure\n");
/* p2c: checkmol.pas, line 7422: Note: character >= 128 encountered [281] */
  if (fg[fg_sulfuric_acid_monoester - 1])
    printf ("Schwefels\344uremonoester\n");
/* p2c: checkmol.pas, line 7423: Note: character >= 128 encountered [281] */
  if (fg[fg_sulfuric_acid_diester - 1])
    printf ("Schwefels\344urediester\n");
/* p2c: checkmol.pas, line 7424: Note: character >= 128 encountered [281] */
  if (fg[fg_sulfuric_acid_amide_ester - 1])
    printf ("Schwefels\344ureamidester\n");
/* p2c: checkmol.pas, line 7425: Note: character >= 128 encountered [281] */
  if (fg[fg_sulfuric_acid_amide - 1])
    printf ("Schwefels\344ureamid\n");
/* p2c: checkmol.pas, line 7426: Note: character >= 128 encountered [281] */
  if (fg[fg_sulfuric_acid_diamide - 1])
    printf ("Schwefels\344urediamid\n");
/* p2c: checkmol.pas, line 7427: Note: character >= 128 encountered [281] */
  if (fg[fg_sulfuryl_halide - 1])
    printf ("Sulfurylhalogenid\n");
  /*  if fg[fg_sulfonic_acid_deriv]            then writeln('Sulfonsäure-Derivat '); */
  if (fg[fg_sulfonic_acid - 1])
    printf ("Sulfons\344ure\n");
/* p2c: checkmol.pas, line 7430: Note: character >= 128 encountered [281] */
  if (fg[fg_sulfonic_acid_ester - 1])
    printf ("Sulfons\344ureester\n");
/* p2c: checkmol.pas, line 7431: Note: character >= 128 encountered [281] */
  if (fg[fg_sulfonamide - 1])
    printf ("Sulfonamid\n");
  if (fg[fg_sulfonyl_halide - 1])
    printf ("Sulfonylhalogenid\n");
  if (fg[fg_sulfone - 1])
    printf ("Sulfon\n");
  if (fg[fg_sulfoxide - 1])
    printf ("Sulfoxid\n");
  /*  if fg[fg_sulfinic_acid_deriv]            then writeln('Sulfinsäure-Derivat'); */
  if (fg[fg_sulfinic_acid - 1])
    printf ("Sulfins\344ure\n");
/* p2c: checkmol.pas, line 7437: Note: character >= 128 encountered [281] */
  if (fg[fg_sulfinic_acid_ester - 1])
    printf ("Sulfins\344ureester\n");
/* p2c: checkmol.pas, line 7438: Note: character >= 128 encountered [281] */
  if (fg[fg_sulfinic_acid_halide - 1])
    printf ("Sulfins\344urehalogenid\n");
/* p2c: checkmol.pas, line 7439: Note: character >= 128 encountered [281] */
  if (fg[fg_sulfinic_acid_amide - 1])
    printf ("Sulfins\344ureamid\n");
/* p2c: checkmol.pas, line 7440: Note: character >= 128 encountered [281] */
  /*  if fg[fg_sulfenic_acid_deriv]            then writeln('Sulfensäure-Derivat'); */
  if (fg[fg_sulfenic_acid - 1])
    printf ("Sulfens\344ure\n");
/* p2c: checkmol.pas, line 7442: Note: character >= 128 encountered [281] */
  if (fg[fg_sulfenic_acid_ester - 1])
    printf ("Sulfens\344ureester\n");
/* p2c: checkmol.pas, line 7443: Note: character >= 128 encountered [281] */
  if (fg[fg_sulfenic_acid_halide - 1])
    printf ("Sulfens\344urehalogenid\n");
/* p2c: checkmol.pas, line 7444: Note: character >= 128 encountered [281] */
  if (fg[fg_sulfenic_acid_amide - 1])
    printf ("Sulfens\344ureamid\n");
/* p2c: checkmol.pas, line 7445: Note: character >= 128 encountered [281] */
  if (fg[fg_thiol - 1])
    printf ("Thiol (Sulfanyl-Verbindung, Mercaptan)\n");
  if (fg[fg_alkylthiol - 1])
    printf ("Alkylthiol\n");
  if (fg[fg_arylthiol - 1])
    printf ("Arylthiol\n");
  /*  if fg[fg_phosphoric_acid_deriv]          then writeln('Phosphorsäure-Derivat'); */
  if (fg[fg_phosphoric_acid - 1])
    printf ("Phosphors\344ure\n");
/* p2c: checkmol.pas, line 7450: Note: character >= 128 encountered [281] */
  if (fg[fg_phosphoric_acid_ester - 1])
    printf ("Phosphors\344ureester\n");
/* p2c: checkmol.pas, line 7451: Note: character >= 128 encountered [281] */
  if (fg[fg_phosphoric_acid_halide - 1])
    printf ("Phosphors\344urehalogenid\n");
/* p2c: checkmol.pas, line 7452: Note: character >= 128 encountered [281] */
  if (fg[fg_phosphoric_acid_amide - 1])
    printf ("Phosphors\344ureamid\n");
/* p2c: checkmol.pas, line 7453: Note: character >= 128 encountered [281] */
  /*  if fg[fg_thiophosphoric_acid_deriv]      then writeln('Thiophosphorsäure-Derivat'); */
  if (fg[fg_thiophosphoric_acid - 1])
    printf ("Thiophosphors\344ure\n");
/* p2c: checkmol.pas, line 7455: Note: character >= 128 encountered [281] */
  if (fg[fg_thiophosphoric_acid_ester - 1])
    printf ("Thiophosphors\344ureester\n");
/* p2c: checkmol.pas, line 7456: Note: character >= 128 encountered [281] */
  if (fg[fg_thiophosphoric_acid_halide - 1])
    printf ("Thiophosphors\344urehalogenid\n");
/* p2c: checkmol.pas, line 7457: Note: character >= 128 encountered [281] */
  if (fg[fg_thiophosphoric_acid_amide - 1])
    printf ("Thiophosphors\344ureamid\n");
/* p2c: checkmol.pas, line 7458: Note: character >= 128 encountered [281] */
  if (fg[fg_phosphonic_acid_deriv - 1])
    printf ("Phosphons\344ure-Derivat \n");
/* p2c: checkmol.pas, line 7459: Note: character >= 128 encountered [281] */
  if (fg[fg_phosphonic_acid - 1])
    printf ("Phosphons\344ure\n");
/* p2c: checkmol.pas, line 7460: Note: character >= 128 encountered [281] */
  if (fg[fg_phosphonic_acid_ester - 1])
    printf ("Phosphons\344ureester\n");
/* p2c: checkmol.pas, line 7461: Note: character >= 128 encountered [281] */
  if (fg[fg_phosphine - 1])
    printf ("Phosphin\n");
  if (fg[fg_phosphinoxide - 1])
    printf ("Phosphinoxid\n");
  if (fg[fg_boronic_acid_deriv - 1])
    printf ("Borons\344ure-Derivat\n");
/* p2c: checkmol.pas, line 7464: Note: character >= 128 encountered [281] */
  if (fg[fg_boronic_acid - 1])
    printf ("Borons\344ure\n");
/* p2c: checkmol.pas, line 7465: Note: character >= 128 encountered [281] */
  if (fg[fg_boronic_acid_ester - 1])
    printf ("Borons\344ureester\n");
/* p2c: checkmol.pas, line 7466: Note: character >= 128 encountered [281] */
  if (fg[fg_alkene - 1])
    printf ("Alken\n");
  if (fg[fg_alkyne - 1])
    printf ("Alkin\n");
  if (fg[fg_aromatic - 1])
    printf ("aromatische Verbindung\n");
  if (fg[fg_heterocycle - 1])
    printf ("heterocyclische Verbindung\n");
  if (fg[fg_alpha_aminoacid - 1])
    printf ("alpha-Aminos\344ure\n");
/* p2c: checkmol.pas, line 7471: Note: character >= 128 encountered [281] */
  if (fg[fg_alpha_hydroxyacid - 1])
    printf ("alpha-Hydroxys\344ure\n");
/* p2c: checkmol.pas, line 7472: Note: character >= 128 encountered [281] */
}


#define sc              ';'


static void
write_fg_code ()
{
  if (fg[fg_cation - 1])
    printf ("000000T2%c", sc);
  if (fg[fg_anion - 1])
    printf ("000000T1%c", sc);
  /*  if fg[fg_carbonyl]                       then write('C2O10000',sc); */
  if (fg[fg_aldehyde - 1])
    printf ("C2O1H000%c", sc);
  if (fg[fg_ketone - 1])
    printf ("C2O1C000%c", sc);
  /*  if fg[fg_thiocarbonyl]                   then write('C2S10000',sc); */
  if (fg[fg_thioaldehyde - 1])
    printf ("C2S1H000%c", sc);
  if (fg[fg_thioketone - 1])
    printf ("C2S1C000%c", sc);
  if (fg[fg_imine - 1])
    printf ("C2N10000%c", sc);
  if (fg[fg_hydrazone - 1])
    printf ("C2N1N000%c", sc);
  if (fg[fg_semicarbazone - 1])
    printf ("C2NNC4ON%c", sc);
  if (fg[fg_thiosemicarbazone - 1])
    printf ("C2NNC4SN%c", sc);
  if (fg[fg_oxime - 1])
    printf ("C2N1OH00%c", sc);
  if (fg[fg_oxime_ether - 1])
    printf ("C2N1OC00%c", sc);
  if (fg[fg_ketene - 1])
    printf ("C3OC0000%c", sc);
  if (fg[fg_ketene_acetal_deriv - 1])
    printf ("C3OCC000%c", sc);
  if (fg[fg_carbonyl_hydrate - 1])
    printf ("C2O2H200%c", sc);
  if (fg[fg_hemiacetal - 1])
    printf ("C2O2HC00%c", sc);
  if (fg[fg_acetal - 1])
    printf ("C2O2CC00%c", sc);
  if (fg[fg_hemiaminal - 1])
    printf ("C2NOHC10%c", sc);
  if (fg[fg_aminal - 1])
    printf ("C2N2CC10%c", sc);
  if (fg[fg_thiohemiaminal - 1])
    printf ("C2NSHC10%c", sc);
  if (fg[fg_thioacetal - 1])
    printf ("C2S2CC00%c", sc);
  if (fg[fg_enamine - 1])
    printf ("C2CNH000%c", sc);
  if (fg[fg_enol - 1])
    printf ("C2COH000%c", sc);
  if (fg[fg_enolether - 1])
    printf ("C2COC000%c", sc);
  if (fg[fg_hydroxy - 1] && hydroxy_generic)
    printf ("O1H00000%c", sc);
  /*  if fg[fg_alcohol]                        then write('O1H0C000',sc); */
  if (fg[fg_prim_alcohol - 1])
    printf ("O1H1C000%c", sc);
  if (fg[fg_sec_alcohol - 1])
    printf ("O1H2C000%c", sc);
  if (fg[fg_tert_alcohol - 1])
    printf ("O1H3C000%c", sc);
  if (fg[fg_1_2_diol - 1])
    printf ("O1H0CO1H%c", sc);
  if (fg[fg_1_2_aminoalcohol - 1])
    printf ("O1H0CN1C%c", sc);
  if (fg[fg_phenol - 1])
    printf ("O1H1A000%c", sc);
  if (fg[fg_1_2_diphenol - 1])
    printf ("O1H2A000%c", sc);
  if (fg[fg_enediol - 1])
    printf ("C2COH200%c", sc);
  if (fg[fg_ether - 1] && ether_generic)
    printf ("O1C00000%c", sc);
  if (fg[fg_dialkylether - 1])
    printf ("O1C0CC00%c", sc);
  if (fg[fg_alkylarylether - 1])
    printf ("O1C0CA00%c", sc);
  if (fg[fg_diarylether - 1])
    printf ("O1C0AA00%c", sc);
  if (fg[fg_thioether - 1])
    printf ("S1C00000%c", sc);
  if (fg[fg_disulfide - 1])
    printf ("S1S1C000%c", sc);
  if (fg[fg_peroxide - 1])
    printf ("O1O1C000%c", sc);
  if (fg[fg_hydroperoxide - 1])
    printf ("O1O1H000%c", sc);
  if (fg[fg_hydrazine - 1])
    printf ("N1N10000%c", sc);
  if (fg[fg_hydroxylamine - 1])
    printf ("N1O1H000%c", sc);
  if (fg[fg_amine - 1] && amine_generic)
    printf ("N1C00000%c", sc);
  /*  if fg[fg_prim_amine]                     then write('N1C10000',sc); */
  if (fg[fg_prim_aliph_amine - 1])
    printf ("N1C1C000%c", sc);
  if (fg[fg_prim_arom_amine - 1])
    printf ("N1C1A000%c", sc);
  /*  if fg[fg_sec_amine]                      then write('N1C20000',sc); */
  if (fg[fg_sec_aliph_amine - 1])
    printf ("N1C2CC00%c", sc);
  if (fg[fg_sec_mixed_amine - 1])
    printf ("N1C2AC00%c", sc);
  if (fg[fg_sec_arom_amine - 1])
    printf ("N1C2AA00%c", sc);
  /*  if fg[fg_tert_amine]                     then write('N1C30000',sc); */
  if (fg[fg_tert_aliph_amine - 1])
    printf ("N1C3CC00%c", sc);
  if (fg[fg_tert_mixed_amine - 1])
    printf ("N1C3AC00%c", sc);
  if (fg[fg_tert_arom_amine - 1])
    printf ("N1C3AA00%c", sc);
  if (fg[fg_quart_ammonium - 1])
    printf ("N1C400T2%c", sc);
  if (fg[fg_n_oxide - 1])
    printf ("N0O10000%c", sc);
  /*  if fg[fg_halogen_deriv]                  then write('XX000000',sc); */
  /* new in v0.2f */
  if (fg[fg_halogen_deriv - 1])
    {
      if (!fg[fg_alkyl_halide - 1] && !fg[fg_aryl_halide - 1] &&
	  !fg[fg_acyl_halide - 1])
	printf ("XX000000%c", sc);
    }
  /*  if fg[fg_alkyl_halide]                   then write('XX00C000',sc); */
  if (fg[fg_alkyl_fluoride - 1])
    printf ("XF00C000%c", sc);
  if (fg[fg_alkyl_chloride - 1])
    printf ("XC00C000%c", sc);
  if (fg[fg_alkyl_bromide - 1])
    printf ("XB00C000%c", sc);
  if (fg[fg_alkyl_iodide - 1])
    printf ("XI00C000%c", sc);
  /*  if fg[fg_aryl_halide]                    then write('XX00A000',sc); */
  if (fg[fg_aryl_fluoride - 1])
    printf ("XF00A000%c", sc);
  if (fg[fg_aryl_chloride - 1])
    printf ("XC00A000%c", sc);
  if (fg[fg_aryl_bromide - 1])
    printf ("XB00A000%c", sc);
  if (fg[fg_aryl_iodide - 1])
    printf ("XI00A000%c", sc);
  if (fg[fg_organometallic - 1])
    printf ("000000MX%c", sc);
  if (fg[fg_organolithium - 1])
    printf ("000000ML%c", sc);
  if (fg[fg_organomagnesium - 1])
    printf ("000000MM%c", sc);
  /*  if fg[fg_carboxylic_acid_deriv]          then write('C3O20000',sc); */
  if (fg[fg_carboxylic_acid - 1])
    printf ("C3O2H000%c", sc);
  if (fg[fg_carboxylic_acid_salt - 1])
    printf ("C3O200T1%c", sc);
  if (fg[fg_carboxylic_acid_ester - 1])
    printf ("C3O2C000%c", sc);
  if (fg[fg_lactone - 1])
    printf ("C3O2CZ00%c", sc);
  /*  if fg[fg_carboxylic_acid_amide]          then write('C3ONC000',sc); */
  if (fg[fg_carboxylic_acid_prim_amide - 1])
    printf ("C3ONC100%c", sc);
  if (fg[fg_carboxylic_acid_sec_amide - 1])
    printf ("C3ONC200%c", sc);
  if (fg[fg_carboxylic_acid_tert_amide - 1])
    printf ("C3ONC300%c", sc);
  if (fg[fg_lactam - 1])
    printf ("C3ONCZ00%c", sc);
  if (fg[fg_carboxylic_acid_hydrazide - 1])
    printf ("C3ONN100%c", sc);
  if (fg[fg_carboxylic_acid_azide - 1])
    printf ("C3ONN200%c", sc);
  if (fg[fg_hydroxamic_acid - 1])
    printf ("C3ONOH00%c", sc);
  if (fg[fg_carboxylic_acid_amidine - 1])
    printf ("C3N2H000%c", sc);
  if (fg[fg_carboxylic_acid_amidrazone - 1])
    printf ("C3NNN100%c", sc);
  if (fg[fg_nitrile - 1])
    printf ("C3N00000%c", sc);
  /*  if fg[fg_acyl_halide]                    then write('C3OXX000',sc); */
  if (fg[fg_acyl_fluoride - 1])
    printf ("C3OXF000%c", sc);
  if (fg[fg_acyl_chloride - 1])
    printf ("C3OXC000%c", sc);
  if (fg[fg_acyl_bromide - 1])
    printf ("C3OXB000%c", sc);
  if (fg[fg_acyl_iodide - 1])
    printf ("C3OXI000%c", sc);
  if (fg[fg_acyl_cyanide - 1])
    printf ("C2OC3N00%c", sc);
  if (fg[fg_imido_ester - 1])
    printf ("C3NOC000%c", sc);
  if (fg[fg_imidoyl_halide - 1])
    printf ("C3NXX000%c", sc);
  /*  if fg[fg_thiocarboxylic_acid_deriv]      then write('C3SO0000',sc); */
  if (fg[fg_thiocarboxylic_acid - 1])
    printf ("C3SOH000%c", sc);
  if (fg[fg_thiocarboxylic_acid_ester - 1])
    printf ("C3SOC000%c", sc);
  if (fg[fg_thiolactone - 1])
    printf ("C3SOCZ00%c", sc);
  if (fg[fg_thiocarboxylic_acid_amide - 1])
    printf ("C3SNH000%c", sc);
  if (fg[fg_thiolactam - 1])
    printf ("C3SNCZ00%c", sc);
  if (fg[fg_imido_thioester - 1])
    printf ("C3NSC000%c", sc);
  if (fg[fg_oxohetarene - 1])
    printf ("C3ONAZ00%c", sc);
  if (fg[fg_thioxohetarene - 1])
    printf ("C3SNAZ00%c", sc);
  if (fg[fg_iminohetarene - 1])
    printf ("C3NNAZ00%c", sc);
  if (fg[fg_orthocarboxylic_acid_deriv - 1])
    printf ("C3O30000%c", sc);
  if (fg[fg_carboxylic_acid_orthoester - 1])
    printf ("C3O3C000%c", sc);
  if (fg[fg_carboxylic_acid_amide_acetal - 1])
    printf ("C3O3NC00%c", sc);
  if (fg[fg_carboxylic_acid_anhydride - 1])
    printf ("C3O2C3O2%c", sc);
  /*  if fg[fg_carboxylic_acid_imide]          then write('C3ONC000',sc); */
  if (fg[fg_carboxylic_acid_unsubst_imide - 1])
    printf ("C3ONCH10%c", sc);
  if (fg[fg_carboxylic_acid_subst_imide - 1])
    printf ("C3ONCC10%c", sc);
  if (fg[fg_co2_deriv - 1])
    printf ("C4000000%c", sc);
  if (fg[fg_carbonic_acid_deriv - 1])
    printf ("C4O30000%c", sc);
  if (fg[fg_carbonic_acid_monoester - 1])
    printf ("C4O3C100%c", sc);
  if (fg[fg_carbonic_acid_diester - 1])
    printf ("C4O3C200%c", sc);
  if (fg[fg_carbonic_acid_ester_halide - 1])
    printf ("C4O3CX00%c", sc);
  if (fg[fg_thiocarbonic_acid_deriv - 1])
    printf ("C4SO0000%c", sc);
  if (fg[fg_thiocarbonic_acid_monoester - 1])
    printf ("C4SOC100%c", sc);
  if (fg[fg_thiocarbonic_acid_diester - 1])
    printf ("C4SOC200%c", sc);
  if (fg[fg_thiocarbonic_acid_ester_halide - 1])
    printf ("C4SOX_00%c", sc);
  if (fg[fg_carbamic_acid_deriv - 1])
    printf ("C4O2N000%c", sc);
  if (fg[fg_carbamic_acid - 1])
    printf ("C4O2NH00%c", sc);
  if (fg[fg_carbamic_acid_ester - 1])
    printf ("C4O2NC00%c", sc);
  if (fg[fg_carbamic_acid_halide - 1])
    printf ("C4O2NX00%c", sc);
  if (fg[fg_thiocarbamic_acid_deriv - 1])
    printf ("C4SN0000%c", sc);
  if (fg[fg_thiocarbamic_acid - 1])
    printf ("C4SNOH00%c", sc);
  if (fg[fg_thiocarbamic_acid_ester - 1])
    printf ("C4SNOC00%c", sc);
  if (fg[fg_thiocarbamic_acid_halide - 1])
    printf ("C4SNXX00%c", sc);
  if (fg[fg_urea - 1])
    printf ("C4O1N200%c", sc);
  if (fg[fg_isourea - 1])
    printf ("C4N2O100%c", sc);
  if (fg[fg_thiourea - 1])
    printf ("C4S1N200%c", sc);
  if (fg[fg_isothiourea - 1])
    printf ("C4N2S100%c", sc);
  if (fg[fg_guanidine - 1])
    printf ("C4N30000%c", sc);
  if (fg[fg_semicarbazide - 1])
    printf ("C4ON2N00%c", sc);
  if (fg[fg_thiosemicarbazide - 1])
    printf ("C4SN2N00%c", sc);
  if (fg[fg_azide - 1])
    printf ("N4N20000%c", sc);
  if (fg[fg_azo_compound - 1])
    printf ("N2N10000%c", sc);
  if (fg[fg_diazonium_salt - 1])
    printf ("N3N100T2%c", sc);
  if (fg[fg_isonitrile - 1])
    printf ("N3C10000%c", sc);
  if (fg[fg_cyanate - 1])
    printf ("C4NO1000%c", sc);
  if (fg[fg_isocyanate - 1])
    printf ("C4NO2000%c", sc);
  if (fg[fg_thiocyanate - 1])
    printf ("C4NS1000%c", sc);
  if (fg[fg_isothiocyanate - 1])
    printf ("C4NS2000%c", sc);
  if (fg[fg_carbodiimide - 1])
    printf ("C4N20000%c", sc);
  if (fg[fg_nitroso_compound - 1])
    printf ("N2O10000%c", sc);
  if (fg[fg_nitro_compound - 1])
    printf ("N4O20000%c", sc);
  if (fg[fg_nitrite - 1])
    printf ("N3O20000%c", sc);
  if (fg[fg_nitrate - 1])
    printf ("N4O30000%c", sc);
  if (fg[fg_sulfuric_acid_deriv - 1])
    printf ("S6O00000%c", sc);
  if (fg[fg_sulfuric_acid - 1])
    printf ("S6O4H000%c", sc);
  if (fg[fg_sulfuric_acid_monoester - 1])
    printf ("S6O4HC00%c", sc);
  if (fg[fg_sulfuric_acid_diester - 1])
    printf ("S6O4CC00%c", sc);
  if (fg[fg_sulfuric_acid_amide_ester - 1])
    printf ("S6O3NC00%c", sc);
  if (fg[fg_sulfuric_acid_amide - 1])
    printf ("S6O3N100%c", sc);
  if (fg[fg_sulfuric_acid_diamide - 1])
    printf ("S6O2N200%c", sc);
  if (fg[fg_sulfuryl_halide - 1])
    printf ("S6O3XX00%c", sc);
  if (fg[fg_sulfonic_acid_deriv - 1])
    printf ("S5O00000%c", sc);
  if (fg[fg_sulfonic_acid - 1])
    printf ("S5O3H000%c", sc);
  if (fg[fg_sulfonic_acid_ester - 1])
    printf ("S5O3C000%c", sc);
  if (fg[fg_sulfonamide - 1])
    printf ("S5O2N000%c", sc);
  if (fg[fg_sulfonyl_halide - 1])
    printf ("S5O2XX00%c", sc);
  if (fg[fg_sulfone - 1])
    printf ("S4O20000%c", sc);
  if (fg[fg_sulfoxide - 1])
    printf ("S2O10000%c", sc);
  if (fg[fg_sulfinic_acid_deriv - 1])
    printf ("S3O00000%c", sc);
  if (fg[fg_sulfinic_acid - 1])
    printf ("S3O2H000%c", sc);
  if (fg[fg_sulfinic_acid_ester - 1])
    printf ("S3O2C000%c", sc);
  if (fg[fg_sulfinic_acid_halide - 1])
    printf ("S3O1XX00%c", sc);
  if (fg[fg_sulfinic_acid_amide - 1])
    printf ("S3O1N000%c", sc);
  if (fg[fg_sulfenic_acid_deriv - 1])
    printf ("S1O00000%c", sc);
  if (fg[fg_sulfenic_acid - 1])
    printf ("S1O1H000%c", sc);
  if (fg[fg_sulfenic_acid_ester - 1])
    printf ("S1O1C000%c", sc);
  if (fg[fg_sulfenic_acid_halide - 1])
    printf ("S1O0XX00%c", sc);
  if (fg[fg_sulfenic_acid_amide - 1])
    printf ("S1O0N100%c", sc);
  /*  if fg[fg_thiol]                          then write('S1H10000',sc); */
  if (fg[fg_alkylthiol - 1])
    printf ("S1H1C000%c", sc);
  if (fg[fg_arylthiol - 1])
    printf ("S1H1A000%c", sc);
  if (fg[fg_phosphoric_acid_deriv - 1])
    printf ("P5O0H000%c", sc);
  if (fg[fg_phosphoric_acid - 1])
    printf ("P5O4H200%c", sc);
  if (fg[fg_phosphoric_acid_ester - 1])
    printf ("P5O4HC00%c", sc);
  if (fg[fg_phosphoric_acid_halide - 1])
    printf ("P5O3HX00%c", sc);
  if (fg[fg_phosphoric_acid_amide - 1])
    printf ("P5O3HN00%c", sc);
  if (fg[fg_thiophosphoric_acid_deriv - 1])
    printf ("P5O0S000%c", sc);
  if (fg[fg_thiophosphoric_acid - 1])
    printf ("P5O3SH00%c", sc);
  if (fg[fg_thiophosphoric_acid_ester - 1])
    printf ("P5O3SC00%c", sc);
  if (fg[fg_thiophosphoric_acid_halide - 1])
    printf ("P5O2SX00%c", sc);
  if (fg[fg_thiophosphoric_acid_amide - 1])
    printf ("P5O2SN00%c", sc);
  if (fg[fg_phosphonic_acid_deriv - 1])
    printf ("P4O30000%c", sc);
  if (fg[fg_phosphonic_acid - 1])
    printf ("P4O3H000%c", sc);
  if (fg[fg_phosphonic_acid_ester - 1])
    printf ("P4O3C000%c", sc);
  if (fg[fg_phosphine - 1])
    printf ("P3000000%c", sc);
  if (fg[fg_phosphinoxide - 1])
    printf ("P2O00000%c", sc);
  if (fg[fg_boronic_acid_deriv - 1])
    printf ("B2O20000%c", sc);
  if (fg[fg_boronic_acid - 1])
    printf ("B2O2H000%c", sc);
  if (fg[fg_boronic_acid_ester - 1])
    printf ("B2O2C000%c", sc);
  if (fg[fg_alkene - 1])
    printf ("000C2C00%c", sc);
  if (fg[fg_alkyne - 1])
    printf ("000C3C00%c", sc);
  if (fg[fg_aromatic - 1])
    printf ("0000A000%c", sc);
  if (fg[fg_heterocycle - 1])
    printf ("0000CZ00%c", sc);
  if (fg[fg_alpha_aminoacid - 1])
    printf ("C3O2HN1C%c", sc);
  if (fg[fg_alpha_hydroxyacid - 1])
    printf ("C3O2HO1H%c", sc);
}

#undef sc


static void
write_fg_binary ()
{
  int i, n;
  char o;

  for (i = 1; i <= max_fg / 8; i++)
    {
      n = 0;
      if (fg[i * 8 - 1])
	n++;
      if (fg[i * 8 - 2])
	n += 2;
      if (fg[i * 8 - 3])
	n += 4;
      if (fg[i * 8 - 4])
	n += 8;
      if (fg[i * 8 - 5])
	n += 16;
      if (fg[i * 8 - 6])
	n += 32;
      if (fg[i * 8 - 7])
	n += 64;
      if (fg[i * 8 - 8])
	n += 128;
      o = (char) n;
      putchar (o);
    }
}


static void
write_fg_bitstring ()
{
  int i;

  for (i = 0; i < max_fg; i++)
    {
      if (fg[i])
	putchar ('1');
      else
	putchar ('0');
    }
}


  /*static void readinputfile (molfilename) char *molfilename;
     {
     /* new version in v0.2g 
     char rline[256];
     char *TEMP;

     molbufindex = 0;
     if (!opt_stdin)
     {
     if (!rfile_is_open)
     {
     assign (rfile, molfilename);
     rewind (rfile);
     rfile_is_open = true;
     }
     /* p2c: checkmol.pas, line 7733: Warning:
     * Don't know how to ASSIGN to a non-explicit file variable [207] 
     *rline = '\0';
     mol_in_queue = false;
     while ((!P_eof (rfile)) && (strpos2 (rline, "$$$$", 1) == 0))
     {
     fgets (rline, 256, rfile);
     TEMP = strchr (rline, '\n');
     if (TEMP != NULL)
     *TEMP = 0;
     /*mol_in_queue := false; 
     if (molbufindex >= max_atoms + max_bonds + 64)
     {
     printf ("Not enough memory for molfile! %12ld\n",
     molbufindex);
     if (rfile != NULL)
     fclose (rfile);
     rfile = NULL;
     _Escape (1);
     }
     molbufindex++;
     strcpy (molbuf[molbufindex - 1], rline);
     if (strpos2 (rline, "$$$$", 1) > 0)
     mol_in_queue = true;
     }
     if (!P_eof (rfile))
     return;
     if (rfile != NULL)
     fclose (rfile);
     rfile = NULL;
     rfile_is_open = false;
     mol_in_queue = false;
     return;
     }
     *rline = '\0';
     mol_in_queue = false;
     while ((!P_eof (stdin)) && (strpos2 (rline, "$$$$", 1) == 0))
     {
     gets (rline);
     if (molbufindex >= max_atoms + max_bonds + 64)
     {
     printf ("Not enough memory!\n");
     _Escape (1);
     }
     molbufindex++;
     strcpy (molbuf[molbufindex - 1], rline);
     if (strpos2 (rline, "$$$$", 1) > 0)
     {
     mol_in_queue = true;
     /* read from standard input 
     }
     }
     } */

static void
readinputfile (char *molfilename)
{
  /* new version in v0.2g */
  char rline[256];
  char *TEMP;

  molbufindex = 0;
  if (!opt_stdin)
    {
      if (!rfile_is_open)
	{
	  rfile = fopen (molfilename, "r");
	  rewind (rfile);
	  rfile_is_open = true;
	}
/* p2c: checkmol.pas, line 7226: Warning:
 * Don't know how to ASSIGN to a non-explicit file variable [207] */
      *rline = '\0';
      mol_in_queue = false;
      while ((fgets (rline, 256, rfile) != NULL)
	     && (strstr (rline, "$$$$") == NULL))
	{
	  TEMP = strchr (rline, '\n');
	  if (TEMP != NULL)
	    *TEMP = 0;
	  /*mol_in_queue := false; */
	  if (molbufindex >= max_atoms + max_bonds + 64)
	    {
	      printf ("Not enough memory for molfile! %d\n", molbufindex);
	      if (rfile != NULL)
		fclose (rfile);
	      rfile = NULL;
	      exit (1);
	    }
	  //molbufindex++;
	  strcpy (molbuf[molbufindex++], rline);
	  if (strstr (rline, "$$$$") != NULL)
	    mol_in_queue = true;
	}
      if (!feof (rfile))
	return;
      if (rfile != NULL)
	fclose (rfile);
      rfile = NULL;
      rfile_is_open = false;
      mol_in_queue = false;
      return;
    }
  *rline = '\0';
  mol_in_queue = false;
  do
    {
      fgets (rline, 256, stdin);
      if (feof (stdin))
	return;
      TEMP = strchr (rline, '\n');
      if (TEMP != NULL)
	*TEMP = '\0';
      if (molbufindex >= max_atoms + max_bonds + 64)
	{
	  printf ("Not enough memory!\n");
	  exit (1);
	}
      //molbufindex++;
      strcpy (molbuf[molbufindex++], rline);
      if (strstr (rline, "$$$$") != NULL)
	{
	  mol_in_queue = true;
	  /* read from standard input */
	}
    }
  while (strstr (rline, "$$$$") == NULL);
}



/* static void copy_mol_to_needle()
{
  int i, j, FORLIM;

  if (n_atoms == 0)
    return;
  /* try 
  ndl_atom = (atom_rec *)safe_malloc(n_atoms * sizeof(atom_rec));
  ndl_bond = (bond_rec *)safe_malloc(n_bonds * sizeof(bond_rec));
  ndl_ring = (ringpath_type *)safe_malloc(sizeof(ringlist));
  ndl_ringprop = (ringprop_rec *)safe_malloc(sizeof(ringprop_type));
  /* except
    on e:Eoutofmemory do
      begin
        writeln('Not enough memory');
        halt(4);
      end;
  end; 
  ndl_n_atoms = n_atoms;
  ndl_n_bonds = n_bonds;
  ndl_n_rings = n_rings;
  ndl_n_heavyatoms = n_heavyatoms;
  ndl_n_heavybonds = n_heavybonds;
  strcpy(ndl_molname, molname);
  ndl_n_Ctot = n_Ctot;
  ndl_n_Otot = n_Otot;
  ndl_n_Ntot = n_Ntot;
  FORLIM = n_atoms;
  for (i = 0; i < FORLIM; i++) {
    strcpy(ndl_atom[i].element, atom[i].element);
    strcpy(ndl_atom[i].atype, atom[i].atype);
    ndl_atom[i].x = atom[i].x;
    ndl_atom[i].y = atom[i].y;
    ndl_atom[i].z = atom[i].z;
    ndl_atom[i].formal_charge = atom[i].formal_charge;
    ndl_atom[i].real_charge = atom[i].real_charge;
    ndl_atom[i].Hexp = atom[i].Hexp;
    ndl_atom[i].Htot = atom[i].Htot;
    ndl_atom[i].neighbor_count = atom[i].neighbor_count;
    ndl_atom[i].ring_count = atom[i].ring_count;
    ndl_atom[i].arom = atom[i].arom;
    ndl_atom[i].stereo_care = atom[i].stereo_care;
    ndl_atom[i].heavy = atom[i].heavy;   /* v0.3l 
    ndl_atom[i].metal = atom[i].metal;   /* v0.3l 
  ndl_atom[i].tag = atom[i].tag;	/* v0.3o 
					   }
					   if (n_bonds > 0) {
					   FORLIM = n_bonds;
					   for (i = 0; i < FORLIM; i++) {
					   ndl_bond[i].a1 = bond[i].a1;
					   ndl_bond[i].a2 = bond[i].a2;
					   ndl_bond[i].btype = bond[i].btype;
					   ndl_bond[i].arom = bond[i].arom;
					   ndl_bond[i].ring_count = bond[i].ring_count;   /* new in v0.3d 
					   ndl_bond[i].topo = bond[i].topo;   /* new in v0.3d 
					   ndl_bond[i].stereo = bond[i].stereo;   /* new in v0.3d 
					   }
					   }
					   if (n_rings > 0) {
					   FORLIM = n_rings;
					   for (i = 0; i < FORLIM; i++) {
					   for (j = 0; j < max_ringsize; j++)
					   ndl_ring[i][j] = ring[i][j];
					   }
					   for (i = 0; i < max_rings; i++) {   /* new in v0.3 
					   ndl_ringprop[i].size = ringprop[i].size;
					   ndl_ringprop[i].arom = ringprop[i].arom;
					   ndl_ringprop[i].envelope = ringprop[i].envelope;
					   }
					   }
					   ndl_molstat.n_QA = molstat.n_QA;
					   ndl_molstat.n_QB = molstat.n_QB;
					   ndl_molstat.n_chg = molstat.n_chg;
					   ndl_molstat.n_C1 = molstat.n_C1;
					   ndl_molstat.n_C2 = molstat.n_C2;
					   ndl_molstat.n_C = molstat.n_C;
					   ndl_molstat.n_CHB1p = molstat.n_CHB1p;
					   ndl_molstat.n_CHB2p = molstat.n_CHB2p;
					   ndl_molstat.n_CHB3p = molstat.n_CHB3p;
					   ndl_molstat.n_CHB4 = molstat.n_CHB4;
					   ndl_molstat.n_O2 = molstat.n_O2;
					   ndl_molstat.n_O3 = molstat.n_O3;
					   ndl_molstat.n_N1 = molstat.n_N1;
					   ndl_molstat.n_N2 = molstat.n_N2;
					   ndl_molstat.n_N3 = molstat.n_N3;
					   ndl_molstat.n_S = molstat.n_S;
					   ndl_molstat.n_SeTe = molstat.n_SeTe;
					   ndl_molstat.n_F = molstat.n_F;
					   ndl_molstat.n_Cl = molstat.n_Cl;
					   ndl_molstat.n_Br = molstat.n_Br;
					   ndl_molstat.n_I = molstat.n_I;
					   ndl_molstat.n_P = molstat.n_P;
					   ndl_molstat.n_B = molstat.n_B;
					   ndl_molstat.n_Met = molstat.n_Met;
					   ndl_molstat.n_X = molstat.n_X;
					   ndl_molstat.n_b1 = molstat.n_b1;
					   ndl_molstat.n_b2 = molstat.n_b2;
					   ndl_molstat.n_b3 = molstat.n_b3;
					   ndl_molstat.n_bar = molstat.n_bar;
					   ndl_molstat.n_C1O = molstat.n_C1O;
					   ndl_molstat.n_C2O = molstat.n_C2O;
					   ndl_molstat.n_CN = molstat.n_CN;
					   ndl_molstat.n_XY = molstat.n_XY;
					   ndl_molstat.n_r3 = molstat.n_r3;
					   ndl_molstat.n_r4 = molstat.n_r4;
					   ndl_molstat.n_r5 = molstat.n_r5;
					   ndl_molstat.n_r6 = molstat.n_r6;
					   ndl_molstat.n_r7 = molstat.n_r7;
					   ndl_molstat.n_r8 = molstat.n_r8;
					   ndl_molstat.n_r9 = molstat.n_r9;
					   ndl_molstat.n_r10 = molstat.n_r10;
					   ndl_molstat.n_r11 = molstat.n_r11;
					   ndl_molstat.n_r12 = molstat.n_r12;
					   ndl_molstat.n_r13p = molstat.n_r13p;
					   ndl_molstat.n_rN = molstat.n_rN;
					   ndl_molstat.n_rN1 = molstat.n_rN1;
					   ndl_molstat.n_rN2 = molstat.n_rN2;
					   ndl_molstat.n_rN3p = molstat.n_rN3p;
					   ndl_molstat.n_rO = molstat.n_rO;
					   ndl_molstat.n_rO1 = molstat.n_rO1;
					   ndl_molstat.n_rO2p = molstat.n_rO2p;
					   ndl_molstat.n_rS = molstat.n_rS;
					   ndl_molstat.n_rX = molstat.n_rX;
					   ndl_molstat.n_rAr = molstat.n_rAr;
					   ndl_molstat.n_rBz = molstat.n_rBz;   /* v0.3l 
					   ndl_molstat.n_br2p = molstat.n_br2p;   /* v0.3n 
					   /* p2c: checkmol.pas, line 7875:
					   * Note: Turbo Pascal conditional compilation directive was ignored [218] 
					   /*$IFDEF extended_molstat
					   /* v0.3m 
					   ndl_molstat.n_psg01 = molstat.n_psg01;
					   ndl_molstat.n_psg02 = molstat.n_psg02;
					   ndl_molstat.n_psg13 = molstat.n_psg13;
					   ndl_molstat.n_psg14 = molstat.n_psg14;
					   ndl_molstat.n_psg15 = molstat.n_psg15;
					   ndl_molstat.n_psg16 = molstat.n_psg16;
					   ndl_molstat.n_psg17 = molstat.n_psg17;
					   ndl_molstat.n_psg18 = molstat.n_psg18;
					   ndl_molstat.n_pstm = molstat.n_pstm;
					   ndl_molstat.n_psla = molstat.n_psla;
					   /*$ENDIF
					   /* make sure some modes can be switched on only by the query file 
					   /* and not by subsequent haystack file(s) 
					   if (ez_flag)   /* new in v0.3f 
					   ez_search = true;
					   if (chir_flag)   /* new in v0.3f 
					   rs_search = true;
					   } */


static void
copy_mol_to_needle ()
{
  //int i, j, FORLIM;

  /*if (n_atoms == 0)
     return; *///If a NoStruct is read, this leads to madness and illegal memory access


  ndl_atom = (atom_rec *) safe_calloc (n_atoms, sizeof (atom_rec));
  ndl_bond = (bond_rec *) safe_calloc (n_bonds, sizeof (bond_rec));
  ndl_ring = (ringpath_type *) safe_calloc (1, sizeof (ringlist));
  ndl_ringprop = (ringprop_rec *) safe_calloc (1, sizeof (ringprop_type));


  ndl_n_atoms = n_atoms;
  ndl_n_bonds = n_bonds;
  ndl_n_rings = n_rings;
  ndl_n_heavyatoms = n_heavyatoms;
  ndl_n_heavybonds = n_heavybonds;
  strcpy (ndl_molname, molname);
  ndl_n_Ctot = n_Ctot;
  ndl_n_Otot = n_Otot;
  ndl_n_Ntot = n_Ntot;
  memcpy (ndl_atom, atom, n_atoms * sizeof (atom_rec));

  if (n_bonds > 0)
    memcpy (ndl_bond, bond, n_bonds * sizeof (bond_rec));

  if (n_rings > 0)
    {
      memcpy (ndl_ring, ring, sizeof (ringlist));
      memcpy (ndl_ringprop, ringprop, sizeof (ringprop_type));
    }

  memcpy (&ndl_molstat, &molstat, sizeof (molstat));


  // make sure some modes can be switched on only by the query file 
  // and not by subsequent haystack file(s) 
  if (ez_flag)			// new in v0.3f 
    ez_search = true;

  if (chir_flag)		// new in v0.3f 
    rs_search = true;

  ndl_querymol = found_querymol;	/* 0.3p */

}

static void
copy_mol_to_tmp ()
{
  //int i, j, FORLIM;

  /*if (n_atoms == 0)
     return; *///If a NoStruct is read, this leads to madness and illegal memory access


  tmp_atom = (atom_rec *) safe_calloc (n_atoms, sizeof (atom_rec));
  tmp_bond = (bond_rec *) safe_calloc (n_bonds, sizeof (bond_rec));
  tmp_ring = (ringpath_type *) safe_calloc (1, sizeof (ringlist));
  tmp_ringprop = (ringprop_rec *) safe_calloc (1, sizeof (ringprop_type));


  tmp_n_atoms = n_atoms;
  tmp_n_bonds = n_bonds;
  tmp_n_rings = n_rings;
  tmp_n_heavyatoms = n_heavyatoms;
  tmp_n_heavybonds = n_heavybonds;
  strcpy (tmp_molname, molname);
  tmp_n_Ctot = n_Ctot;
  tmp_n_Otot = n_Otot;
  tmp_n_Ntot = n_Ntot;
  memcpy (tmp_atom, atom, n_atoms * sizeof (atom_rec));

  if (n_bonds > 0)
    memcpy (tmp_bond, bond, n_bonds * sizeof (bond_rec));

  if (n_rings > 0)
    {
      memcpy (tmp_ring, ring, sizeof (ringlist));
      memcpy (tmp_ringprop, ringprop, sizeof (ringprop_type));
    }

  memcpy (&tmp_molstat, &molstat, sizeof (molstat));


  // make sure some modes can be switched on only by the query file 
  // and not by subsequent haystack file(s) 
  if (ez_flag)			// new in v0.3f 
    ez_search = true;

  if (chir_flag)		// new in v0.3f 
    rs_search = true;

  ndl_querymol = found_querymol;	/* 0.3p */

}

/* static void copy_mol_to_tmp()
{
  int i, j, FORLIM;

  if (n_atoms == 0)
    return;
  /* try 
  tmp_atom = (atom_rec *)safe_malloc(n_atoms * sizeof(atom_rec));
  tmp_bond = (bond_rec *)safe_malloc(n_bonds * sizeof(bond_rec));
  tmp_ring = (ringpath_type *)safe_malloc(sizeof(ringlist));
  tmp_ringprop = (ringprop_rec *)safe_malloc(sizeof(ringprop_type));
  /* except
    on e:Eoutofmemory do
      begin
        writeln('Not enough memory');
        halt(4);
      end;
  end;
  tmp_n_atoms = n_atoms;
  tmp_n_bonds = n_bonds;
  tmp_n_rings = n_rings;
  tmp_n_heavyatoms = n_heavyatoms;
  tmp_n_heavybonds = n_heavybonds;
  strcpy(tmp_molname, molname);
  tmp_n_Ctot = n_Ctot;
  tmp_n_Otot = n_Otot;
  tmp_n_Ntot = n_Ntot;
  FORLIM = n_atoms;
  for (i = 0; i < FORLIM; i++) {
    strcpy(tmp_atom[i].element, atom[i].element);
    strcpy(tmp_atom[i].atype, atom[i].atype);
    tmp_atom[i].x = atom[i].x;
    tmp_atom[i].y = atom[i].y;
    tmp_atom[i].z = atom[i].z;
    tmp_atom[i].formal_charge = atom[i].formal_charge;
    tmp_atom[i].real_charge = atom[i].real_charge;
    tmp_atom[i].Hexp = atom[i].Hexp;
    tmp_atom[i].Htot = atom[i].Htot;
    tmp_atom[i].neighbor_count = atom[i].neighbor_count;
    tmp_atom[i].ring_count = atom[i].ring_count;
    tmp_atom[i].arom = atom[i].arom;
    tmp_atom[i].stereo_care = atom[i].stereo_care;
    tmp_atom[i].heavy = atom[i].heavy;   /* v0.3l 
    tmp_atom[i].metal = atom[i].metal;   /* v0.3l 
    tmp_atom[i].tag = atom[i].tag;   /* v0.3o 
  }
  if (n_bonds > 0) {
    FORLIM = n_bonds;
    for (i = 0; i < FORLIM; i++) {
      tmp_bond[i].a1 = bond[i].a1;
      tmp_bond[i].a2 = bond[i].a2;
      tmp_bond[i].btype = bond[i].btype;
      tmp_bond[i].arom = bond[i].arom;
      tmp_bond[i].ring_count = bond[i].ring_count;   /* new in v0.3d 
      tmp_bond[i].topo = bond[i].topo;   /* new in v0.3d 
      tmp_bond[i].stereo = bond[i].stereo;   /* new in v0.3d 
    }
  }
  if (n_rings > 0) {
    FORLIM = n_rings;
    for (i = 0; i < FORLIM; i++) {
      for (j = 0; j < max_ringsize; j++)
	tmp_ring[i][j] = ring[i][j];
    }
    for (i = 0; i < max_rings; i++) {   /* new in v0.3 
      tmp_ringprop[i].size = ringprop[i].size;
      tmp_ringprop[i].arom = ringprop[i].arom;
      tmp_ringprop[i].envelope = ringprop[i].envelope;
    }
  }
  tmp_molstat.n_QA = molstat.n_QA;
  tmp_molstat.n_QB = molstat.n_QB;
  tmp_molstat.n_chg = molstat.n_chg;
  tmp_molstat.n_C1 = molstat.n_C1;
  tmp_molstat.n_C2 = molstat.n_C2;
  tmp_molstat.n_C = molstat.n_C;
  tmp_molstat.n_CHB1p = molstat.n_CHB1p;
  tmp_molstat.n_CHB2p = molstat.n_CHB2p;
  tmp_molstat.n_CHB3p = molstat.n_CHB3p;
  tmp_molstat.n_CHB4 = molstat.n_CHB4;
  tmp_molstat.n_O2 = molstat.n_O2;
  tmp_molstat.n_O3 = molstat.n_O3;
  tmp_molstat.n_N1 = molstat.n_N1;
  tmp_molstat.n_N2 = molstat.n_N2;
  tmp_molstat.n_N3 = molstat.n_N3;
  tmp_molstat.n_S = molstat.n_S;
  tmp_molstat.n_SeTe = molstat.n_SeTe;
  tmp_molstat.n_F = molstat.n_F;
  tmp_molstat.n_Cl = molstat.n_Cl;
  tmp_molstat.n_Br = molstat.n_Br;
  tmp_molstat.n_I = molstat.n_I;
  tmp_molstat.n_P = molstat.n_P;
  tmp_molstat.n_B = molstat.n_B;
  tmp_molstat.n_Met = molstat.n_Met;
  tmp_molstat.n_X = molstat.n_X;
  tmp_molstat.n_b1 = molstat.n_b1;
  tmp_molstat.n_b2 = molstat.n_b2;
  tmp_molstat.n_b3 = molstat.n_b3;
  tmp_molstat.n_bar = molstat.n_bar;
  tmp_molstat.n_C1O = molstat.n_C1O;
  tmp_molstat.n_C2O = molstat.n_C2O;
  tmp_molstat.n_CN = molstat.n_CN;
  tmp_molstat.n_XY = molstat.n_XY;
  tmp_molstat.n_r3 = molstat.n_r3;
  tmp_molstat.n_r4 = molstat.n_r4;
  tmp_molstat.n_r5 = molstat.n_r5;
  tmp_molstat.n_r6 = molstat.n_r6;
  tmp_molstat.n_r7 = molstat.n_r7;
  tmp_molstat.n_r8 = molstat.n_r8;
  tmp_molstat.n_r9 = molstat.n_r9;
  tmp_molstat.n_r10 = molstat.n_r10;
  tmp_molstat.n_r11 = molstat.n_r11;
  tmp_molstat.n_r12 = molstat.n_r12;
  tmp_molstat.n_r13p = molstat.n_r13p;
  tmp_molstat.n_rN = molstat.n_rN;
  tmp_molstat.n_rN1 = molstat.n_rN1;
  tmp_molstat.n_rN2 = molstat.n_rN2;
  tmp_molstat.n_rN3p = molstat.n_rN3p;
  tmp_molstat.n_rO = molstat.n_rO;
  tmp_molstat.n_rO1 = molstat.n_rO1;
  tmp_molstat.n_rO2p = molstat.n_rO2p;
  tmp_molstat.n_rS = molstat.n_rS;
  tmp_molstat.n_rX = molstat.n_rX;
  tmp_molstat.n_rAr = molstat.n_rAr;
  tmp_molstat.n_rBz = molstat.n_rBz;   /* v0.3l 
  tmp_molstat.n_br2p = molstat.n_br2p;   /* v0.3n 
/* p2c: checkmol.pas, line 8022:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] 
  /*$IFDEF extended_molstat
  /* v0.3m 
  tmp_molstat.n_psg01 = molstat.n_psg01;
  tmp_molstat.n_psg02 = molstat.n_psg02;
  tmp_molstat.n_psg13 = molstat.n_psg13;
  tmp_molstat.n_psg14 = molstat.n_psg14;
  tmp_molstat.n_psg15 = molstat.n_psg15;
  tmp_molstat.n_psg16 = molstat.n_psg16;
  tmp_molstat.n_psg17 = molstat.n_psg17;
  tmp_molstat.n_psg18 = molstat.n_psg18;
  tmp_molstat.n_pstm = molstat.n_pstm;
  tmp_molstat.n_psla = molstat.n_psla;
  /*$ENDIF
  /* make sure some modes can be switched on only by the query file 
  /* and not by subsequent haystack file(s) 
  if (ez_flag)   /* new in v0.3f 
    ez_search = true;
  if (chir_flag)   /* new in v0.3f 
    rs_search = true;
} */

static void
copy_tmp_to_mol ()
{
  //int i, j, FORLIM;

  /*if (n_atoms == 0)
     return; *///If a NoStruct is read, this leads to madness and illegal memory access


  atom = (atom_rec *) safe_calloc (n_atoms, sizeof (atom_rec));
  bond = (bond_rec *) safe_calloc (n_bonds, sizeof (bond_rec));
  ring = (ringpath_type *) safe_calloc (1, sizeof (ringlist));
  ringprop = (ringprop_rec *) safe_calloc (1, sizeof (ringprop_type));


  n_atoms = tmp_n_atoms;
  n_bonds = tmp_n_bonds;
  n_rings = tmp_n_rings;
  n_heavyatoms = tmp_n_heavyatoms;
  n_heavybonds = tmp_n_heavybonds;
  strcpy (molname, tmp_molname);
  n_Ctot = tmp_n_Ctot;
  n_Otot = tmp_n_Otot;
  n_Ntot = tmp_n_Ntot;
  memcpy (atom, tmp_atom, tmp_n_atoms * sizeof (atom_rec));

  if (tmp_n_bonds > 0)
    memcpy (bond, tmp_bond, tmp_n_bonds * sizeof (bond_rec));

  if (tmp_n_rings > 0)
    {
      memcpy (ring, tmp_ring, sizeof (ringlist));
      memcpy (ringprop, tmp_ringprop, sizeof (ringprop_type));
    }

  memcpy (&molstat, &tmp_molstat, sizeof (tmp_molstat));


  // make sure some modes can be switched on only by the query file 
  // and not by subsequent haystack file(s) 
  if (ez_flag)			// new in v0.3f 
    ez_search = true;

  if (chir_flag)		// new in v0.3f 
    rs_search = true;

}

/*static void copy_tmp_to_mol()
{
  int i, j, FORLIM;

  if (tmp_n_atoms == 0)
    return;
  n_atoms = tmp_n_atoms;
  n_bonds = tmp_n_bonds;
  n_rings = tmp_n_rings;
  n_heavyatoms = tmp_n_heavyatoms;
  n_heavybonds = tmp_n_heavybonds;
  strcpy(molname, tmp_molname);
  n_Ctot = tmp_n_Ctot;
  n_Otot = tmp_n_Otot;
  n_Ntot = tmp_n_Ntot;
  /* try 
  atom = (atom_rec *)safe_malloc(n_atoms * sizeof(atom_rec));
  bond = (bond_rec *)safe_malloc(n_bonds * sizeof(bond_rec));
  ring = (ringpath_type *)safe_malloc(sizeof(ringlist));
  ringprop = (ringprop_rec *)safe_malloc(sizeof(ringprop_type));
  FORLIM = tmp_n_atoms;
  /* except
    on e:Eoutofmemory do
      begin
        writeln('Not enough memory');
        halt(4);
      end;
  end; 
  for (i = 0; i < FORLIM; i++) {
    strcpy(atom[i].element, tmp_atom[i].element);
    strcpy(atom[i].atype, tmp_atom[i].atype);
    atom[i].x = tmp_atom[i].x;
    atom[i].y = tmp_atom[i].y;
    atom[i].z = tmp_atom[i].z;
    atom[i].formal_charge = tmp_atom[i].formal_charge;
    atom[i].real_charge = tmp_atom[i].real_charge;
    atom[i].Hexp = tmp_atom[i].Hexp;
    atom[i].Htot = tmp_atom[i].Htot;
    atom[i].neighbor_count = tmp_atom[i].neighbor_count;
    atom[i].ring_count = tmp_atom[i].ring_count;
    atom[i].arom = tmp_atom[i].arom;
    atom[i].stereo_care = tmp_atom[i].stereo_care;
    atom[i].heavy = tmp_atom[i].heavy;   /* v0.3l 
    atom[i].metal = tmp_atom[i].metal;   /* v0.3l 
    atom[i].tag = tmp_atom[i].tag;   /* v0.3o 
  }
  if (tmp_n_bonds > 0) {
    FORLIM = tmp_n_bonds;
    for (i = 0; i < FORLIM; i++) {
      bond[i].a1 = tmp_bond[i].a1;
      bond[i].a2 = tmp_bond[i].a2;
      bond[i].btype = tmp_bond[i].btype;
      bond[i].arom = tmp_bond[i].arom;
      bond[i].ring_count = tmp_bond[i].ring_count;   /* new in v0.3d 
      bond[i].topo = tmp_bond[i].topo;   /* new in v0.3d 
      bond[i].stereo = tmp_bond[i].stereo;   /* new in v0.3d 
    }
  }
  if (tmp_n_rings > 0) {
    FORLIM = tmp_n_rings;
    for (i = 0; i < FORLIM; i++) {
      for (j = 0; j < max_ringsize; j++)
	ring[i][j] = tmp_ring[i][j];
    }
    for (i = 0; i < max_rings; i++) {   /* new in v0.3 
      ringprop[i].size = tmp_ringprop[i].size;
      ringprop[i].arom = tmp_ringprop[i].arom;
      ringprop[i].envelope = tmp_ringprop[i].envelope;
    }
  }
  molstat.n_QA = tmp_molstat.n_QA;
  molstat.n_QB = tmp_molstat.n_QB;
  molstat.n_chg = tmp_molstat.n_chg;
  molstat.n_C1 = tmp_molstat.n_C1;
  molstat.n_C2 = tmp_molstat.n_C2;
  molstat.n_C = tmp_molstat.n_C;
  molstat.n_CHB1p = tmp_molstat.n_CHB1p;
  molstat.n_CHB2p = tmp_molstat.n_CHB2p;
  molstat.n_CHB3p = tmp_molstat.n_CHB3p;
  molstat.n_CHB4 = tmp_molstat.n_CHB4;
  molstat.n_O2 = tmp_molstat.n_O2;
  molstat.n_O3 = tmp_molstat.n_O3;
  molstat.n_N1 = tmp_molstat.n_N1;
  molstat.n_N2 = tmp_molstat.n_N2;
  molstat.n_N3 = tmp_molstat.n_N3;
  molstat.n_S = tmp_molstat.n_S;
  molstat.n_SeTe = tmp_molstat.n_SeTe;
  molstat.n_F = tmp_molstat.n_F;
  molstat.n_Cl = tmp_molstat.n_Cl;
  molstat.n_Br = tmp_molstat.n_Br;
  molstat.n_I = tmp_molstat.n_I;
  molstat.n_P = tmp_molstat.n_P;
  molstat.n_B = tmp_molstat.n_B;
  molstat.n_Met = tmp_molstat.n_Met;
  molstat.n_X = tmp_molstat.n_X;
  molstat.n_b1 = tmp_molstat.n_b1;
  molstat.n_b2 = tmp_molstat.n_b2;
  molstat.n_b3 = tmp_molstat.n_b3;
  molstat.n_bar = tmp_molstat.n_bar;
  molstat.n_C1O = tmp_molstat.n_C1O;
  molstat.n_C2O = tmp_molstat.n_C2O;
  molstat.n_CN = tmp_molstat.n_CN;
  molstat.n_XY = tmp_molstat.n_XY;
  molstat.n_r3 = tmp_molstat.n_r3;
  molstat.n_r4 = tmp_molstat.n_r4;
  molstat.n_r5 = tmp_molstat.n_r5;
  molstat.n_r6 = tmp_molstat.n_r6;
  molstat.n_r7 = tmp_molstat.n_r7;
  molstat.n_r8 = tmp_molstat.n_r8;
  molstat.n_r9 = tmp_molstat.n_r9;
  molstat.n_r10 = tmp_molstat.n_r10;
  molstat.n_r11 = tmp_molstat.n_r11;
  molstat.n_r12 = tmp_molstat.n_r12;
  molstat.n_r13p = tmp_molstat.n_r13p;
  molstat.n_rN = tmp_molstat.n_rN;
  molstat.n_rN1 = tmp_molstat.n_rN1;
  molstat.n_rN2 = tmp_molstat.n_rN2;
  molstat.n_rN3p = tmp_molstat.n_rN3p;
  molstat.n_rO = tmp_molstat.n_rO;
  molstat.n_rO1 = tmp_molstat.n_rO1;
  molstat.n_rO2p = tmp_molstat.n_rO2p;
  molstat.n_rS = tmp_molstat.n_rS;
  molstat.n_rX = tmp_molstat.n_rX;
  molstat.n_rAr = tmp_molstat.n_rAr;
  molstat.n_rBz = tmp_molstat.n_rBz;   /* v0.3l 
  molstat.n_br2p = tmp_molstat.n_br2p;   /* v0.3n 
/* p2c: checkmol.pas, line 8169:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
  /*$IFDEF extended_molstat
     molstat.n_psg01 = tmp_molstat.n_psg01;
     molstat.n_psg02 = tmp_molstat.n_psg02;
     molstat.n_psg13 = tmp_molstat.n_psg13;
     molstat.n_psg14 = tmp_molstat.n_psg14;
     molstat.n_psg15 = tmp_molstat.n_psg15;
     molstat.n_psg16 = tmp_molstat.n_psg16;
     molstat.n_psg17 = tmp_molstat.n_psg17;
     molstat.n_psg18 = tmp_molstat.n_psg18;
     molstat.n_pstm = tmp_molstat.n_pstm;
     molstat.n_psla = tmp_molstat.n_psla;
     /*$ENDIF */
  /* make sure some modes can be switched on only by the query file */
  /* and not by subsequent haystack file(s) 
     if (ez_flag)
     ez_search = true;
     if (chir_flag)
     rs_search = true;
     } */





static void
get_ringstat (r_id)
     int r_id;
{
  int i, j;
  ringpath_type testring;
  int ring_size, a_ref;
  str2 elem;
  int nN = 0, nO = 0, nS = 0, nX = 0;

  if (r_id < 1 || r_id > n_rings)
    return;
  memset (testring, 0, sizeof (ringpath_type));
  ring_size = ringprop[r_id - 1].size;	/* v0.3j */
  for (j = 0; j < ring_size; j++)	/* v0.3j */
    testring[j] = ring[r_id - 1][j];
/* p2c: checkmol.pas, line 8238:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
#ifdef reduced_SAR
  if (ring_size <= 2 || ringprop[r_id - 1].envelope != false)
    /* v0.3n: ignore envelope rings */
    return;
#else
  if (ring_size <= 2)
    return;
#endif
  for (i = 0; i < ring_size; i++)
    {
      a_ref = testring[i];
      strcpy (elem, atom[a_ref - 1].element);
      if (strcmp (elem, "C ") && strcmp (elem, "A "))
	{
	  nX++;			/* general heteroatom count */
	  if (!strcmp (elem, "N "))
	    nN++;
	  if (!strcmp (elem, "O "))
	    nO++;
	  if (!strcmp (elem, "S "))
	    nS++;
	}
    }
  if (nN > 0)
    {
      molstat.n_rN++;
      if (nN == 1)
	molstat.n_rN1++;
      if (nN == 2)
	molstat.n_rN2++;
      if (nN > 2)
	molstat.n_rN3p++;
    }
  if (nO > 0)
    {
      molstat.n_rO++;
      if (nO == 1)
	molstat.n_rO1++;
      if (nO == 2)
	molstat.n_rO2p++;
    }
  if (nS > 0)
    molstat.n_rS++;
  if (nX > 0)
    molstat.n_rX++;
  /* general ringsize descriptors; v0.3m */
  switch (ring_size)
    {

    case 3:
      molstat.n_r3++;
      break;

    case 4:
      molstat.n_r4++;
      break;

    case 5:
      molstat.n_r5++;
      break;

    case 6:
      molstat.n_r6++;
      break;

    case 7:
      molstat.n_r7++;
      break;

    case 8:
      molstat.n_r8++;
      break;

    case 9:
      molstat.n_r9++;
      break;

    case 10:
      molstat.n_r10++;
      break;

    case 11:
      molstat.n_r11++;
      break;

    case 12:
      molstat.n_r12++;
      break;

    default:
      molstat.n_r13p++;
      break;
    }				/* end v0.3m        */
}


static void
get_molstat ()
{
  int i;
  str2 elem;
  str3 atype;
  int a1, a2;
  str2 a1el, a2el;
  char btype;
  int hbc;
  int n_b2formal = 0;		/* new in v0.2e */
  int FORLIM;

  if (n_atoms == 0)
    return;
  FORLIM = n_atoms;
  for (i = 0; i < FORLIM; i++)
    {
      if (atom[i].heavy)
	{
	  strcpy (elem, atom[i].element);
	  strcpy (atype, atom[i].atype);
	  if (!strcmp (atype, "C1 "))
	    molstat.n_C1++;
	  if (!strcmp (atype, "C2 ") || !strcmp (atype, "CAR"))
	    molstat.n_C2++;
	  if (!strcmp (elem, "C "))
	    molstat.n_C++;
	  if (!strcmp (atype, "O2 "))
	    molstat.n_O2++;
	  if (!strcmp (atype, "O3 "))
	    molstat.n_O3++;
	  if (!strcmp (atype, "N1 "))
	    molstat.n_N1++;
	  if (!strcmp (atype, "N2 ") || !strcmp (atype, "NAR") ||
	      !strcmp (atype, "NAM") && atom[i].arom == true)
	    /* v0.3n */
	    molstat.n_N2++;
	  if (!strcmp (atype, "N3 ") || !strcmp (atype, "NPL") ||
	      !strcmp (atype, "N3+") ||
	      !strcmp (atype, "NAM") && atom[i].arom == false)
	    /* v0.3n */
	    molstat.n_N3++;
	  if (!strcmp (elem, "A "))	/* query atom */
	    molstat.n_QA++;
	  if (!strcmp (elem, "Q "))	/* query atom */
	    molstat.n_QA++;
	  if (!strcmp (elem, "X "))	/* query atom */
	    molstat.n_QA++;
	  if (!strcmp (elem, "S "))
	    molstat.n_S++;
	  if (!strcmp (elem, "SE"))
	    molstat.n_SeTe++;
	  if (!strcmp (elem, "TE"))
	    molstat.n_SeTe++;
	  if (!strcmp (elem, "F "))
	    molstat.n_F++;
	  if (!strcmp (elem, "CL"))
	    molstat.n_Cl++;
	  if (!strcmp (elem, "BR"))
	    molstat.n_Br++;
	  if (!strcmp (elem, "I "))
	    molstat.n_I++;
	  if (!strcmp (elem, "P "))
	    molstat.n_P++;
	  if (!strcmp (elem, "B "))
	    molstat.n_B++;
	  /* check for known metals */
	  if (atom[i].metal)	/* v0.3l */
	    molstat.n_Met++;
	  /* still missing: unknown elements */

	  /* check number of heteroatom bonds per C atom */
	  if (!strcmp (elem, "C "))
	    {
	      hbc = raw_hetbond_count (i + 1);
	      /* new in v0.2j (replaces hetbond_count) */
	      if (hbc >= 1)
		molstat.n_CHB1p++;
	      if (hbc >= 2)
		molstat.n_CHB2p++;
	      if (hbc >= 3)
		molstat.n_CHB3p++;
	      if (hbc == 4)
		molstat.n_CHB4++;
	    }
	  if (atom[i].formal_charge != 0)
	    {
	      molstat.n_chg++;
	      //n_charges++;
	    }
	  if (atom[i].nucleon_number != 0)
	    {
	      molstat.n_iso++;
	    }
	  if (atom[i].radical_type != 0)
	    {
	      molstat.n_rad++;
	    }
	  /* check for "other" elements;  v0.3l */
	  if (!atom[i].metal && strcmp (elem, "C ") && strcmp (elem, "N ")
	      && strcmp (elem, "O ") && strcmp (elem, "S ")
	      && strcmp (elem, "SE") && strcmp (elem, "TE")
	      && strcmp (elem, "P ") && strcmp (elem, "B ")
	      && strcmp (elem, "A ") && strcmp (elem, "Q "))
	    molstat.n_X++;
	  /*(elem = 'F ') or (elem = 'CL') or (elem = 'BR') or (elem = 'I ') or  (* leave halogens as type X, v0.3m */
/* p2c: checkmol.pas, line 8353:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
	  /*$IFDEF extended_molstat */
	  if (!strcmp (elem, "LI") || !strcmp (elem, "NA")
	      || !strcmp (elem, "K ") || !strcmp (elem, "RB")
	      || !strcmp (elem, "CS") || !strcmp (elem, "FR"))
	    molstat.n_psg01++;
	  if (!strcmp (elem, "BE") || !strcmp (elem, "MG")
	      || !strcmp (elem, "CA") || !strcmp (elem, "SR")
	      || !strcmp (elem, "BA") || !strcmp (elem, "RA"))
	    molstat.n_psg02++;
	  if (!strcmp (elem, "B ") || !strcmp (elem, "AL")
	      || !strcmp (elem, "GA") || !strcmp (elem, "IN")
	      || !strcmp (elem, "TL"))
	    molstat.n_psg13++;
	  if (!strcmp (elem, "C ") || !strcmp (elem, "SI")
	      || !strcmp (elem, "GE") || !strcmp (elem, "SN")
	      || !strcmp (elem, "PB"))
	    molstat.n_psg14++;
	  if (!strcmp (elem, "N ") || !strcmp (elem, "P ")
	      || !strcmp (elem, "AS") || !strcmp (elem, "SB")
	      || !strcmp (elem, "BI"))
	    molstat.n_psg15++;
	  if (!strcmp (elem, "O ") || !strcmp (elem, "S ")
	      || !strcmp (elem, "SE") || !strcmp (elem, "TE")
	      || !strcmp (elem, "PO"))
	    molstat.n_psg16++;
	  if (!strcmp (elem, "F ") || !strcmp (elem, "CL")
	      || !strcmp (elem, "BR") || !strcmp (elem, "I ")
	      || !strcmp (elem, "AT"))
	    molstat.n_psg17++;
	  if (!strcmp (elem, "HE") || !strcmp (elem, "NE")
	      || !strcmp (elem, "AR") || !strcmp (elem, "KR")
	      || !strcmp (elem, "XE") || !strcmp (elem, "RN"))
	    molstat.n_psg18++;
	  if (!strcmp (elem, "SC") || !strcmp (elem, "Y ")
	      || !strcmp (elem, "LU") || !strcmp (elem, "LR")
	      || !strcmp (elem, "TI") || !strcmp (elem, "ZR")
	      || !strcmp (elem, "HF") || !strcmp (elem, "RF")
	      || !strcmp (elem, "V ") || !strcmp (elem, "NB")
	      || !strcmp (elem, "TA") || !strcmp (elem, "DB")
	      || !strcmp (elem, "CR") || !strcmp (elem, "MO")
	      || !strcmp (elem, "W ") || !strcmp (elem, "SG")
	      || !strcmp (elem, "MN") || !strcmp (elem, "TC")
	      || !strcmp (elem, "RE") || !strcmp (elem, "BH")
	      || !strcmp (elem, "FE") || !strcmp (elem, "RU")
	      || !strcmp (elem, "OS") || !strcmp (elem, "HS")
	      || !strcmp (elem, "CO") || !strcmp (elem, "RH")
	      || !strcmp (elem, "IR") || !strcmp (elem, "MT")
	      || !strcmp (elem, "NI") || !strcmp (elem, "PD")
	      || !strcmp (elem, "PT") || !strcmp (elem, "DS")
	      || !strcmp (elem, "CU") || !strcmp (elem, "AG")
	      || !strcmp (elem, "AU") || !strcmp (elem, "RG")
	      || !strcmp (elem, "ZN") || !strcmp (elem, "CD")
	      || !strcmp (elem, "HG"))
/* p2c: checkmol.pas, line 8439: 
 * Note: Line breaker spent 0.0 seconds, 5000 tries on line 10035 [251] */
	    molstat.n_pstm++;
	  if (!strcmp (elem, "LA") || !strcmp (elem, "CE")
	      || !strcmp (elem, "PR") || !strcmp (elem, "ND")
	      || !strcmp (elem, "PM") || !strcmp (elem, "SM")
	      || !strcmp (elem, "EU") || !strcmp (elem, "GD")
	      || !strcmp (elem, "TB") || !strcmp (elem, "DY")
	      || !strcmp (elem, "HO") || !strcmp (elem, "ER")
	      || !strcmp (elem, "TM") || !strcmp (elem, "YB")
	      || !strcmp (elem, "AC") || !strcmp (elem, "TH")
	      || !strcmp (elem, "PA") || !strcmp (elem, "U ")
	      || !strcmp (elem, "NP") || !strcmp (elem, "PU")
	      || !strcmp (elem, "AM") || !strcmp (elem, "CM")
	      || !strcmp (elem, "BK") || !strcmp (elem, "CF")
	      || !strcmp (elem, "ES") || !strcmp (elem, "FM")
	      || !strcmp (elem, "MD") || !strcmp (elem, "NO"))
/* p2c: checkmol.pas, line 8439: 
 * Note: Line breaker spent 0.0 seconds, 5000 tries on line 10048 [251] */
	    molstat.n_psla++;
	  /*$ENDIF */
	}			/* is heavy */
    }				/* atoms */
  if (n_bonds > 0)
    {
      FORLIM = n_bonds;
      for (i = 0; i < FORLIM; i++)
	{
	  a1 = bond[i].a1;
	  a2 = bond[i].a2;
	  strcpy (a1el, atom[a1 - 1].element);
	  strcpy (a2el, atom[a2 - 1].element);
	  btype = bond[i].btype;
	  if (bond[i].arom)
	    molstat.n_bar++;
	  else
	    {
	      if (btype == 'S' && atom[a1 - 1].heavy && atom[a2 - 1].heavy)
		molstat.n_b1++;
	      if (btype == 'D')
		molstat.n_b2++;
	      if (btype == 'T')
		molstat.n_b3++;
	    }
	  /* v0.3n: ignore bonds to (explicit) hydrogens */
	  if (!strcmp (a1el, "C ") && !strcmp (a2el, "O ") ||
	      !strcmp (a1el, "O ") && !strcmp (a2el, "C "))
	    {
	      if (btype == 'S')
		molstat.n_C1O++;
	      if (btype == 'D')
		molstat.n_C2O++;
	    }
	  if (!strcmp (a1el, "C ") && !strcmp (a2el, "N ") ||
	      !strcmp (a1el, "N ") && !strcmp (a2el, "C "))
	    molstat.n_CN++;
	  if (strcmp (a1el, "C ") && atom[a1 - 1].heavy
	      && strcmp (a2el, "C ") && atom[a2 - 1].heavy)
	    molstat.n_XY++;
	  /* new in v0.3n: number of bonds belonging to more than one ring */
	  if (bond[i].ring_count > 1)
	    molstat.n_br2p++;
	}
    }				/* bonds */
  if (n_rings <= 0)
    {
      return;
    }				/* rings */
  /* v0.3n */
  n_countablerings = 0;		/* v0.3n */
  FORLIM = n_rings;
  for (i = 1; i <= FORLIM; i++)
    {
      if (ringprop[i - 1].envelope == false)	/* v0.3n */
	n_countablerings++;
      if (is_arene (i) && ringprop[i - 1].envelope == false)
	{			/* v0.3n: ignore envelope rings */
	  molstat.n_rAr++;
	  if ((ringprop[i - 1].size == 6) && (is_heterocycle (i) == false))
	    /* v0.3l */
	    molstat.n_rBz++;
	}
      get_ringstat (i);
      if (ringprop[i - 1].arom == true && ringprop[i - 1].envelope == false)
	/* new in v0.3n; replaces assignment below */
	n_b2formal++;
    }
  /*n_b2formal := n_rar;  (* new in v0.2e; adds 1 formal double bond for each aromatic ring */
  /* in order to allow an isolated double bond in the needle */
  /* to be matched as a ring fragment of an aromatic ring */
  if (n_b2formal > molstat.n_bar / 2)
    n_b2formal = molstat.n_bar / 2;
  molstat.n_b2 += n_b2formal;
}


static void
fix_ssr_ringcounts ()
{
  /* new in v0.3n */
  /* if SAR -> SSR fallback happens, set some molstat values */
  /* to a maximum (ring counts for various ring sizes); */
  /* this should be necessary only for ring sizes which */
  /* are a) too large for the SSR (depending on ssr_vringsize) */
  /* and b) which are likely to contain "envelope rings" */
  /* (size 6 and above) */
  /*  if (molstat.n_r3 = 0) then molstat.n_r3 := max_rings; */
  /*  if (molstat.n_r4 = 0) then molstat.n_r4 := max_rings; */
  /*  if (molstat.n_r5 = 0) then molstat.n_r5 := max_rings; */
  if (molstat.n_r6 == 0)
    molstat.n_r6 = max_rings;
  if (molstat.n_r7 == 0)
    molstat.n_r7 = max_rings;
  if (molstat.n_r8 == 0)
    molstat.n_r8 = max_rings;
  if (molstat.n_r9 == 0)
    molstat.n_r9 = max_rings;
  if (molstat.n_r10 == 0)
    molstat.n_r10 = max_rings;
  if (molstat.n_r11 == 0)
    molstat.n_r11 = max_rings;
  if (molstat.n_r12 == 0)
    molstat.n_r12 = max_rings;
  if (molstat.n_r13p == 0)
    molstat.n_r13p = max_rings;
}



static void
write_molstat ()
{
  if (auto_ssr)			/* v0.3n */
    fix_ssr_ringcounts ();
  printf ("n_atoms:%d;", n_heavyatoms);
  /* count only non-H atoms (some molfiles contain explicit H's) */
  if (n_bonds > 0)		/* count only bonds between non-H atoms */
    printf ("n_bonds:%d;", n_heavybonds);
/* p2c: checkmol.pas, line 8471:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
#ifdef REDUCED_SAR
  if (n_rings > 0)		/* changed to non-envelope rings in v0.3n */
    printf ("n_rings:%d;", n_countablerings);
#else
  if (n_rings > 0)		/* changed to non-envelope rings in v0.3n */
    printf ("n_rings:%d;", n_rings);
#endif
  /*      if n_QA    > 0 then write('n_QA:',n_QA,';'); */
  /*      if n_QB    > 0 then write('n_QB:',n_QB,';'); */
  if (molstat.n_chg > 0)	/* 0.3x */
    printf ("n_chg:%d;", molstat.n_chg);

  if (molstat.n_C1 > 0)
    printf ("n_C1:%d;", molstat.n_C1);
  if (molstat.n_C2 > 0)
    printf ("n_C2:%d;", molstat.n_C2);
  /* requirement of a given number of sp3 carbons might be too restrictive, */
  /* so we use the total number of carbons instead  (initially used variable n_C3 is now n_C) */
  if (molstat.n_C > 0)
    printf ("n_C:%d;", molstat.n_C);
  if (molstat.n_CHB1p > 0)
    printf ("n_CHB1p:%d;", molstat.n_CHB1p);
  if (molstat.n_CHB2p > 0)
    printf ("n_CHB2p:%d;", molstat.n_CHB2p);
  if (molstat.n_CHB3p > 0)
    printf ("n_CHB3p:%d;", molstat.n_CHB3p);
  if (molstat.n_CHB4 > 0)
    printf ("n_CHB4:%d;", molstat.n_CHB4);
  if (molstat.n_O2 > 0)
    printf ("n_O2:%d;", molstat.n_O2);
  if (molstat.n_O3 > 0)
    printf ("n_O3:%d;", molstat.n_O3);
  if (molstat.n_N1 > 0)
    printf ("n_N1:%d;", molstat.n_N1);
  if (molstat.n_N2 > 0)
    printf ("n_N2:%d;", molstat.n_N2);
  if (molstat.n_N3 > 0)
    printf ("n_N3:%d;", molstat.n_N3);
  if (molstat.n_S > 0)
    printf ("n_S:%d;", molstat.n_S);
  if (molstat.n_SeTe > 0)
    printf ("n_SeTe:%d;", molstat.n_SeTe);
  if (molstat.n_F > 0)
    printf ("n_F:%d;", molstat.n_F);
  if (molstat.n_Cl > 0)
    printf ("n_Cl:%d;", molstat.n_Cl);
  if (molstat.n_Br > 0)
    printf ("n_Br:%d;", molstat.n_Br);
  if (molstat.n_I > 0)
    printf ("n_I:%d;", molstat.n_I);
  if (molstat.n_P > 0)
    printf ("n_P:%d;", molstat.n_P);
  if (molstat.n_B > 0)
    printf ("n_B:%d;", molstat.n_B);
  if (molstat.n_Met > 0)
    printf ("n_Met:%d;", molstat.n_Met);
  if (molstat.n_X > 0)
    printf ("n_X:%d;", molstat.n_X);
  if (molstat.n_b1 > 0)
    printf ("n_b1:%d;", molstat.n_b1);
  if (molstat.n_b2 > 0)
    printf ("n_b2:%d;", molstat.n_b2);
  if (molstat.n_b3 > 0)
    printf ("n_b3:%d;", molstat.n_b3);
  if (molstat.n_bar > 0)
    printf ("n_bar:%d;", molstat.n_bar);
  if (molstat.n_C1O > 0)
    printf ("n_C1O:%d;", molstat.n_C1O);
  if (molstat.n_C2O > 0)
    printf ("n_C2O:%d;", molstat.n_C2O);
  if (molstat.n_CN > 0)
    printf ("n_CN:%d;", molstat.n_CN);
  if (molstat.n_XY > 0)
    printf ("n_XY:%d;", molstat.n_XY);
  if (molstat.n_r3 > 0)
    printf ("n_r3:%d;", molstat.n_r3);
  if (molstat.n_r4 > 0)
    printf ("n_r4:%d;", molstat.n_r4);
  if (molstat.n_r5 > 0)
    printf ("n_r5:%d;", molstat.n_r5);
  if (molstat.n_r6 > 0)
    printf ("n_r6:%d;", molstat.n_r6);
  if (molstat.n_r7 > 0)
    printf ("n_r7:%d;", molstat.n_r7);
  if (molstat.n_r8 > 0)
    printf ("n_r8:%d;", molstat.n_r8);
  if (molstat.n_r9 > 0)
    printf ("n_r9:%d;", molstat.n_r9);
  if (molstat.n_r10 > 0)
    printf ("n_r10:%d;", molstat.n_r10);
  if (molstat.n_r11 > 0)
    printf ("n_r11:%d;", molstat.n_r11);
  if (molstat.n_r12 > 0)
    printf ("n_r12:%d;", molstat.n_r12);
  if (molstat.n_r13p > 0)
    printf ("n_r13p:%d;", molstat.n_r13p);
  if (molstat.n_rN > 0)
    printf ("n_rN:%d;", molstat.n_rN);
  if (molstat.n_rN1 > 0)
    printf ("n_rN1:%d;", molstat.n_rN1);
  if (molstat.n_rN2 > 0)
    printf ("n_rN2:%d;", molstat.n_rN2);
  if (molstat.n_rN3p > 0)
    printf ("n_rN3p:%d;", molstat.n_rN3p);
  if (molstat.n_rO > 0)
    printf ("n_rO:%d;", molstat.n_rO);
  if (molstat.n_rO1 > 0)
    printf ("n_rO1:%d;", molstat.n_rO1);
  if (molstat.n_rO2p > 0)
    printf ("n_rO2p:%d;", molstat.n_rO2p);
  if (molstat.n_rS > 0)
    printf ("n_rS:%d;", molstat.n_rS);
  if (molstat.n_rX > 0)
    printf ("n_rX:%d;", molstat.n_rX);
  if (molstat.n_rAr > 0)
    printf ("n_rar:%d;", molstat.n_rAr);
/* p2c: checkmol.pas, line 8532:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
  /*$IFDEF extended_molstat */
  if (molstat.n_rBz > 0)
    printf ("n_rbz:%d;", molstat.n_rBz);
  if (molstat.n_br2p > 0)
    printf ("n_br2p:%d;", molstat.n_br2p);
  if (molstat.n_psg01 > 0)
    printf ("n_psg01:%d;", molstat.n_psg01);
  if (molstat.n_psg02 > 0)
    printf ("n_psg02:%d;", molstat.n_psg02);
  if (molstat.n_psg13 > 0)
    printf ("n_psg13:%d;", molstat.n_psg13);
  if (molstat.n_psg14 > 0)
    printf ("n_psg14:%d;", molstat.n_psg14);
  if (molstat.n_psg15 > 0)
    printf ("n_psg15:%d;", molstat.n_psg15);
  if (molstat.n_psg16 > 0)
    printf ("n_psg16:%d;", molstat.n_psg16);
  if (molstat.n_psg17 > 0)
    printf ("n_psg17:%d;", molstat.n_psg17);
  if (molstat.n_psg18 > 0)
    printf ("n_psg18:%d;", molstat.n_psg18);
  if (molstat.n_pstm > 0)
    printf ("n_pstm:%d;", molstat.n_pstm);
  if (molstat.n_psla > 0)
    printf ("n_psla:%d;", molstat.n_psla);
  if (molstat.n_iso > 0)
    printf ("n_iso:%d;", molstat.n_iso);
  if (molstat.n_rad > 0)
    printf ("n_rad:%d;", molstat.n_rad);
  /*$ENDIF */
  putchar ('\n');
}


static void
write_molstat_X ()
{
  if (auto_ssr)			/* v0.3n */
    fix_ssr_ringcounts ();
/* p2c: checkmol.pas, line 8556:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
#ifdef REDUCED_SAR
  printf ("%d,", n_heavyatoms);
  printf ("%d,", n_heavybonds);
  printf ("%d,", n_countablerings);
  /* v0.3n: n_rings =?> n_countablerings */
#else
  printf ("%d,", n_heavyatoms);
  printf ("%d,", n_heavybonds);
  printf ("%d,", n_rings);	/* v0.3n: n_rings ==> n_countablerings */
#endif
  printf ("%d,", molstat.n_QA);
  printf ("%d,", molstat.n_QB);

  /* 0.3x */
  printf ("%d,", molstat.n_chg);


  printf ("%d,", molstat.n_C1);
  printf ("%d,", molstat.n_C2);
  printf ("%d,", molstat.n_C);
  printf ("%d,", molstat.n_CHB1p);
  printf ("%d,", molstat.n_CHB2p);
  printf ("%d,", molstat.n_CHB3p);
  printf ("%d,", molstat.n_CHB4);
  printf ("%d,", molstat.n_O2);
  printf ("%d,", molstat.n_O3);
  printf ("%d,", molstat.n_N1);
  printf ("%d,", molstat.n_N2);
  printf ("%d,", molstat.n_N3);
  printf ("%d,", molstat.n_S);
  printf ("%d,", molstat.n_SeTe);
  printf ("%d,", molstat.n_F);
  printf ("%d,", molstat.n_Cl);
  printf ("%d,", molstat.n_Br);
  printf ("%d,", molstat.n_I);
  printf ("%d,", molstat.n_P);
  printf ("%d,", molstat.n_B);
  printf ("%d,", molstat.n_Met);
  printf ("%d,", molstat.n_X);
  printf ("%d,", molstat.n_b1);
  printf ("%d,", molstat.n_b2);
  printf ("%d,", molstat.n_b3);
  printf ("%d,", molstat.n_bar);
  printf ("%d,", molstat.n_C1O);
  printf ("%d,", molstat.n_C2O);
  printf ("%d,", molstat.n_CN);
  printf ("%d,", molstat.n_XY);
  printf ("%d,", molstat.n_r3);
  printf ("%d,", molstat.n_r4);
  printf ("%d,", molstat.n_r5);
  printf ("%d,", molstat.n_r6);
  printf ("%d,", molstat.n_r7);
  printf ("%d,", molstat.n_r8);
  printf ("%d,", molstat.n_r9);
  printf ("%d,", molstat.n_r10);
  printf ("%d,", molstat.n_r11);
  printf ("%d,", molstat.n_r12);
  printf ("%d,", molstat.n_r13p);
  printf ("%d,", molstat.n_rN);
  printf ("%d,", molstat.n_rN1);
  printf ("%d,", molstat.n_rN2);
  printf ("%d,", molstat.n_rN3p);
  printf ("%d,", molstat.n_rO);
  printf ("%d,", molstat.n_rO1);
  printf ("%d,", molstat.n_rO2p);
  printf ("%d,", molstat.n_rS);
  printf ("%d,", molstat.n_rX);
  printf ("%d", molstat.n_rAr);
/* p2c: checkmol.pas, line 8579:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
  /*$IFDEF extended_molstat */
  printf (",%d", molstat.n_rBz);
  printf (",%d", molstat.n_br2p);
  printf (",%d", molstat.n_psg01);
  printf (",%d", molstat.n_psg02);
  printf (",%d", molstat.n_psg13);
  printf (",%d", molstat.n_psg14);
  printf (",%d", molstat.n_psg15);
  printf (",%d", molstat.n_psg16);
  printf (",%d", molstat.n_psg17);
  printf (",%d", molstat.n_psg18);
  printf (",%d", molstat.n_pstm);
  printf (",%d", molstat.n_psla);
  printf (",%d", molstat.n_iso);
  printf (",%d\n", molstat.n_rad);
  /*$ENDIF */
}


/* routines for substructure matching */


static int
find_ndl_ref_atom ()
{
  int Result, i;
  int score = -1, index = 0;
  int n_nb, n_hc, FORLIM;

  /* finds a characteristic atom in the needle molecule, */
  /* i.e., one with as many substituents as possible and */
  /* with as many heteroatom substitutents as possible; */
  /* added in v0.2d: make sure that reference atom is a heavy atom */
  /* and not (accidentally) an explicit hydrogen; */
  /* new in v0.3d: special treatment in case of E/Z geometry search */
  /* to ensure that the entire A-B=C-D fragment is enclosed in one */
  /* matchpath, regardless where the recursive search starts; */
  /* refined in v0.3f: exclude only alkene-C as reference atoms */
  /* added in v0.3o: needle atom must be "tagged" in order to be */
  /* selected (prevents unconnected fragments from being overlooked) */
  if (ndl_n_atoms == 0)
    return Result;
  if (ez_search && ndl_n_heavyatoms > 2)
    {
      FORLIM = ndl_n_atoms;
      for (i = 1; i <= FORLIM; i++)
	{			/* ignore sp2-carbons if not aromatic */
	  /*if ((ndl_atom^[i].atype <> 'C2 ') or (ndl_atom^[i].arom = true)) then */
	  if (ndl_alkene_C (i) == false && ndl_atom[i - 1].tag)
	    {			/* v0.3o */
	      n_nb = ndl_atom[i - 1].neighbor_count;
	      n_hc = ndl_hetatom_count (i);
	      if (n_nb * 11 + n_hc * 7 > score && ndl_atom[i - 1].heavy)
		{
		  /* v0.3j */
		  index = i;
		  score = n_nb * 11 + n_hc * 7;	/* changed in v0.3j */
		}
	    }
	}
    }
  /* it is possible that no suitable reference atom has been found here */
  /* (e.g., with "pure" polyenes), so we need a fallback option anyway */
  if (index == 0)
    {
      ez_search = false;	/* just in case it was true */
      opt_geom = false;		/* just in case it was true */
      FORLIM = ndl_n_atoms;
      for (i = 1; i <= FORLIM; i++)
	{
	  n_nb = ndl_atom[i - 1].neighbor_count;
	  n_hc = ndl_hetatom_count (i);
	  if (n_nb * 11 + n_hc * 7 > score && ndl_atom[i - 1].heavy &&
	      ndl_atom[i - 1].tag)
	    {			/* v0.3j */
	      index = i;
	      score = n_nb * 11 + n_hc * 7;	/* changed in v0.3j */
	    }
	  /* v0.3o */
	}
    }
  /* now index must be > 0 in any case (except for H2, or all tags have been cleared) */
  if (index == 0)		/* just to be sure... */
    index++;
  return index;
}


static void
cv_init ()
{
  /* new in v0.3j */
  int i;

  if (cv == NULL)
    return;
  memset (cv, 0, sizeof (connval_type));

  for (i = 0; i < ndl_n_atoms; i++)
    cv[i].def = ndl_atom[i].neighbor_count;
}


static int
cv_count ()
{
  /* new in v0.3j, modified in v0.3m */
  int i, j;
  int cvlist[max_atoms];
  int cvdef;
  boolean isnew;
  int entries = 0;
  int FORLIM;

  if (cv == NULL)
    return 0;
  memset (cvlist, 0, sizeof (int) * max_atoms);
  FORLIM = ndl_n_atoms;
  for (i = 0; i < FORLIM; i++)
    {
      if (ndl_atom[i].heavy == true)
	{
	  cvdef = cv[i].def;
	  isnew = true;
	  if (entries > 0)
	    {
	      for (j = 0; j < entries; j++)
		{
		  if (cvlist[j] == cvdef)
		    isnew = false;
		}
	    }
	  if (isnew)
	    {
	      entries++;
	      cvlist[entries - 1] = cvdef;
	    }
	  /* now we have a list of unique connection values */
	}
    }
  return entries;
}


static int
cv_iterate (n_cv_prev)
     int n_cv_prev;
{
  /* new in v0.3j, modified in v0.3m */
  int Result, i, j;
  neighbor_rec nb;
  int nnb, nsum, n_cv, FORLIM;

  if (cv == NULL || ndl_n_atoms == 0)
    return Result;
  FORLIM = ndl_n_atoms;
  /* update the connection values (Morgan algorithm) */

  memset (nb, 0, sizeof (neighbor_rec));

  for (i = 1; i <= FORLIM; i++)
    {
      if (ndl_atom[i - 1].heavy == true)
	{
	  get_ndl_neighbors (nb, i);
	  nnb = ndl_atom[i - 1].neighbor_count;
	  nsum = 0;
	  if (nnb > 0)
	    {
	      for (j = 0; j < nnb; j++)
		{
		  if (ndl_atom[nb[j] - 1].heavy == true)
		    nsum += cv[nb[j] - 1].def;
		}
	    }
	  cv[i - 1].tmp = nsum;
	}
    }
  n_cv = cv_count ();
  if (n_cv > n_cv_prev)
    {
      FORLIM = ndl_n_atoms;
      for (i = 0; i < FORLIM; i++)
	cv[i].def = cv[i].tmp;
    }
  return n_cv;
}


static int
find_ndl_ref_atom_cv ()
{
  /* new in v0.3j, modified in v0.3m */
  int Result, i;
  int res = 1, it = 0;
  int n_cv;
  int n_cv_prev = 0;
  boolean finished = false;
  int cvmax = 0;
  int FORLIM;

  if (ndl_n_atoms == 0)
    return 0;
  /* try */
  cv = (connval_rec *) safe_malloc (sizeof (connval_type));
  /* except
     on e:Eoutofmemory do
     begin
     res := find_ndl_ref_atom;
     $IFDEF debug
     debugoutput('memory allocation for connection values failed, reverting to standard procedure');
     $ENDIF
     end;
     end; */
  cv_init ();
  do
    {
      it++;			/* iteration counter (a safeguard against infinite loops) */
      n_cv = cv_iterate (n_cv_prev);
      if (n_cv <= n_cv_prev)
	finished = true;
      n_cv_prev = n_cv;
    }
  while (!(finished || it > 10000));
  FORLIM = ndl_n_atoms;
  /* now that we have canonical connection values (Morgan algorithm), */
  /* pick the atom with the highest value */
  /* added in v0.3o: atom must be "tagged" */
  for (i = 1; i <= FORLIM; i++)
    {
      /*writeln('cv for atom ',i,': ',cv^[i].def); */
      if (((cv[i - 1].def > cvmax) && (ndl_alkene_C (i) == false ||
				       ez_search == false))
	  && ndl_atom[i - 1].tag)
	{			/* v0.3o */
	  cvmax = cv[i - 1].def;
	  res = i;
	}
    }
  Result = res;
  /* try */
  if (cv != NULL)
    {
      free (cv);
      cv = NULL;
    }
  /* except
     on e:Einvalidpointer do begin end;
     end; */
  return Result;
}


static boolean
atomtypes_OK_strict (ndl_a, hst_a)
     int ndl_a, hst_a;
{
  /* new in v0.2f */
  str2 ndl_el;
  str3 ndl_atype;
  str2 hst_el;
  str3 hst_atype;
  int ndl_nbc, hst_nbc, ndl_Hexp, hst_Htot;
  boolean res = false;

  strcpy (ndl_el, ndl_atom[ndl_a - 1].element);
  strcpy (ndl_atype, ndl_atom[ndl_a - 1].atype);
  ndl_nbc = ndl_atom[ndl_a - 1].neighbor_count;
  ndl_Hexp = ndl_atom[ndl_a - 1].Hexp;
  strcpy (hst_el, atom[hst_a - 1].element);
  strcpy (hst_atype, atom[hst_a - 1].atype);
  hst_nbc = atom[hst_a - 1].neighbor_count;
  hst_Htot = atom[hst_a - 1].Htot;
  /* v0.3o: formal charges must be the same */

  if (ndl_atom[ndl_a - 1].formal_charge != atom[hst_a - 1].formal_charge)
    return false;

  /* v0.3x: isotope nucleon numbers must be the same */

  if (ndl_atom[ndl_a - 1].nucleon_number != atom[hst_a - 1].nucleon_number)
    return false;

  /* v0.3x: radicals must be the same */

  if (ndl_atom[ndl_a - 1].radical_type != atom[hst_a - 1].radical_type)
    return false;

  if (!strcmp (ndl_atype, hst_atype))
    res = true;
  else
    {
      if (!strcmp (ndl_el, hst_el) && ndl_atom[ndl_a - 1].arom &&
	  atom[hst_a - 1].arom)
	res = true;
      if (ndl_querymol
	  && (ndl_atom[ndl_a - 1].q_arom && atom[hst_a - 1].arom))
	res = true;		/* 0.3 p */
    }
  if (!strcmp (ndl_el, "A ") && atom[hst_a - 1].heavy)
    res = true;
  if (!strcmp (ndl_el, "Q "))
    {
      if (atom[hst_a - 1].heavy && strcmp (hst_el, "C "))
	res = true;
    }
  if (!strcmp (ndl_el, "X "))
    {
      if (!strcmp (hst_el, "F ") || !strcmp (hst_el, "CL") ||
	  !strcmp (hst_el, "BR") || !strcmp (hst_el, "I ")
	  || !strcmp (hst_el, "AT"))
	res = true;
    }
  /* if needle atom has more substituents than haystack atom ==> no match */
  if (ndl_nbc > hst_nbc)
    res = false;
  /* check for explicit hydrogens */
  if (ndl_Hexp > hst_Htot)
    res = false;
/* p2c: checkmol.pas, line 8859:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
  /*$IFDEF debug */
  /* if res then debugoutput('atom types OK ('+inttostr(ndl_a)+'/'+inttostr(hst_a)+')')
     else debugoutput('atom types not OK ('+inttostr(ndl_a)+':'+ndl_atype+'/'+inttostr(hst_a)+':'+hst_atype+')'); */
  /*$ENDIF */
  /* new in v0.3m: in "fingerprint mode", also query atom symbols must match */
  if (opt_fp)
    {
      if (strcmp (ndl_el, hst_el))
	res = false;
    }
  return res;
}


static boolean
atomtypes_OK (ndl_a, hst_a)
     int ndl_a, hst_a;
{
  str2 ndl_el, hst_el;
  int ndl_nbc, hst_nbc, ndl_Hexp, hst_Htot;
  boolean res = false;

  if (ndl_a < 1 || ndl_a > ndl_n_atoms || hst_a < 1 || hst_a > n_atoms)
    return false;
  /* check for opposite charges;  v0.3l, refined in v0.3o, 0.3x */
  /* except in strict mode, matching pairs of charged+uncharged atoms  */
  /* are tolerated (this is a feature, not a bug) */
  if (opt_chg)
    {
      if (ndl_atom[ndl_a - 1].formal_charge != atom[hst_a - 1].formal_charge)
	return false;
    }
//  else
//    {
//      if (ndl_atom[ndl_a - 1].formal_charge != 0 &&
//        atom[hst_a - 1].formal_charge != 0 &&
//        ndl_atom[ndl_a - 1].formal_charge != atom[hst_a - 1].formal_charge)
//      return false;
//    }
//
//  /* v0.3x: isotopes must be the same */
  if (opt_iso)
    {
      if (ndl_atom[ndl_a - 1].nucleon_number !=
	  atom[hst_a - 1].nucleon_number)
	return false;
    }
//  else
//    {
//      if (ndl_atom[ndl_a - 1].nucleon_number != 0 &&
//        atom[hst_a - 1].nucleon_number != 0 &&
//        ndl_atom[ndl_a - 1].nucleon_number !=
//        atom[hst_a - 1].nucleon_number)
//      return false;
//    }
//
//  /* v0.3x: radicals must be the same */
  if (opt_rad)
    {
      if (ndl_atom[ndl_a - 1].radical_type != atom[hst_a - 1].radical_type)
	return false;
    }
//  else
//    {
//      if (ndl_atom[ndl_a - 1].radical_type != 0 &&
//        atom[hst_a - 1].radical_type != 0 &&
//        ndl_atom[ndl_a - 1].radical_type != atom[hst_a - 1].radical_type)
//      return false;
//    }

  /* in exact mode, check if (disconnected) fragment is already tagged; v0.3o */
  if (opt_exact && atom[hst_a - 1].tag == true)
    {
/* p2c: checkmol.pas, line 8899:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
      /*$IFDEF debug */
      /* debugoutput('fragmnet already tagged at '+inttostr(hst_a)); */
      /*$ENDIF */
      return false;
    }
  if (opt_strict)		/* new in v0.2f */
    return (atomtypes_OK_strict (ndl_a, hst_a));
  strcpy (ndl_el, ndl_atom[ndl_a - 1].element);
  ndl_nbc = ndl_atom[ndl_a - 1].neighbor_count;
  ndl_Hexp = ndl_atom[ndl_a - 1].Hexp;
  strcpy (hst_el, atom[hst_a - 1].element);
  hst_nbc = atom[hst_a - 1].neighbor_count;
  hst_Htot = atom[hst_a - 1].Htot;
  if (!strcmp (ndl_el, hst_el))	/* very simplified... */
    res = true;
  if (!strcmp (ndl_el, "A ") && atom[hst_a - 1].heavy)
    res = true;
  if (!strcmp (ndl_el, "Q "))
    {
      if (atom[hst_a - 1].heavy && strcmp (hst_el, "C "))
	res = true;
    }
  if (!strcmp (ndl_el, "X "))
    {
      if (!strcmp (hst_el, "F ") || !strcmp (hst_el, "CL") ||
	  !strcmp (hst_el, "BR") || !strcmp (hst_el, "I ")
	  || !strcmp (hst_el, "AT"))
	res = true;
    }
  /* v0.3o: in exact mode, check for identical neighbor_count */
  if (opt_exact)
    {
      if (ndl_nbc != hst_nbc)
	{
	  res = false;
/* p2c: checkmol.pas, line 8934:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
	  /*$IFDEF debug */
	  //debugoutput
	  //  ("exact match failed: different number of neighbor atoms");
	  /*$ENDIF */
	}
    }
  /* if needle atom has more substituents than haystack atom ==> no match */
  if (ndl_nbc > hst_nbc)
    res = false;
  /* check for explicit hydrogens */
  if (ndl_Hexp > hst_Htot)
    res = false;
/* p2c: checkmol.pas, line 8943:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
  /*$IFDEF debug */
  /* if res then debugoutput('atom types OK ('+inttostr(ndl_a)+'/'+inttostr(hst_a)+')')
     else debugoutput('atom types not OK ('+inttostr(ndl_a)+'/'+inttostr(hst_a)+')'); */
  /*$ENDIF */
  return res;
}


static boolean
bondtypes_OK_strict (ndl_b, hst_b)
     int ndl_b, hst_b;
{
  boolean ndl_arom, hst_arom;
  char ndl_btype, hst_btype;
  int ndl_rc;			/* new in v0.3d */
  int hst_rc;			/* new in v0.3d */
  int ndl_btopo;		/* new in v0.3d */
  boolean res = false;
/* p2c: checkmol.pas, line 8960:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
  /*$IFDEF debug */
  char na[256], ha[256];
  char tstr[256];

  /*$ENDIF */
/* p2c: checkmol.pas, line 8966:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
  /*$IFDEF debug */
  *tstr = '\0';			/* for debugging purposes only */
  /*$ENDIF */
  ndl_arom = ndl_bond[ndl_b - 1].arom;
  ndl_btype = ndl_bond[ndl_b - 1].btype;
  ndl_rc = ndl_bond[ndl_b - 1].ring_count;
  ndl_btopo = ndl_bond[ndl_b - 1].topo;
  hst_arom = bond[hst_b - 1].arom;
  hst_btype = bond[hst_b - 1].btype;
  hst_rc = bond[hst_b - 1].ring_count;
/* p2c: checkmol.pas, line 8976:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
  /*$IFDEF debug */
  /*if (ndl_arom)
    strcpy (na, "(ar)");
  else
    *na = '\0';
  if (hst_arom)
    strcpy (ha, "(ar)");
  else
    *ha = '\0';*/
  /*$ENDIF */
  if (ndl_arom == true && hst_arom == true)
    res = true;
  if (ndl_arom == false && hst_arom == false)
    {
      if (ndl_btype == hst_btype)
	res = true;
      if (ndl_btype == 'l' && (hst_btype == 'S' || hst_btype == 'D'))
	res = true;
      if (ndl_btype == 's' && hst_btype == 'S')
	res = true;
      if (ndl_btype == 'd' && hst_btype == 'D')
	res = true;
    }
  /* a little exception: */
  if (ndl_arom == false && hst_arom == true)
    {
      if (ndl_btype == 'A')
	res = true;
      if (ndl_btype == 's' || ndl_btype == 'd')
	res = true;
      if (ndl_bond[ndl_b - 1].q_arom)
	res = true;		/* 0.3p */
    }
  if (ndl_btype == 'a')
    res = true;
  /* new in v0.3d: strict comparison of topology (and even ring_count!) */
  if (ndl_btopo < btopo_always_any || ndl_btopo == btopo_exact_rc)
    {
      if (ndl_rc != hst_rc)
	{
	  res = false;		/* this excludes further ring annulations as well as */
/* p2c: checkmol.pas, line 9001:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
	  /*$IFDEF debug */
	  /* open-chains query structures to be found in rings */
	  /*
	     tstr := ' ringcount mismatch ('+inttostr(ndl_rc)+'/'+inttostr(hst_rc)+')';   */
	  /*$ENDIF */
	}
    }
  else
    {
      if (ndl_btopo == btopo_excess_rc && hst_rc <= ndl_rc)
	{
	  res = false;
/* p2c: checkmol.pas, line 9010:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
	  /*$IFDEF debug */
	  /* tstr := ' ringcount mismatch ('+inttostr(ndl_rc)+'/'+inttostr(hst_rc)+')'; */
	  /*$ENDIF */
	}
    }
/* p2c: checkmol.pas, line 9015:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
  /*$IFDEF debug */
  /* if res then debugoutput('bond types OK ('+inttostr(ndl_b)+':'+ndl_btype+na+'/'+inttostr(hst_b)+':'+hst_btype+ha+')') else
     debugoutput('bond types not OK ('+inttostr(ndl_b)+':'+ndl_btype+na+'/'+inttostr(hst_b)+':'+hst_btype+ha+tstr+')'); */
  /*$ENDIF */
  return res;
}


static boolean
bondtypes_OK (ndl_b, hst_b)
     int ndl_b, hst_b;
{
  boolean ndl_arom, hst_arom;
  char ndl_btype, hst_btype;
  int ndl_rc;			/* new in v0.3d */
  int hst_rc;			/* new in v0.3d */
  int ndl_btopo;		/* new in v0.3d */
  boolean res = false;
/* p2c: checkmol.pas, line 9032:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
  /*$IFDEF debug */
  char na[256], ha[256];
  char tstr[256];
  /*$ENDIF */
  int a1, a2;
  str2 a1_el, a2_el;

  if (ndl_b < 1 || ndl_b > ndl_n_bonds || hst_b < 1 || hst_b > n_bonds)
    return false;
  if (opt_strict)		/* new in v0.2f */
    return (bondtypes_OK_strict (ndl_b, hst_b));
/* p2c: checkmol.pas, line 9051:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
  /*$IFDEF debug */
  *tstr = '\0';			/* for debug purposes only */
  /*$ENDIF */
  ndl_arom = ndl_bond[ndl_b - 1].arom;
  ndl_btype = ndl_bond[ndl_b - 1].btype;
  hst_arom = bond[hst_b - 1].arom;
  hst_btype = bond[hst_b - 1].btype;
  ndl_rc = ndl_bond[ndl_b - 1].ring_count;
  hst_rc = bond[hst_b - 1].ring_count;
  ndl_btopo = ndl_bond[ndl_b - 1].topo;
/* p2c: checkmol.pas, line 9061:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
  /*$IFDEF debug */
  //if (ndl_arom)
//    strcpy (na, "(ar)");
//  else
//    *na = '\0';
//  if (hst_arom)
//    strcpy (ha, "(ar)");
//  else
//    *ha = '\0';
  /*$ENDIF */
  if (ndl_arom == true && hst_arom == true)
    res = true;
  if (ndl_arom == false && hst_arom == false)
    {
      if (ndl_btype == hst_btype)
	res = true;
      if (ndl_btype == 'l' && (hst_btype == 'S' || hst_btype == 'D'))
	res = true;
      if (ndl_btype == 's' && hst_btype == 'S')
	res = true;
      if (ndl_btype == 'd' && hst_btype == 'D')
	res = true;
    }
  /* a little exception: */
  if (ndl_arom == false && hst_arom == true)
    {
      if (ndl_btype == 'A')
	res = true;
      if (ndl_btype == 's' || ndl_btype == 'd')
	res = true;
      if (ndl_btype == 'D')
	{			/* added in 0.2d: do not accept C=O etc. as C-O/arom */
	  a1 = ndl_bond[ndl_b - 1].a1;
	  a2 = ndl_bond[ndl_b - 1].a2;
	  strcpy (a1_el, ndl_atom[a1 - 1].element);
	  strcpy (a2_el, ndl_atom[a2 - 1].element);
	  if (strcmp (a1_el, "O ") && strcmp (a2_el, "O ")
	      && strcmp (a1_el, "S ") && strcmp (a2_el, "S ")
	      && strcmp (a1_el, "SE") && strcmp (a2_el, "SE")
	      && strcmp (a1_el, "TE") && strcmp (a2_el, "TE"))
	    res = true;
	}
      if (ndl_bond[ndl_b - 1].q_arom)
	res = true;		/* 0.3p */
    }
  if (ndl_btype == 'a')
    res = true;
  /* new in v0.3d: obey topology requirements in query structure */
  if (ndl_btopo != btopo_any && ndl_btopo != btopo_always_any)
    {
  if (ndl_btopo == btopo_ring && hst_rc == 0)
    res = false;
  if (ndl_btopo == btopo_chain && hst_rc > 0)
    res = false;
  if (ndl_btopo == btopo_excess_rc && hst_rc <= ndl_rc)
    res = false;
  if (ndl_btopo == btopo_exact_rc && hst_rc != ndl_rc)
    res = false;
}    
/* p2c: checkmol.pas, line 9098:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
  /*$IFDEF debug */
  /* if res = false then tstr := ' bond topology mismatch '+inttostr(ndl_rc)+'/'+inttostr(hst_rc); */
  /*$ENDIF */
/* p2c: checkmol.pas, line 9102:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
  /*$IFDEF debug */
  /*
     if res then debugoutput('bond types OK ('+inttostr(ndl_b)+':'+ndl_btype+na+'/'+inttostr(hst_b)+':'+hst_btype+ha+')') else
     debugoutput('bond types not OK ('+inttostr(ndl_b)+':'+ndl_btype+na+'/'+inttostr(hst_b)+':'+hst_btype+ha+tstr+')'); */
  /*$ENDIF */
  return res;
}


static boolean
matrix_OK (m, ndl_dim, hst_dim)
boolean (*m)[max_neighbors];
     int ndl_dim, hst_dim;
{
  /* new, recursive version in v0.2i: can handle up to max_neighbors substituents */
  boolean mr = false;
  matchmatrix lm;
  int i, ii, j, lndl_dim, lhst_dim;

  if (ndl_dim < 1 || ndl_dim > max_neighbors || hst_dim < 1 ||
      hst_dim > max_neighbors || ndl_dim > hst_dim)
    return false;
  if (ndl_dim == 1)
    {
      for (i = 0; i < hst_dim; i++)
	{
	  if (m[0][i])
	    mr = true;
	}
      return mr;
    }
  for (i = 1; i <= hst_dim; i++)
    {
      if (m[0][i - 1])
	{
	  /* write remaining fields into a new matchmatrix which is smaller by 1x1 */
	  memset (lm, false, sizeof (matchmatrix));
	  for (j = 2; j <= ndl_dim; j++)
	    {
	      lhst_dim = 0;
	      for (ii = 1; ii <= hst_dim; ii++)
		{
		  if (ii != i)
		    {
		      lhst_dim++;
		      lm[j - 2][lhst_dim - 1] = m[j - 1][ii - 1];
		    }
		}
	    }
	  lndl_dim = ndl_dim - 1;
	  if (matrix_OK (lm, lndl_dim, lhst_dim))
	    {			/* recursive call to matrix_OK */
	      return true;
	      /* stop any further work immediately */
	    }
	}
    }
  return false;
}


static boolean
is_flat (angle_deg)
     double angle_deg;
{
  /* new in v0.3j */
  if (fabs (angle_deg) > 5 && fabs (angle_deg) < 175)
    return false;
  else
    return true;
}


static boolean
chirality_OK (ndl_cp, hst_cp)
     int *ndl_cp, *hst_cp;
{
  boolean res = true;
  double ndl_ct, hst_ct, ndl_ct_deg, hst_ct_deg;
  p_3d np1, np2, np3, np4, hp1, hp2, hp3, hp4;
  int level = 0;
  int i;
  boolean up = false, down = false, updown = false;
  int ta1, ta2, ta3, ta4, ba1, ba2, FORLIM;

  /* fill temporary atom variables */
  ta1 = ndl_cp[0];		/* this is the central atom */
  ta2 = ndl_cp[1];
  ta3 = ndl_cp[2];
  ta4 = ndl_cp[3];
  /* first, get the central atom of the needle */
  np2.x = ndl_atom[ta1 - 1].x;
  np2.y = ndl_atom[ta1 - 1].y;
  np2.z = ndl_atom[ta1 - 1].z;
  /* next, do the same for all 3 substituent atoms */
  np1.x = ndl_atom[ta2 - 1].x;
  np1.y = ndl_atom[ta2 - 1].y;
  np1.z = ndl_atom[ta2 - 1].z;
  np3.x = ndl_atom[ta3 - 1].x;
  np3.y = ndl_atom[ta3 - 1].y;
  np3.z = ndl_atom[ta3 - 1].z;
  np4.x = ndl_atom[ta4 - 1].x;
  np4.y = ndl_atom[ta4 - 1].y;
  np4.z = ndl_atom[ta4 - 1].z;
  /* now check all needle bonds if we should care about up/down bonds */
  if (ndl_n_bonds > 0)
    {
      FORLIM = ndl_n_bonds;
      for (i = 0; i < FORLIM; i++)
	{
	  if (ndl_bond[i].stereo == bstereo_up ||
	      ndl_bond[i].stereo == bstereo_down)
	    {
	      ba1 = ndl_bond[i].a1;
	      ba2 = ndl_bond[i].a2;
	      if (ba1 == ta1 && ndl_bond[i].stereo == bstereo_up)
		{
		  up = true;
		  if (ba2 == ta2 || ba2 == ta3 || ba2 == ta4)
		    {
		      updown = true;
		      if (ba2 == ta2)
			np1.z += 0.8;
		      if (ba2 == ta3)
			np3.z += 0.8;
		      if (ba2 == ta4)
			np4.z += 0.8;
		    }
		  else
		    level++;
		}
	      if (ba1 == ta1 && ndl_bond[i].stereo == bstereo_down)
		{
		  down = true;
		  if (ba2 == ta2 || ba2 == ta3 || ba2 == ta4)
		    {
		      updown = true;
		      if (ba2 == ta2)
			np1.z -= 0.8;
		      if (ba2 == ta3)
			np3.z -= 0.8;
		      if (ba2 == ta4)
			np4.z -= 0.8;
		    }
		  else
		    level--;
		}
	      if (ba2 == ta1 && ndl_bond[i].stereo == bstereo_up)
		{
		  down = true;
		  if (ba1 == ta2 || ba1 == ta3 || ba1 == ta4)
		    {
		      updown = true;
		      if (ba1 == ta2)
			np1.z -= 0.8;
		      if (ba1 == ta3)
			np3.z -= 0.8;
		      if (ba1 == ta4)
			np4.z -= 0.8;
		    }
		  else
		    level--;
		}
	      if (ba2 == ta1 && ndl_bond[i].stereo == bstereo_down)
		{
		  up = true;
		  if (ba1 == ta2 || ba1 == ta3 || ba1 == ta4)
		    {
		      updown = true;
		      if (ba1 == ta2)
			np1.z += 0.8;
		      if (ba1 == ta3)
			np3.z += 0.8;
		      if (ba1 == ta4)
			np4.z += 0.8;
		    }
		  else
		    level++;
		}
	    }
	}			/* for i ... */
      if (updown == false && level != 0)
	{
	  if (level > 0)
	    np2.z += 0.3;
	  if (level < 0)
	    np2.z -= 0.3;
	}
      else
	{
	  if (up)
	    np2.z += 0.1;
	  if (down)
	    np2.z -= 0.1;
	}
    }
  /* fill temporary atom variables again */
  ta1 = hst_cp[0];
  ta2 = hst_cp[1];
  ta3 = hst_cp[2];
  ta4 = hst_cp[3];
  /* then, get the central atom of the haystack */
  hp2.x = atom[ta1 - 1].x;
  hp2.y = atom[ta1 - 1].y;
  hp2.z = atom[ta1 - 1].z;
  /* next, do the same for all 3 substituent atoms */
  hp1.x = atom[ta2 - 1].x;
  hp1.y = atom[ta2 - 1].y;
  hp1.z = atom[ta2 - 1].z;
  hp3.x = atom[ta3 - 1].x;
  hp3.y = atom[ta3 - 1].y;
  hp3.z = atom[ta3 - 1].z;
  hp4.x = atom[ta4 - 1].x;
  hp4.y = atom[ta4 - 1].y;
  hp4.z = atom[ta4 - 1].z;
  /* now check all haystack bonds if we should care about up/down bonds */
  level = 0;
  updown = false;
  up = false;
  down = false;
  if (n_bonds > 0)
    {
      FORLIM = n_bonds;
      for (i = 0; i < FORLIM; i++)
	{
	  if (bond[i].stereo == bstereo_up || bond[i].stereo == bstereo_down)
	    {
	      ba1 = bond[i].a1;
	      ba2 = bond[i].a2;
	      if (ba1 == ta1 && bond[i].stereo == bstereo_up)
		{
		  up = true;
		  if (ba2 == ta2 || ba2 == ta3 || ba2 == ta4)
		    {
		      updown = true;
		      if (ba2 == ta2)
			hp1.z += 0.8;
		      if (ba2 == ta3)
			hp3.z += 0.8;
		      if (ba2 == ta4)
			hp4.z += 0.8;
		    }
		  else
		    level++;
		}
	      if (ba1 == ta1 && bond[i].stereo == bstereo_down)
		{
		  down = true;
		  if (ba2 == ta2 || ba2 == ta3 || ba2 == ta4)
		    {
		      updown = true;
		      if (ba2 == ta2)
			hp1.z -= 0.8;
		      if (ba2 == ta3)
			hp3.z -= 0.8;
		      if (ba2 == ta4)
			hp4.z -= 0.8;
		    }
		  else
		    level--;
		}
	      if (ba2 == ta1 && bond[i].stereo == bstereo_up)
		{
		  down = true;
		  if (ba1 == ta2 || ba1 == ta3 || ba1 == ta4)
		    {
		      updown = true;
		      if (ba1 == ta2)
			hp1.z -= 0.8;
		      if (ba1 == ta3)
			hp3.z -= 0.8;
		      if (ba1 == ta4)
			hp4.z -= 0.8;
		    }
		  else
		    level--;
		}
	      if (ba2 == ta1 && bond[i].stereo == bstereo_down)
		{
		  up = true;
		  if (ba1 == ta2 || ba1 == ta3 || ba1 == ta4)
		    {
		      updown = true;
		      if (ba1 == ta2)
			hp1.z += 0.8;
		      if (ba1 == ta3)
			hp3.z += 0.8;
		      if (ba1 == ta4)
			hp4.z += 0.8;
		    }
		  else
		    level++;
		}
	    }
	}			/* for i ... */
      if (updown == false && level != 0)
	{
	  if (level > 0)
	    hp2.z += 0.3;
	  if (level < 0)
	    hp2.z -= 0.3;
	}
      else
	{
	  if (up)
	    hp2.z += 0.1;
	  if (down)
	    hp2.z -= 0.1;
	}
    }
  /* get the pseudo-torsion angles */
  ndl_ct = ctorsion (np1, np2, np3, np4);
  hst_ct = ctorsion (hp1, hp2, hp3, hp4);
  ndl_ct_deg = radtodeg (ndl_ct);
  hst_ct_deg = radtodeg (hst_ct);
  /* now do a plausibility check and finally check the sense */
  /* (clockwise or counterclockwise) */
  /*
     if (abs(ndl_ct_deg) > 5) and (abs(ndl_ct_deg) < 175) and
     (abs(hst_ct_deg) > 5) and (abs(hst_ct_deg) < 175) and
     (ndl_ct_deg * hst_ct_deg < 0) then res := false;
   */
  if (((!is_flat (ndl_ct_deg)) && (!is_flat (hst_ct_deg))) &&
      ndl_ct_deg * hst_ct_deg < 0)
    res = false;
  if (rs_strict)
    {
      if (((is_flat (ndl_ct_deg) && (!is_flat (hst_ct_deg))) |
	   (is_flat (hst_ct_deg) && (!is_flat (ndl_ct_deg)))) ||
	  ndl_ct_deg * hst_ct_deg < 0)
	res = false;
    }
  return res;
}


static boolean
ndl_maybe_chiral (na)
     int na;
{
  /* new in v0.3h */
  boolean res = false;
  str2 el;
  str3 at;
  int n_nb;

  strcpy (el, ndl_atom[na - 1].element);
  strcpy (at, ndl_atom[na - 1].atype);
  n_nb = ndl_atom[na - 1].neighbor_count;
  if (!strcmp (at, "C3 ") && n_nb > 2)
    res = true;
  if (!strcmp (el, "N "))
    {
      if (!strcmp (at, "N3+") && n_nb == 4)
	res = true;
    }
  if (!strcmp (el, "S "))
    {				/* sulfoxide */
      if ((n_nb == 3) && (ndl_hetatom_count (na) == 1))
	res = true;
    }
  if (strcmp (el, "P ") && strcmp (el, "AS"))	/* "As" added in v0.3j */
    return res;
  if (n_nb > 3)			/* are we missing something here? */
    res = true;
  if (ndl_hetatom_count (na) >= 2)	/* v0.3m; ignore phosphates etc. */
    res = false;
  return res;
}


static boolean
is_matching (ndl_xmp, hst_xmp)
     int *ndl_xmp, *hst_xmp;
{
  boolean Result;
  int i, j, k, l, m, ndl_n_nb, n_nb, ndl_a, hst_a;
  int ndl_b = 0, hst_b = 0, prev_ndl_a = 0, prev_hst_a = 0;
  int next_ndl_a, next_hst_a;
  neighbor_rec ndl_nb, hst_nb;
  matchmatrix mm;
  int ndl_mp_len, hst_mp_len;
  matchpath_type ndl_mp, hst_mp;
  boolean emptyline, res, ndl_cis, hst_cis;
  int na1, na2, na3, na4;	/* v0.3d */
  int ha1, ha2, ha3, ha4;	/* atom variables for E/Z check */
  int prev_ndl_b;
  int prev_hst_b;
  p_3d p1, p2, p3, p4;
  /*hst_torsion, ndl_torsion : double; */
  chirpath_type ncp, hcp;
  int n_hits, n_singlehits;
/* p2c: checkmol.pas, line 9433:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
  /*$IFDEF debug */
  //char tmpstr[256];

  /*$ENDIF */
  /* initialize local matchpath variables */
  //memset (ndl_mp, 0, sizeof (matchpath_type));
  //memset (hst_mp, 0, sizeof (matchpath_type));
  /* copy content of external variables into local ones */
  memcpy (ndl_mp, ndl_xmp, sizeof (matchpath_type));
  memcpy (hst_mp, hst_xmp, sizeof (matchpath_type));

  /*for (i = 0; i < max_matchpath_length; i++)
     {
     ndl_mp[i] = ndl_xmp[i];
     hst_mp[i] = hst_xmp[i];
     } */


  ndl_mp_len = matchpath_length (ndl_mp);
  hst_mp_len = matchpath_length (hst_mp);
  if (ndl_mp_len != hst_mp_len)
    {
      /* this should never happen.... */
/* p2c: checkmol.pas, line 9451:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
      /*$IFDEF debug */
      //debugoutput ("needle and haystack matchpaths are of different length");
      /*$ENDIF */
      return false;
    }
  ndl_a = ndl_mp[ndl_mp_len - 1];
  hst_a = hst_mp[hst_mp_len - 1];
  ndl_atom[ndl_a - 1].tag = false;
  /* new in v0.3o: mark the last needle atom as "visited" */
  if (ndl_mp_len > 1)
    {
      prev_ndl_a = ndl_mp[ndl_mp_len - 2];
      prev_hst_a = hst_mp[hst_mp_len - 2];
    }
  /* if geometry checking is on, check it here */
  if (ez_search == true && ndl_mp_len > 3)
    {
      na1 = ndl_mp[ndl_mp_len - 1];
      na2 = ndl_mp[ndl_mp_len - 2];
      na3 = ndl_mp[ndl_mp_len - 3];
      na4 = ndl_mp[ndl_mp_len - 4];
      ha1 = hst_mp[hst_mp_len - 1];
      ha2 = hst_mp[hst_mp_len - 2];
      ha3 = hst_mp[hst_mp_len - 3];
      ha4 = hst_mp[hst_mp_len - 4];
      prev_ndl_b = get_ndl_bond (na2, na3);
      prev_hst_b = get_bond (ha2, ha3);
      if (ndl_bond[prev_ndl_b - 1].btype == 'D' &&
	  bond[prev_hst_b - 1].arom == false
	  && (ndl_bond[prev_ndl_b - 1].stereo !=
	      bstereo_double_either && bond[prev_hst_b - 1].stereo !=
	      bstereo_double_either)
	  /* 0.3x always match if needle and/or haystack bond is double_either */
	  &&
	  (!strcmp (atom[ha2 - 1].element, "C ")
	   || !strcmp (atom[ha2 - 1].element, "N "))
	  && (!strcmp (atom[ha3 - 1].element, "C ")
	      || !strcmp (atom[ha3 - 1].element, "N ")))
	{			/* v0.3g; check C=C, C=N, N=N bonds */
	  p1.x = atom[ha1 - 1].x;
	  p1.y = atom[ha1 - 1].y;
	  p1.z = atom[ha1 - 1].z;
	  p2.x = atom[ha2 - 1].x;
	  p2.y = atom[ha2 - 1].y;
	  p2.z = atom[ha2 - 1].z;
	  p3.x = atom[ha3 - 1].x;
	  p3.y = atom[ha3 - 1].y;
	  p3.z = atom[ha3 - 1].z;
	  p4.x = atom[ha4 - 1].x;
	  p4.y = atom[ha4 - 1].y;
	  p4.z = atom[ha4 - 1].z;
	  hst_cis = is_cis (p1, p2, p3, p4);
	  /*hst_torsion := torsion(p1,p2,p3,p4); */
	  p1.x = ndl_atom[na1 - 1].x;
	  p1.y = ndl_atom[na1 - 1].y;
	  p1.z = ndl_atom[na1 - 1].z;
	  p2.x = ndl_atom[na2 - 1].x;
	  p2.y = ndl_atom[na2 - 1].y;
	  p2.z = ndl_atom[na2 - 1].z;
	  p3.x = ndl_atom[na3 - 1].x;
	  p3.y = ndl_atom[na3 - 1].y;
	  p3.z = ndl_atom[na3 - 1].z;
	  p4.x = ndl_atom[na4 - 1].x;
	  p4.y = ndl_atom[na4 - 1].y;
	  p4.z = ndl_atom[na4 - 1].z;
	  /*ndl_torsion := torsion(p1,p2,p3,p4); */
	  ndl_cis = is_cis (p1, p2, p3, p4);
	  if (ndl_cis != hst_cis)
	    {
/* p2c: checkmol.pas, line 9501:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
	      /*$IFDEF debug */
	      //debugoutput ("E/Z geometry mismatch");
	      /*$ENDIF */
	      return false;
	    }
	}
    }				/* end of E/Z geometry check */
  /* check whatever can be checked as early as now: */
  /* e.g. different elements or more substituents on needle atom than on haystack */
  if (!atomtypes_OK (ndl_a, hst_a))
    return false;
  /* positive scenarios, e.g. one-atom fragments  (v0.3o) */
  if (atom[hst_a - 1].neighbor_count == 0 &&
      ndl_atom[ndl_a - 1].neighbor_count == 0)
    {
      if (!atomtypes_OK (ndl_a, hst_a))
	return false;
/* p2c: checkmol.pas: Note: Eliminated unused assignment statement [338] */
      atom[hst_a - 1].tag = true;
      return true;
    }
  /* and other possibilities: */
  ndl_b = get_ndl_bond (prev_ndl_a, ndl_a);
  hst_b = get_bond (prev_hst_a, hst_a);
/* p2c: checkmol.pas, line 9529:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
  /*$IFDEF debug */
  /* debugoutput('Now checking atoms '+inttostr(ndl_a)+'/'+inttostr(hst_a)+', bonds '+inttostr(ndl_b)+'/'+inttostr(hst_b)); */
  /*$ENDIF */
  if (ndl_b > 0 && hst_b > 0)
    {
      /* do a quick check if bond types match */
      if (!bondtypes_OK (ndl_b, hst_b))
	{
/* p2c: checkmol.pas, line 9537:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
	  /*$IFDEF debug */
	  /*
	     debugoutput('  failed match of bonds '+inttostr(ndl_b)+'/'+inttostr(hst_b)); */
	  /*$ENDIF */
	  return false;
	}
    }
  /* a) we reached the end of our needle fragment (and atom/bond types match) */
  if ((ndl_atom[ndl_a - 1].neighbor_count == 1) && atomtypes_OK (ndl_a,
								 hst_a) &&
      bondtypes_OK (ndl_b, hst_b))
    {
      return true;
/* p2c: checkmol.pas, line 9549:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
      /*$IFDEF debug */
      /* debugoutput('  ==> end of needle fragment at atom '+inttostr(ndl_a)+' (match)'); */
      /*$ENDIF */
    }
  /* a.1) haystack fragment forms a ring, but needle does not;  v0.3m */
  if ((matchpath_pos (ndl_a, ndl_mp) == matchpath_length (ndl_mp)) &&
      (matchpath_pos (hst_a, hst_mp) < matchpath_length (hst_mp)))
    {
      return false;
/* p2c: checkmol.pas, line 9559:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
      /*$IFDEF debug */
      /*
         debugoutput('  haystack forms a ring and needle does not at '+inttostr(hst_a));
         {$ENDIF */
    }
  /* b) a ring is formed (ndl_a is already in the path) and atom/bond types match */
  if ((matchpath_pos (ndl_a, ndl_mp) > 0) &&
      (matchpath_pos (ndl_a, ndl_mp) < matchpath_length (ndl_mp)))
    {
      if ((matchpath_pos (ndl_a, ndl_mp) == matchpath_pos (hst_a, hst_mp)) &&
	  atomtypes_OK (ndl_a, hst_a) && bondtypes_OK (ndl_b, hst_b))
	{
	  /* 1st chirality check */
	  if (!((matchpath_pos (ndl_a, ndl_mp) > 1 && (rs_search ||
						       ndl_atom[ndl_a -
								1].
						       stereo_care)) &&
		ndl_maybe_chiral (ndl_a)))
	    {			/* new in v0.3h */
	      return true;
	    }			/* end of 1st chirality check */
	  na1 = ndl_a;		/* the (potential) chiral center (v0.3f) */
	  na2 = ndl_mp[matchpath_pos (ndl_a, ndl_mp) - 2];
	  na3 = ndl_mp[matchpath_pos (ndl_a, ndl_mp)];
	  na4 = ndl_mp[matchpath_length (ndl_mp) - 2];
	  ha1 = hst_a;
	  ha2 = hst_mp[matchpath_pos (hst_a, hst_mp) - 2];
	  ha3 = hst_mp[matchpath_pos (hst_a, hst_mp)];
	  ha4 = hst_mp[matchpath_length (hst_mp) - 2];
	  memset (ncp, 0, sizeof (chirpath_type));
	  memset (hcp, 0, sizeof (chirpath_type));
	  ncp[0] = na1;
	  ncp[1] = na2;
	  ncp[2] = na3;
	  ncp[3] = na4;
	  hcp[0] = ha1;
	  hcp[1] = ha2;
	  hcp[2] = ha3;
	  hcp[3] = ha4;
	  if (!chirality_OK (ncp, hcp))
	    {
/* p2c: checkmol.pas, line 9589:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
	      /*$IFDEF debug */
	      //debugoutput ("chirality check failed at ring junction");
	      /*$ENDIF */
	      return false;
	    }
/* p2c: checkmol.pas, line 9596:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
	  /*$IFDEF debug */
	  //debugoutput ("chirality check succeeded at ring junction");
	  /*$ENDIF */
	  return true;
/* p2c: checkmol.pas, line 9602:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
	  /*$IFDEF debug */
	  /* debugoutput('matchpath forms ring at: '+inttostr(ndl_a)+' (match)'); */
	  /*$ENDIF */
	}
      else
	{
	  return false;
/* p2c: checkmol.pas, line 9609:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
	  /*$IFDEF debug */
	  /*
	     debugoutput('matchpath forms ring at: '+inttostr(ndl_a)+' (no match)'); */
	  /*$ENDIF */
	}
    }
  /* in all other cases, do the hard work: */
  /* first, get all heavy-atom neighbors of needle and haystack; */
  /* at the beginning of the search, this means all neighbors, then it means */
  /* all but the previous atom (where we came from) */
  memset (ndl_nb, 0, sizeof (neighbor_rec));
  memset (hst_nb, 0, sizeof (neighbor_rec));

  if (matchpath_length (ndl_mp) == 1)
    {
      ndl_n_nb = ndl_atom[ndl_a - 1].neighbor_count;
      n_nb = atom[hst_a - 1].neighbor_count;
      get_ndl_neighbors (ndl_nb, ndl_a);
      get_neighbors (hst_nb, hst_a);
    }
  else
    {
      ndl_n_nb = ndl_atom[ndl_a - 1].neighbor_count - 1;
      n_nb = atom[hst_a - 1].neighbor_count - 1;
      get_ndl_nextneighbors (ndl_nb, ndl_a, prev_ndl_a);
      get_nextneighbors (hst_nb, hst_a, prev_hst_a);
    }
  /* v0.3o: mark all neighbor atoms as "visited" */
  for (i = 0; i < ndl_n_nb; i++)
    ndl_atom[ndl_nb[i] - 1].tag = false;
  /* now that the neighbor-arrays are filled, get all */
  /* combinations of matches recursively; */
  /* first, initialize the match matrix */
  memset (mm, false, sizeof (matchmatrix));	/* new in v0.2i */
  /* make sure there are not too many neighbors (max. max_neighbors)   */
  if (ndl_n_nb > max_neighbors || n_nb > max_neighbors)
    {				/* updated in v0.2i */
/* p2c: checkmol.pas, line 9644:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
      /*$IFDEF debug */
      //debugoutput ("too many neighbors - exiting");
      /*$ENDIF */
      return false;
    }
  /* check if matchpath is not already filled up */
  if (matchpath_length (ndl_mp) == max_matchpath_length)
    {
/* p2c: checkmol.pas, line 9653:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
      /*$IFDEF debug */
      //debugoutput ("matchpath too int - exiting");
      /*$ENDIF */
      return false;
    }
  /* next, check which chain of the needle matches which chain of the haystack  */
  for (i = 0; i < ndl_n_nb; i++)
    {
      emptyline = true;
      next_ndl_a = ndl_nb[i];
      for (j = 0; j < n_nb; j++)
	{
	  next_hst_a = hst_nb[j];
	  ndl_mp[ndl_mp_len] = next_ndl_a;
	  hst_mp[hst_mp_len] = next_hst_a;
	  if (is_matching (ndl_mp, hst_mp))
	    {			/* recursive function call */

	      if (max_match_recursion_depth != 0
		  && ++recursion_depth > max_match_recursion_depth)
		{
#ifndef MAKE_SHARED_LIBRARY
		  if (opt_verbose)
#endif
		    printf
		      ("Warning: max. number of match recursions (%i) reached, reverting to non-exhaustive match\n",
		       max_match_recursion_depth);
		  //n_rings = max_rings;
		  return true;
		}

	      mm[i][j] = true;
	      emptyline = false;
	    }
	}
      /* if a needle substituent does not match any of the haystack substituents, */
      /* stop any further work immediately */
      if (emptyline)
	return false;
    }
  /* finally, check the content of the matrix */
  res = matrix_OK (mm, ndl_n_nb, n_nb);
  /* optional: chirality check */
  if (!((res && (rs_search || ndl_atom[ndl_a - 1].stereo_care)) &&
	ndl_maybe_chiral (ndl_a)))
    return res;
  /* first, we have to clean up the match matrix in order to remove */
  /* "impossible" multiple matches (new in v0.3h) */
  for (i = 1; i <= 3; i++)
    {
      for (j = 1; j <= max_neighbors; j++)
	{			/* haystack dimension */
	  n_hits = 0;
	  l = 0;
	  for (k = 1; k <= max_neighbors; k++)
	    {			/* needle dimension */
	      if (mm[k - 1][j - 1])
		{
		  n_hits++;
		  l = k;
		}
	    }
	  if (n_hits == 1)
	    {			/* a unique match ==> kick out any other match at this pos. */
	      for (m = 1; m <= max_neighbors; m++)
		{
		  if (m != j)
		    mm[l - 1][m - 1] = false;
		}
	    }
	}
    }
  /* end of match matrix clean-up */
  if (prev_ndl_a > 0)
    {
      n_singlehits = 1;
      ncp[1] = prev_ndl_a;
      hcp[1] = prev_hst_a;
    }
  else
    n_singlehits = 0;
  ncp[0] = ndl_a;
  hcp[0] = hst_a;
  i = 0;
  l = 0;
  while (n_singlehits < 3 && i < 4)
    {
      i++;
      n_hits = 0;
      for (k = 1; k <= n_nb; k++)
	{
	  if (mm[i - 1][k - 1])
	    {
	      n_hits++;
	      l = k;
	    }
	}
      if (n_hits == 1)
	{
	  n_singlehits++;
	  ncp[n_singlehits] = ndl_nb[i - 1];
	  hcp[n_singlehits] = hst_nb[l - 1];
	}
    }
  if (n_singlehits != 3)
    return res;
  if (!chirality_OK (ncp, hcp))
    {
/* p2c: checkmol.pas, line 9749:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
      /*$IFDEF debug */
      //debugoutput ("chirality check failed");
      /*$ENDIF */
      res = false;
    }
  else
    {
/* p2c: checkmol.pas, line 9755:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
      /*$IFDEF debug */
      //debugoutput ("chirality check OK");
      /*$ENDIF */
    }
  return res;
/* p2c: checkmol.pas, line 9762:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
  /*$IFDEF debug */
  /* if res then tmpstr := ' MATCH' else tmpstr := ' NO MATCH';
     debugoutput('result for atoms '+inttostr(ndl_a)+'/'+inttostr(hst_a)+', bonds '+inttostr(ndl_b)+'/'+inttostr(hst_b)+':'+tmpstr); */
  /*$ENDIF */
}


static boolean
quick_match ()
{
  /* added in v0.2c */
  int i;
  boolean res = true;
  str3 ndl_atype;
  str2 ndl_el;			/* v0.3l */
  int ndl_chg;			/* v0.3l */
  int ndl_rad;			/* v0.3x */
  int ndl_iso;			/* v0.3x */


  if ((ez_search || rs_search) && ndl_n_heavyatoms > 3)
    /* v0.3f, v0.3m, v0.3o */
    return false;
  if (ndl_n_atoms < 1 || n_atoms < 1 || ndl_n_atoms > n_atoms ||
      ndl_n_bonds > n_bonds)
    {				/* just to be sure... */
/* p2c: checkmol.pas, line 9786:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
      /*$IFDEF debug */
/* p2c: checkmol.pas: Note: Eliminated unused assignment statement [338] */
      //debugoutput (" ==> quick_match failed");
      /*$ENDIF */
      return false;
    }

  if (ndl_n_heavyatoms > 1)
    {
      for (i = 0; i < ndl_n_atoms; i++)
	{
	  /*if atom^[i].atype <> ndl_atom^[i].atype then res := false;    (* changed in */
	  if (strcmp (atom[i].element, ndl_atom[i].element))	/* v0.2k */
	    return false;
	  //  if (atom[i].formal_charge != ndl_atom[i].formal_charge) /* v0.3o */
	  //res = false;


	  if (opt_chg)
	    {
	      if (ndl_atom[i].formal_charge != atom[i].formal_charge)
		return false;
	    }
/*  else
    {
      if (ndl_atom[i].formal_charge != 0 &&
	  atom[i].formal_charge != 0 &&
	  ndl_atom[i].formal_charge != atom[i].formal_charge)
	return false;
    } */

	  /* v0.3x: isotopes must be the same */
	  if (opt_iso)
	    {
	      if (ndl_atom[i].nucleon_number != atom[i].nucleon_number)
		return false;
	    }
/*  else
    {
      if (ndl_atom[i].nucleon_number != 0 &&
	  atom[i].nucleon_number != 0 &&
	  ndl_atom[i].nucleon_number !=
	  atom[i].nucleon_number)
	return false;
    }*/

	  /* v0.3x: radicals must be the same */
	  if (opt_rad)
	    {
	      if (ndl_atom[i].radical_type != atom[i].radical_type)
		return false;
	    }
/*  else
    {
      if (ndl_atom[i].radical_type != 0 &&
	  atom[i].radical_type != 0 &&
	  ndl_atom[i].radical_type != atom[i].radical_type)
	return false;
    }*/

	}
/* p2c: checkmol.pas, line 9798:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
      /*$IFDEF debug */
      //if (res)
      //debugoutput (" ==> quick_match: atoms OK");
      //else
      //  debugoutput (" ==> quick_match: atoms not OK");
      /*$ENDIF */
      if (ndl_n_bonds > 0)
	{

	  for (i = 0; i < ndl_n_bonds; i++)
	    {
	      if (ndl_bond[i].a1 != bond[i].a1 || ndl_bond[i].a2 != bond[i].a2
		  || ndl_bond[i].btype != bond[i].btype)
		return false;
	    }
	}
/* p2c: checkmol.pas, line 9810:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
      /*$IFDEF debug */
      //if (res)
      //  debugoutput (" ==> quick_match: bonds OK");
      //else
      //  debugoutput (" ==> quick_match: bonds not OK");
      /*$ENDIF */
      /* added in v0.2d: special case: needle contains only one heavy atom; refined in v0.3l, v0.3o */
    }
  else
    {

      /* first, find out the element and atom type of the only heavy atom       */
      for (i = 0; i < ndl_n_atoms; i++)
	{
	  if (ndl_atom[i].heavy)
	    {
	      //strcpy (ndl_atype, ndl_atom[i].atype);
	      strcpy (ndl_el, ndl_atom[i].element);	/* v0.3l */
	      ndl_chg = ndl_atom[i].formal_charge;	/* v0.3l */
	      ndl_iso = ndl_atom[i].nucleon_number;	/* 0.3x */
	      ndl_rad = ndl_atom[i].radical_type;	/* 0.3x */
	    }
	}

      for (i = 0; i < n_atoms; i++)
	{			/* v0.3l, v0.3o */
	  if (		//	!strcmp (atom[i].atype, ndl_atype) && 
	       !strcmp (atom[i].element, ndl_el))
	    {
	    

	      if (opt_chg || opt_strict)
		{
		  if (ndl_chg != atom[i].formal_charge)
		    return false;
		}


	      if (opt_iso || opt_strict)
		{
		  if (ndl_iso != atom[i].nucleon_number)
		    return false;
		}



	      if (opt_rad || opt_strict)
		{
		  if (ndl_rad != atom[i].radical_type)
		    return false;
		}
		  return true;
	    } else {
	        res=false;
	    }    
	}
    }
/* p2c: checkmol.pas, line 9828:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
  /*$IFDEF debug */
  //if (res)
  //  debugoutput (" ==> quick_match succeeded");
  //else
  //  debugoutput (" ==> quick_match failed (2)");
  /*$ENDIF */
  return res;
}

static void
perform_match ()
{
  int i = 0;
  int j;
  /*ndl_ref_atom : integer;  (* since v0.3j as a global variable */
  int ndl_n_nb, ndl_n_hc, n_nb, n_hc;
  boolean qm;			/* v0.3l */
  /* check for NoStruct (0 atoms);  v0.3l */
  if (n_atoms == 0 || ndl_n_atoms == 0)
    {
      matchresult = false;
/* p2c: checkmol.pas, line 9849:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
      /*$IFDEF debug */
      //debugoutput ("NoStruct encountered - aborted match routine");
      /*$ENDIF */
      return;
    }
  /* if we perform an exact match, needle and haystack must have */
  /* the same number of atoms, bonds, and rings */
  if (opt_exact && opt_iso)	/* 0.3x */
    {
      if (n_heavyatoms != ndl_n_heavyatoms
	  || n_heavybonds != ndl_n_heavybonds)
	{
	  matchresult = false;
/* p2c: checkmol.pas, line 9861:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
	  /*$IFDEF debug */
	  //debugoutput ("different number of heavy atoms and/or bonds");
	  /*$ENDIF */
	  //return;
	}
    }

  /* have a quick look if needle and haystack are identical molfiles */
  qm = quick_match ();		/* v0.3l */
  if (qm)
    {
      matchresult = true;
      clear_ndl_atom_tags ();	/* v0.3o */
      return;
    }
  /* if we have only one heavy atom and quick_match fails, return "false";  v0.3l */
  if (ndl_n_heavyatoms == 1)
    {
      matchresult = false;
      return;
    }
/* p2c: checkmol.pas, line 9881:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
  /*$IFDEF debug */
  /* debugoutput('needle reference atom: '+inttostr(ndl_ref_atom)+' ('+ndl_atom^[ndl_ref_atom].atype+')'); */
  /*$ENDIF */
  ndl_n_nb = ndl_atom[ndl_ref_atom - 1].neighbor_count;
  ndl_n_hc = ndl_hetatom_count (ndl_ref_atom);
/* p2c: checkmol.pas, line 9886:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
  /*$IFDEF debug */
  /* debugoutput('neighbor atoms: '+inttostr(ndl_n_nb)+'  heteroatom neighbors: '+inttostr(ndl_n_hc)); */
  /*$ENDIF */
  matchresult = false;
  for (j = 0; j < max_matchpath_length; j++)
    {
      ndl_matchpath[j] = 0;
      hst_matchpath[j] = 0;
    }
  ndl_matchpath[0] = ndl_ref_atom;
  while (i < n_atoms && matchresult == false)
    {
      i++;
      n_nb = atom[i - 1].neighbor_count;
      n_hc = hetatom_count (i);
      if (n_nb >= ndl_n_nb && n_hc >= ndl_n_hc)
	{
/* p2c: checkmol.pas, line 9904:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
	  /*$IFDEF debug */
	  /* debugoutput('trying atom '+inttostr(i)+'; neighbor atoms: '+inttostr(n_nb)+' heteroatom neighbors: '+inttostr(n_hc)); */
	  /*$ENDIF */

	  recursion_depth = 0;
	  hst_matchpath[0] = i;
	  matchresult = is_matching (ndl_matchpath, hst_matchpath);
/* p2c: checkmol.pas, line 9909:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
	  /*$IFDEF debug */
	  /* if matchresult then debugoutput('matching atom in haystack: '+inttostr(i)+' ('+atom^[i].atype+')'); */
	  /*$ENDIF */
	  if (matchresult)	/* v0.3o; mark this fragment as matched */
	    atom[i - 1].tag = true;
	}
    }
}


static void
clear_rings ()
{
  int i, FORLIM;
  n_rings = 0;
  memset (ring, 0, sizeof (ringlist));
  for (i = 0; i < max_rings; i++)
    {				/* new in v0.3 */
      ringprop[i].size = 0;
      ringprop[i].arom = false;
      ringprop[i].envelope = false;
    }
  if (n_atoms > 0)
    {
      FORLIM = n_atoms;
      for (i = 0; i < FORLIM; i++)
	atom[i].ring_count = 0;
    }
  if (n_bonds > 0)
    {
      FORLIM = n_bonds;
      for (i = 0; i < FORLIM; i++)
	bond[i].ring_count = 0;
    }
}


static int
ring_lastpos (s)
     int *s;
{
  int i, rc;
  int rlp = 0;
  int FORLIM;
  if (n_rings <= 0)
    return rlp;
  FORLIM = n_rings;
  for (i = 1; i <= FORLIM; i++)
    {
      rc = ringcompare (s, ring[i - 1]);
      if (rc_identical (rc))
	rlp = i;
    }
  return rlp;
}


static void
remove_redundant_rings ()
{
  int i, j, k, rlp;
  ringpath_type tmp_path;
  int FORLIM, FORLIM1;
  if (n_rings < 2)
    return;
  FORLIM = n_rings;
  for (i = 1; i < FORLIM; i++)
    {
      memcpy (tmp_path, ring[i - 1], sizeof (ringpath_type));
      rlp = ring_lastpos (tmp_path);
      while (rlp > i)
	{
	  FORLIM1 = n_rings;
/* p2c: checkmol.pas, line 9970:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
	  /*$IFDEF debug */
	  /* debugoutput('removing redundant ring: '+inttostr(rlp)+' (identical to ring '+inttostr(i)+')'); */
	  /*$ENDIF */
	  for (j = rlp; j < FORLIM1; j++)
	    {
	      memcpy (ring[j - 1], ring[j], sizeof (ringpath_type));
	      ringprop[j - 1].size = ringprop[j].size;	/* new in v0.3 */
	      ringprop[j - 1].arom = ringprop[j].arom;
	      ringprop[j - 1].envelope = ringprop[j].envelope;
	    }
	  for (k = 0; k < max_ringsize; k++)
	    ring[n_rings - 1][k] = 0;
	  n_rings--;
	  rlp = ring_lastpos (tmp_path);
	}
    }
}


static int
count_aromatic_rings ()
{
  int i;
  int n = 0;
  int FORLIM;
  if (n_rings <= 0)
    return n;
  FORLIM = n_rings;
  for (i = 0; i < FORLIM; i++)
    {
      if (ringprop[i].arom)
	n++;
    }
  return n;
}


static void
chk_envelopes ()
{
  /* new in v0.3d */
  /* checks if a ring completely contains one or more other rings */
  int a, i, j, k, l, pl, pli;
  boolean found_atom, found_all_atoms, found_ring;
  int FORLIM;
  if (n_rings < 2)
    return;
  FORLIM = n_rings;
  for (i = 1; i < FORLIM; i++)
    {
      found_ring = false;
      j = 0;
      pli = ringprop[i].size;	/* path_length(ring^[i]); */
      while (j < i && found_ring == false)
	{
	  j++;
	  found_all_atoms = true;
	  pl = ringprop[j - 1].size;	/* path_length(ring^[j]); */
	  for (k = 0; k < pl; k++)
	    {
	      found_atom = false;
	      a = ring[j - 1][k];
	      for (l = 0; l < pli; l++)
		{
		  if (ring[i][l] == a)
		    found_atom = true;
		}
	      if (found_atom == false)
		found_all_atoms = false;
	    }
	  if (found_all_atoms)
	    found_ring = true;
	}
      if (found_ring)
	ringprop[i].envelope = true;
    }
}


static void
update_ringcount ()
{
  int i, j, a1, a2, b, pl, FORLIM;
  if (n_rings <= 0)
    return;
  chk_envelopes ();
  FORLIM = n_rings;
  for (i = 0; i < FORLIM; i++)
    {
      if (ringprop[i].envelope == false)
	{
	  pl = ringprop[i].size;	/* path_length(ring^[i]);  (* v0.3d */
	  a2 = ring[i][pl - 1];
	  for (j = 0; j < pl; j++)
	    {
	      a1 = ring[i][j];
	      atom[a1 - 1].ring_count++;
	      b = get_bond (a1, a2);
	      bond[b - 1].ring_count++;
	      a2 = a1;
	    }
	}
    }
}


static boolean
normalize_ionic_bonds ()
{
  /* v0.3k */
  /* changed from a procedure into a function in v0.3m */
  int i, a1, a2, fc1, fc2;
  char bt;
  boolean res = false;		/* v0.3m */
  int FORLIM;
  /* v0.3m */
  if (n_bonds == 0)
    return false;
  FORLIM = n_bonds;
  for (i = 0; i < FORLIM; i++)
    {
      a1 = bond[i].a1;
      a2 = bond[i].a2;
      bt = bond[i].btype;
      fc1 = atom[a1 - 1].formal_charge;
      fc2 = atom[a2 - 1].formal_charge;
      if (fc1 * fc2 == -1 && (bt == 'S' || bt == 'D'))
	{
	  atom[a1 - 1].formal_charge = 0;
	  atom[a2 - 1].formal_charge = 0;
	  if (!strcmp (atom[a1 - 1].atype, "N3+"))	/* v0.3m */
	    strcpy (atom[a1 - 1].atype, "N3 ");
	  if (!strcmp (atom[a2 - 1].atype, "N3+"))	/* v0.3m */
	    strcpy (atom[a2 - 1].atype, "N3 ");
	  if (bt == 'D')
	    bond[i].btype = 'T';
	  if (bt == 'S')
	    bond[i].btype = 'D';
	  res = true;		/* v0.3m */
	}
    }
  return res;			/* v0.3m (return true if any change was made */
}

static void
chk_wildcard_rings ()		// new in v0.3p
// checks if there are any wildcard atom types or bond types
// in a ring of the needle; if yes ==> set the q_arom flag in the
// atom and bond record of all ring members in order to perform the 
// match a bit more generously
{

  int i, j, rs;
  int a1, a2, b;
  boolean wcr;
  str3 at;
  char bt;

  if (ndl_querymol == false)
    return;
  if (ndl_n_rings == 0)
    return;
  // now look for any not-yet-aromatic rings which contain a wildcard
  for (i = 0; i < ndl_n_rings; i++)
    {
      wcr = false;
      if (ndl_ringprop[i].arom == false)
	{
	  rs = ndl_ringprop[i].size;
	  a2 = ndl_ring[i][rs];
	  for (j = 0; j < rs; j++)
	    {
	      a1 = ndl_ring[i][j];
	      b = get_ndl_bond (a1, a2);
	      strcpy (at, ndl_atom[a1].atype);
	      bt = ndl_bond[b].btype;
	      if (!strcmp (at, "A  ") || !strcmp (at, "Q  "))
		wcr = true;
	      if (bt == 'l' || bt == 's' || bt == 'd' || bt == 'a')
		wcr = true;
	      a2 = a1;
	    }
	  if (wcr)
	    {			// if yes, flag all atoms and bonds in this ring as "potentially" aromatic
	      // {$IFDEF debug}
	      // debugoutput('wildcard ring found');
	      // {$ENDIF}
	      a2 = ndl_ring[i][rs];
	      for (j = 0; j < rs; j++)
		{
		  a1 = ndl_ring[i][j];
		  b = get_ndl_bond (a1, a2);
		  strcpy (at, ndl_atom[a1].atype);
		  bt = ndl_bond[b].btype;
		  ndl_atom[a1].q_arom = true;
		  ndl_bond[b].q_arom = true;
		  a2 = a1;
		}
	    }
	}
    }
  // and now undo this flagging for all rings which contain no wildcard
  for (i = 0; i < ndl_n_rings; i++)
    {
      wcr = false;
      rs = ndl_ringprop[i].size;
      a2 = ndl_ring[i][rs];
      for (j = 0; j < rs; j++)
	{
	  a1 = ndl_ring[i][j];
	  b = get_ndl_bond (a1, a2);
	  strcpy (at, ndl_atom[a1].atype);
	  bt = ndl_bond[b].btype;
	  if (!strcmp (at, "A  ") || !strcmp (at, "Q  "))
	    wcr = true;
	  if (bt == 'l' || bt == 's' || bt == 'd' || bt == 'a')
	    wcr = true;
	  a2 = a1;
	}
      if (!wcr)
	{			// if yes, unflag all atoms and bonds in this ring
	  a2 = ndl_ring[i][rs];
	  for (j = 0; j < rs; j++)
	    {
	      a1 = ndl_ring[i][j];
	      b = get_ndl_bond (a1, a2);
	      strcpy (at, ndl_atom[a1].atype);
	      bt = ndl_bond[b].btype;
	      ndl_atom[a1].q_arom = false;
	      ndl_bond[b].q_arom = false;
	      a2 = a1;
	    }
	}
    }
  // some further refinement would be necessary here in order to unflag everything
  // which contains a wildcard but which definitely cannot be aromatic
}

#ifndef MAKE_SHARED_LIBRARY

int
main (int argc, char *argv[])
{				/* main routine */
  char STR1[256], STR6[256];
  int FORLIM;
  /* progmode = pmMatchMol */
  rfile = NULL;
  strcpy (progname, argv[0]);
  strncpy (STR1, progname, 253);
  if (strstr (STR1, "matchmol") != NULL)
    progmode = pmMatchMol;
  else
    {
      strncpy (STR6, progname, 253);
      if (strstr (STR6, "checkmol") == NULL)
	{
	  printf ("THOU SHALLST NOT RENAME ME!\n");
	  exit (9);
	}
      progmode = pmCheckMol;
    }
  if (argc == 1)
    {
      show_usage ();
      exit (1);
    }
  init_globals ();
  init_molstat (&molstat);
  parse_args (argc, argv);
  if (ringsearch_mode == rs_sar)
    max_vringsize = max_ringsize;
  else
    max_vringsize = ssr_vringsize;
  /* v0.3n (was: 10) */
  /*if opt_verbose then writeln(progname+' v',version,'  N. Haider 2003-2007'); */
  if (progmode == pmMatchMol)
    {
      left_trim (ndl_molfilename);
      left_trim (molfilename);
      if ((*molfilename == '\0' || *ndl_molfilename == '\0'
	   || argc < 3) && !opt_stdin)
	{

	  show_usage ();
	  exit (2);		/* new in v0.2k */
	}
      if (!(file_exists (ndl_molfilename)) && !opt_stdin)
	{			/*not  fileexists(ndl_molfilename) REPLACE!!! 
				   printf("2");
				   /* p2c: checkmol.pas, line 10128:
				   * Warning: Expected an expression, found a ')' [227] */
	  if (strlen (ndl_molfilename) > 1 && ndl_molfilename[0] == '-')
	    show_usage ();
	  else
	    printf ("file %s not found!\n", ndl_molfilename);
	  /* new in v0.2k */
	  exit (2);
	}
    }

  if (!(file_exists (molfilename)) && !opt_stdin)
    {				/*not  fileexists(ndl_molfilename) REPLACE!!! */
/* p2c: checkmol.pas, line 10128:
 * Warning: Expected an expression, found a ')' [227] */

      if (strlen (molfilename) > 1 && molfilename[0] == '-')
	show_usage ();
      else
	printf ("file %s not found!\n", molfilename);
      /* new in v0.2k */
      exit (2);
    }

  /* read the first molecule and process it; if we are in "matchmol" mode, */
  /* this is the "needle" */
  if (progmode == pmMatchMol)
    readinputfile (ndl_molfilename);
  else
    readinputfile (molfilename);
  li = 1;			/* initialize line pointer for input buffer */
  get_filetype (filetype, ndl_molfilename);
  if (!strcmp (filetype, "unknown"))
    {
      printf ("unknown query file format!\n");
      if (!opt_verbose)
	exit (3);
      printf ("===========================================\n");
      FORLIM = molbufindex;
      for (i = 1; i <= FORLIM; i++)
	puts (molbuf[i - 1]);
      exit (3);
    }
  mol_OK = true;		/* added in v0.2i */
  if (!strcmp (filetype, "alchemy"))
    read_molfile (ndl_molfilename);
  if (!strcmp (filetype, "sybyl"))
    read_mol2file (ndl_molfilename);
  if (!strcmp (filetype, "mdl"))
    read_MDLmolfile (ndl_molfilename);
  count_neighbors ();
  if (!mol_OK || n_atoms < 1)
    {				/* v0.3g; check if this is a valid query structure */
      printf ("invalid molecule\n");
      exit (3);
    }
  if (!found_arominfo || progmode == pmCheckMol)
    {				/* added in v0.2b/0.2c */
/* p2c: checkmol.pas, line 10172:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
      /*$IFDEF debug 
         if (!found_arominfo)
         debugoutput
         ("no aromaticity information found - checking myself...");
         else
         debugoutput ("performing full aromaticity check");
         /* new in v0.3d */
      /*$ENDIF */
      chk_ringbonds ();
      if (ringsearch_mode == rs_ssr)
	remove_redundant_rings ();
      if (n_rings >= max_rings)
	{
	  if (opt_verbose)
	    printf
	      ("Warning: max. number of rings (%i) reached, reverting to SSR search\n",
	       max_rings);
	  ringsearch_mode = rs_ssr;
	  auto_ssr = true;	/* v0.3n */
	  clear_rings ();
	  max_vringsize = ssr_vringsize;	/* v0.3n (was: 10) */
	  chk_ringbonds ();
	  remove_redundant_rings ();
	}
      update_ringcount ();
      /* new in v0.3k: if output is a molfile, leave the original */
      /* representation of N-oxides, S-oxides, nitro groups, etc. */
      /* unchanged (ionic or non-ionic), in any other case make covalent bonds */
      if (!opt_xmdlout)		/* v0.3k */
	normalize_ionic_bonds ();
      update_atypes ();
      update_Htotal ();		/* added in v0.3 */
      chk_arom ();
      if (ringsearch_mode == rs_ssr)
	{			/* new in v0.3 */
	  do
	    {
	      prev_n_ar = count_aromatic_rings ();
	      chk_arom ();
	      n_ar = count_aromatic_rings ();
	    }
	  while (prev_n_ar - n_ar != 0);
	}
    }
  else
    {				/* v0.3k */
/* p2c: checkmol.pas, line 10206:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
      /*$IFDEF debug 
         debugoutput ("found aromaticity information in input file");
         /*$ENDIF */
      if (!opt_xmdlout)
	normalize_ionic_bonds ();
      update_atypes ();		/* added in v0.2f */
      update_Htotal ();		/* end v0.2b snippet */
    }
  if (progmode == pmCheckMol)
    {
      if (opt_verbose)
	write_mol ();
      get_molstat ();
      if (opt_molstat)
	{
	  if (opt_molstat_X)
	    write_molstat_X ();
	  else
	    write_molstat ();
	}
      else
	{
	  if (found_querymol)
	    {
	      printf ("input structure contains query atom or query bond!\n");
	      exit (1);
	    }
	  chk_functionalgroups ();
	  if (opt_none)
	    opt_text = true;
	  if (opt_text)
	    write_fg_text ();
	  if (opt_text_de)
	    write_fg_text_de ();
	  if (opt_code)
	    write_fg_code ();
	  if (opt_bin)
	    write_fg_binary ();
	  if (opt_bitstring)
	    write_fg_bitstring ();
	  if (opt_xmdlout)
	    write_MDLmolfile ();
	}
      /*if opt_verbose   then write_mol; */
      zap_molecule ();
    }
  else
    {
      /* now transfer all data to the "needle" set of variables, except for "fingerprint" mode */
      if (!opt_fp)
	{			/* v0.3m */
	  copy_mol_to_needle ();
	  //chk_wildcard_rings (); /* 0.3p */
	  set_ndl_atom_tags ();	/* v0.3o */
	  if (opt_verbose)
	    write_needle_mol ();
	  if (rs_strict)	/* v0.3j */
	    ndl_ref_atom = find_ndl_ref_atom_cv ();
	  else
	    ndl_ref_atom = find_ndl_ref_atom ();
	}
      else
	{
	  copy_mol_to_tmp ();	/* v0.3m */
	  if (opt_verbose)
	    printf ("1st molecule stored in buffer: %s\n", tmp_molname);
	}
      /* next, read the "haystack" file and process it */
      li = 1;
      mol_count = 0;
      fpdecimal = 0;		/* v0.3m */
      fpindex = 0;		/* v0.3m */
      do
	{
	  /* new in v0.3i: reset ringsearch_mode to its initial value */
	  /* for each new molecule */
	  ringsearch_mode = opt_rs;
	  if (ringsearch_mode == rs_sar)
	    max_vringsize = max_ringsize;
	  else
	    max_vringsize = ssr_vringsize;
	  /* v0.3n (was: 10) */
	  readinputfile (molfilename);
	  li = 1;
	  get_filetype (filetype, molfilename);
	  if (strcmp (filetype, "unknown"))
	    {
	      found_arominfo = false;	/* added in v0.2b */
	      mol_OK = true;	/* added in v0.2i */
	      if (!strcmp (filetype, "alchemy"))
		read_molfile (molfilename);
	      if (!strcmp (filetype, "sybyl"))
		read_mol2file (molfilename);
	      if (!strcmp (filetype, "mdl"))
		read_MDLmolfile (molfilename);
	      mol_count++;
	      fpindex++;
	      count_neighbors ();
	      /*if (not mol_OK) or (n_atoms < 1) then writeln(mol_count,':no valid structure found') else */
	      if (!mol_OK || n_atoms < 1
		  && !(opt_fp && fpformat == fpf_decimal))
		printf ("%i:F\n", mol_count);
	      else
		{
		  if (opt_exact
		      && (n_Ctot != ndl_n_Ctot || n_Otot != ndl_n_Otot
			  || n_Ntot != ndl_n_Ntot))
		    {		/* new in v0.3g */
		      if (!opt_molout && !(opt_fp && fpformat == fpf_decimal))
			printf ("%i:F\n", mol_count);
		    }
		  else
		    {
		      if (!found_arominfo || opt_strict && tmfmismatch)
			{	/* added in v0.3m */
/* p2c: checkmol.pas, line 10294:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
			  /*$IFDEF debug 
			     debugoutput
			     ("no aromaticity information found (or tweak mismatch) - checking myself...");
			     /*$ENDIF */
			  chk_ringbonds ();
			  if (ringsearch_mode == rs_ssr)
			    remove_redundant_rings ();
			  if (n_rings == max_rings)
			    {
			      if (opt_verbose)
				printf
				  ("Warning: max. number of rings (%i) reached, reverting to SSR search\n",
				   max_rings);
			      ringsearch_mode = rs_ssr;
			      clear_rings ();
			      max_vringsize = ssr_vringsize;	/* v0.3n (was: 10) */
			      chk_ringbonds ();
			      remove_redundant_rings ();
			    }
			  update_ringcount ();
			  update_atypes ();
			  update_Htotal ();	/* added in v0.3 */
			  chk_arom ();
			  if (ringsearch_mode == rs_ssr)
			    {	/* new in v0.3 */
			      do
				{
				  prev_n_ar = count_aromatic_rings ();
				  chk_arom ();
				  n_ar = count_aromatic_rings ();
				}
			      while (prev_n_ar - n_ar != 0);
			    }
			}
		      else
			{	/* added in v0.2f */
/* p2c: checkmol.pas, line 10322:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
			  /*$IFDEF debug 
			     debugoutput
			     ("found aromaticity information in input file");
			     /*$ENDIF */
			  if (opt_strict)
			    update_atypes ();
			  update_Htotal ();
			}
		      init_molstat (&ndl_molstat);
		      if (normalize_ionic_bonds ())
			/* new in v0.3k, modified in v0.3m */
			update_atypes ();
		      if (opt_verbose && !opt_fp)
			write_mol ();
		      /* if in "fingerprint mode", exchange needle and haystack */
		      if (opt_fp)
			{	/* v0.3m */
			  zap_needle ();
			  copy_mol_to_needle ();
			  //chk_wildcard_rings (); /* 0.3p */
			  zap_molecule ();
			  copy_tmp_to_mol ();
			  if (opt_verbose)
			    write_needle_mol ();
			  if (rs_strict)	/* v0.3j */
			    ndl_ref_atom = find_ndl_ref_atom_cv ();
			  else
			    ndl_ref_atom = find_ndl_ref_atom ();
			  if (opt_verbose)
			    write_mol ();
			}	/* v0.3m */
		      /* now that we have both molecules, perform the comparison */
		      /* v0.3o: takes care of disconnected fragment... */
		      clear_atom_tags ();
		      set_ndl_atom_tags ();
		      matchsummary = true;
		      perform_match ();
		      matchsummary = matchresult;
		      if (count_tagged_ndl_heavyatoms () > 0
			  && matchsummary == true)
			{
			  do
			    {
			      if (rs_strict)
				ndl_ref_atom = find_ndl_ref_atom_cv ();
			      else
				ndl_ref_atom = find_ndl_ref_atom ();
			      perform_match ();
			      if (matchresult == false)
				matchsummary = false;
			    }
			  while (count_tagged_ndl_heavyatoms () != 0 &&
				 matchsummary != false);
			}
		      /* end of disconnected-fragment matching (v0.3o) */
		      if (matchsummary == true)
			{	/* v0.3o */
			  if (opt_molout)
			    {
			      FORLIM = molbufindex;
			      for (i = 1; i <= FORLIM; i++)
				puts (molbuf[i - 1]);
			    }
			  else
			    {
			      if (!opt_fp)	/* inttostr(mol_count) REPLACE!!!, */
				printf ("%i:T\n", mol_count);
			      else
				{
				  if (ndl_n_heavyatoms == n_heavyatoms &&
				      ndl_n_heavybonds == n_heavybonds)
				    fp_exacthit = true;
				  else
				    fp_exacthit = false;
				  if (fp_exacthit)
				    fp_exactblock = true;
				  if (fpformat == fpf_boolean)
				    {
				      if (fp_exacthit)	/* inttostr(mol_count), REPACE!!! */
					printf ("%i:TX\n", mol_count);
				      else
					printf ("%i:T\n", mol_count);
				    }
				  /* inttostr(mol_count), REPLACE!!! */
				  if (fpformat == fpf_decimal)
				    {
				      fpincrement = 1;
				      FORLIM = fpindex;
				      for (i = 1; i <= FORLIM; i++)
					fpincrement <<= 1;
				      fpdecimal += fpincrement;
				    }
				}
			    }
			}
		      else
			{
			  if (!
			      (opt_molout || opt_fp
			       && fpformat == fpf_decimal))
			    /* inttostr(mol_count), REPLACE!!! */
			    printf ("%i:F\n", mol_count);
			}
		      if (opt_fp && fpformat == fpf_decimal
			  && fpindex == fp_blocksize)
			{
			  if (fp_exactblock)
			    fpdecimal++;
			  printf ("%i\n", fpdecimal);
			  fpindex = 0;
			  fpdecimal = 0;
			  fp_exactblock = false;
			}
		      zap_molecule ();
		      molbufindex = 0;
		    }
		}
	    }
	  else
	    {
	      /* v0.3l */
	      /* mol_OK */
	      printf ("%i:unknown file format\n", mol_count);
	    }
	}
      while (mol_in_queue != false);
      /* if filetype <> 'unknown' */
      if (opt_fp && fpformat == fpf_decimal && fpindex > 0)
	{
	  if (fp_exactblock)
	    fpdecimal++;
	  printf ("%i\n", fpdecimal);
	}
      zap_needle ();
      if (rfile_is_open)
	{			/* new in v0.2g */
	  if (rfile != NULL)
	    fclose (rfile);
	  rfile = NULL;
	}
    }
  if (rfile != NULL)
    fclose (rfile);
  exit (0);
}

#else

static void
init_globals_dll (void)
{

//printf("init_globals_dll\n");

  int i;
  opt_verbose = false;
  opt_debug = false;
  opt_stdin = false;
  opt_text = false;
  opt_code = false;
  opt_bin = false;
  opt_bitstring = false;
  opt_molout = false;
  opt_molstat = false;
  opt_molstat_X = false;
  opt_xmdlout = false;
  opt_fp = false;		/* new in v0.3m */
  /*cm_mdlmolfile   := false; */
  found_arominfo = false;
  found_querymol = false;
  ndl_querymol = false;
  opt_rs = rs_sar;		/* v0.3i */
  ringsearch_mode = opt_rs;
  rfile_is_open = false;	/* new in v0.2g */
  ez_flag = false;		/* new in v0.3f */
  chir_flag = false;		/* new in v0.3f */
  n_Ctot = 0;
  n_Otot = 0;
  n_Ntot = 0;			/* new in v0.3g */
  //for (i = 0; i < max_fg; i++)
  //  fg[i] = false;
  memset (fg, 0, sizeof (fglist));

  if (!yet_initialized)
    {
      molbuf = (void *) safe_malloc (sizeof (molbuftype));
      opt_exact = false;
      opt_strict = false;	/* new in v0.2f */
      opt_metalrings = false;	/* new in v0.3 */
      opt_geom = false;		/* new in v0.3d */
      opt_chiral = false;	/* new in v0.3f */
      opt_iso = false;		/* new in v0.3x */
      opt_chg = false;		/* new in v0.3x */
      opt_rad = false;		/* new in v0.3x */
      ez_search = false;	/* new in v0.3d */
      rs_search = false;	/* new in v0.3f */
      rs_strict = false;	/* new in v0.3j */
      ndl_n_Ctot = 0;
      ndl_n_Otot = 0;
      ndl_n_Ntot = 0;		/* new in v0.3g */
      yet_initialized = true;
    }

  ether_generic = false;	/* v0.3j */
  amine_generic = false;	/* v0.3j */
  hydroxy_generic = false;	/* v0.3j */
  fpformat = fpf_decimal;	/* v0.3m */
  fpindex = 0;			/* v0.3m */
  fp_exacthit = false;		/* v0.3m */
  fp_exactblock = false;	/* v0.3m */
  tmfcode = 0;			/* v0.3m */
  tmfmismatch = false;		/* v0.3m */
  auto_ssr = false;
  recursion_depth = 0;
}

static void
mm_init_mol (void)
{

//printf("mm_init_mol\n");
  init_globals_dll ();
  init_molstat (&molstat);
  if (opt_rs_dll == RPA_DEFAULT)
    {
      ringsearch_mode = opt_rs;
      //printf("DEFAULT: %i\n",ringsearch_mode);
    }
  else
    {
      ringsearch_mode = opt_rs_dll;
    }
  //printf("RPA: %i\n",ringsearch_mode);

  if (ringsearch_mode == rs_sar)
    {
      max_vringsize = max_ringsize;
    }
  else
    {
      max_vringsize = ssr_vringsize;
    }
  zap_molecule ();
  molbufindex = 0;
  mol_count = 0;
//printf("mm_init_mol\n");
}

static void
mm_elab_mol (boolean checkmol_mode, boolean normalize_ionic_bnds)
{
//printf("mm_elab_mol\n");

  li = 1;			// initialize line pointer for input buffer
  get_filetype (filetype, ndl_molfilename);
  if (strcmp (filetype, "unknown") == 0)
    {
      //messagebox (0,'Error in mm_ElabMol: Unknown file format','MATCHMOLDLL ERROR',0);
      exit (3);
    }

  if (checkmol_mode == true)
    progmode = pmCheckMol;
  else
    progmode = pmMatchMol;
  if (strcmp (filetype, "alchemy") == 0)
    read_molfile (ndl_molfilename);
  if (strcmp (filetype, "sybyl") == 0)
    read_mol2file (ndl_molfilename);
  if (strcmp (filetype, "mdl") == 0)
    read_MDLmolfile (ndl_molfilename);
  if (checkmol_mode)
    {
      if (found_querymol)
	{
	  printf
	    ("Warning: Input structure contains query atom or query bond.\n");
	}
    }

  count_neighbors ();
  if (!found_arominfo || checkmol_mode || opt_strict)
    {
      //printf("No arom found or checkmol mode\n");
      chk_ringbonds ();
      if (ringsearch_mode == rs_ssr)
	remove_redundant_rings ();
      if (n_rings >= max_rings)
	{

	  printf
	    ("Warning: max. number of rings (%i) reached, reverting to SSR search\n",
	     max_rings);
	  ringsearch_mode = rs_ssr;
	  auto_ssr = true;
	  clear_rings ();
	  max_vringsize = ssr_vringsize;
	  chk_ringbonds ();
	  remove_redundant_rings ();
	}

      update_ringcount ();
      if (normalize_ionic_bnds)	/* v0.3k */
	normalize_ionic_bonds ();
      update_atypes ();
      update_Htotal ();
      chk_arom ();
      if (ringsearch_mode == rs_ssr)
	{			/* new in v0.3 */
	  do
	    {
	      prev_n_ar = count_aromatic_rings ();
	      chk_arom ();
	      n_ar = count_aromatic_rings ();
	    }
	  while (prev_n_ar - n_ar != 0);
	}
    }
  else
    {
      if (normalize_ionic_bnds)	/* v0.3k  */
	normalize_ionic_bonds ();
      //if (opt_strict)
      update_atypes ();
      update_Htotal ();
    }



//printf("mm_elab_mol\n");
}

DLLEXPORT void
mm_set_current_mol_as_query (void)
{
//printf("mm_set_current_mol_as_query\n");
  zap_needle ();
//mm_ElabMol;
  copy_mol_to_needle ();
  //chk_wildcard_rings (); /* 0.3p */
  set_ndl_atom_tags ();		/* v0.3o */
  if (opt_geom)			/* v0.3d */
    ez_search = true;
  else if (!ez_flag && ez_search)
    ez_search = false;		//chir_flag initialized in read_MDLmolfile, otherwise we lose that info
  if (opt_chiral)
    {				/* v0.3f */
      rs_search = true;

      //printf("%i\n",rs_search);
    }
  else if (!chir_flag && rs_search)
    {
      rs_search = false;	//chir_flag initialized in read_MDLmolfile, otherwise we lose that info
      //printf("%i\n",rs_search);
    }
  if (opt_chiral && opt_strict && opt_exact)	/* new in v0.3j */
    rs_strict = true;
  else
    rs_strict = false;
  /* if (rs_strict)              /* v0.3j 
     ndl_ref_atom = find_ndl_ref_atom_cv ();
     //ndl_ref_atom = find_ndl_ref_atom ();
     else
     ndl_ref_atom = find_ndl_ref_atom (); */



  molbufindex = 0;
  mol_count = 0;
//printf("mm_set_current_mol_as_query\n");
}

DLLEXPORT int
mm_get_rings (void)
{
  return n_rings;
}

DLLEXPORT void
xm_version (char *buffer)
{
  buffer[0] = '\0';
  strncpy (buffer, version, 255);
}

DLLEXPORT int
mm_get_atom_ring (int atom_number)
{

  int i, j, a1, pl;
  int ret = 0;
  a1 = atom[atom_number].ring_count;
  if (n_rings > 0)
    {
      for (i = 1; i < n_rings; i++)
	{

	  pl = path_length (ring[i]);
//          a2 := ring^[i,pl];
	  for (j = 1; j < pl; j++)
	    {

	      a1 = ring[i][j];
	      if (atom_number == a1)
		ret = i;
//
//              inc(atom^[a1].ring_count);
//              b := get_bond(a1,a2);
//              inc(bond^[b].ring_count);
//              a2 := a1;
	    }
	}
    }
  return ret;
}

static void
mm_read_input_line (char *st)
{
//printf("mm_read_input_line_in\n");
//var
//yyy:pchar;

  mol_in_queue = false;
  if (molbufindex < (max_atoms + max_bonds + slack))
    {

//yyy:=Pchar(IntToStr(molbufindex));
//messagebox (0,yyy,'',0);
//printf("%i\n",molbufindex);
//printf("B:%s\n",st);
      strcpy (molbuf[molbufindex++], st);
//printf("%x %x\n",&molbuf,molbuf);
//printf("%s\n",molbuf[molbufindex-1]);
      //  molbufindex++;
    }
  else
    {
      //messagebox(0,'Error in mm_Readinputline; memory problem','ERROR',0);
      printf ("Not enough memory for molfile! %i\n", molbufindex);
      exit (1);
    }
//printf("mm_read_input_line_out\n");
}


static void
mm_set_mol_dll (const char *st, boolean checkmol_mode,
		boolean normalize_ionic_bnds)
{
//printf("mm_set_mol\n");
//printf("%s\n",st);
  char bb;
  char aa;
  int i;
  int k;
  int J;
  int spt = 0;
  char tt[256];
  char bb10 = '\n';
  char bb13 = '\r';
  char bb0 = '\0';
  int lenst;
//char d[256];
  lenst = strlen (st);
//tt=(char*)safe_malloc(256*sizeof(char));
  tt[0] = '\0';
//messagebox(0,st,'',0);
  mm_init_mol ();
  for (i = spt; i < lenst; i++)
    {
      bb = st[i];
      if ((bb == bb10) || (i == lenst))
	{
	  J = 0;
	  // d:='';
	  for (k = spt; k < i; k++)
	    {

	      aa = st[k];
	      if ((aa != bb10) && (aa != bb13))
		{
		  //d:=d+aa;
		  tt[J] = aa;
		  J++;
		}
	    }
	  tt[J] = bb0;
	  spt = i;
	  //messagebox (0,tt,tt,0);
//printf("A:%s\n",tt);
	  mm_read_input_line (tt);
	}
    }
//free(tt);
  mm_elab_mol (checkmol_mode, normalize_ionic_bnds);
//printf("mm_set_mol\n");
}

DLLEXPORT void
cm_set_mol (const char *st, int normalize_ionic_bnds)
{
  mm_set_mol_dll (st, true, (normalize_ionic_bnds != FEATURE_OFF));
}

DLLEXPORT void
mm_set_mol (const char *st)
{
  mm_set_mol_dll (st, false, true);
}

DLLEXPORT void
xm_set_strict_typing (int strict_typing)
{
  if (!yet_initialized)
    init_globals_dll ();
  opt_strict = (strict_typing != FEATURE_OFF);
  //opt_strict=false; //This never worked right and is harmful
}

DLLEXPORT void
mm_set_r_s_check (int r_s_check)
{
  if (!yet_initialized)
    init_globals_dll ();
  opt_chiral = (r_s_check != FEATURE_OFF);
}

DLLEXPORT void
mm_set_e_z_check (int e_z_check)
{
  if (!yet_initialized)
    init_globals_dll ();
  opt_geom = (e_z_check != FEATURE_OFF);
}

DLLEXPORT void
mm_set_chg_check (int chg_check)
{
  if (!yet_initialized)
    init_globals_dll ();
  opt_chg = (chg_check != FEATURE_OFF);
}

DLLEXPORT void
mm_set_iso_check (int iso_check)
{
  if (!yet_initialized)
    init_globals_dll ();
  opt_iso = (iso_check != FEATURE_OFF);
}

DLLEXPORT void
mm_set_rad_check (int rad_check)
{
  if (!yet_initialized)
    init_globals_dll ();
  opt_rad = (rad_check != FEATURE_OFF);
}

DLLEXPORT void
mm_set_exact_match (int exact)
{
  if (!yet_initialized)
    init_globals_dll ();
  opt_exact = (exact != FEATURE_OFF);
}

DLLEXPORT int
mm_match ()
{
  mol_count = 1;
//     mm_ElabMol;
/*printf("%i\n",opt_exact);
printf("%i\n",n_Ctot);
printf("%i\n",ndl_n_Ctot);
printf("%i\n",n_Otot);
printf("%i\n",ndl_n_Otot);
printf("%i\n",n_Ntot);
printf("%i\n",ndl_n_Ntot);*/
  if (opt_exact
      && (n_Ctot != ndl_n_Ctot || n_Otot != ndl_n_Otot
	  || n_Ntot != ndl_n_Ntot))
    return 0;
  init_molstat (&ndl_molstat);
  //perform_match ();
  //---------------------------------------------------- 0.3o
  if (rs_strict)		/* v0.3j */
    ndl_ref_atom = find_ndl_ref_atom_cv ();
  //ndl_ref_atom = find_ndl_ref_atom ();
  else
    ndl_ref_atom = find_ndl_ref_atom ();
  clear_atom_tags ();
  set_ndl_atom_tags ();
  matchsummary = true;
  perform_match ();
  matchsummary = matchresult;
  if (count_tagged_ndl_heavyatoms () > 0 && matchsummary == true)
    {
      do
	{
	  if (rs_strict)
	    ndl_ref_atom = find_ndl_ref_atom_cv ();
	  else
	    ndl_ref_atom = find_ndl_ref_atom ();
	  perform_match ();
	  if (matchresult == false)
	    matchsummary = false;
	}
      while (count_tagged_ndl_heavyatoms () != 0 && matchsummary != false);
    }

  //-----------------------------------------------------

  //mol_count = 0;

  //molbufindex = 0;

  //return matchresult ? 1 : 0;

  return matchsummary ? 1 : 0;
}

//-------------------------

DLLEXPORT void
xm_set_ring_perception_algorithm (int algo)
{
  switch (algo)
    {
    case RPA_SAR:
      opt_rs_dll = rs_sar;
      break;
    case RPA_SSR:
      opt_rs_dll = rs_ssr;
      break;
    default:
      opt_rs_dll = RPA_DEFAULT;
      break;
    }
  //printf("RPA_SET: %i\n",opt_rs_dll);
}

static void
write_molstat_X_dll (char *out_buffer)
{
  char tmp_buf[256];
  out_buffer[0] = '\0';
  if (auto_ssr)			/* v0.3n */
    fix_ssr_ringcounts ();
  sprintf (tmp_buf, "%d,", n_heavyatoms);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", n_heavybonds);
  strcat (out_buffer, tmp_buf);
#ifdef REDUCED_SAR
  sprintf (tmp_buf, "%d,", n_countablerings);
  strcat (out_buffer, tmp_buf);
#else
  sprintf (tmp_buf, "%d,", n_rings);
  strcat (out_buffer, tmp_buf);
#endif
  sprintf (tmp_buf, "%d,", molstat.n_QA);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_QB);
  strcat (out_buffer, tmp_buf);
  //if (opt_chg)
//    {                         /* 0.3x */
//      sprintf (tmp_buf, "%d,", molstat.n_chg);
//    }
//  else
//    {
  sprintf (tmp_buf, "%d,", molstat.n_chg);
  //   }
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_C1);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_C2);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_C);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_CHB1p);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_CHB2p);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_CHB3p);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_CHB4);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_O2);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_O3);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_N1);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_N2);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_N3);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_S);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_SeTe);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_F);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_Cl);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_Br);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_I);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_P);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_B);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_Met);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_X);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_b1);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_b2);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_b3);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_bar);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_C1O);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_C2O);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_CN);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_XY);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_r3);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_r4);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_r5);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_r6);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_r7);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_r8);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_r9);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_r10);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_r11);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_r12);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_r13p);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_rN);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_rN1);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_rN2);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_rN3p);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_rO);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_rO1);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_rO2p);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_rS);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_rX);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_rAr);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_rBz);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_br2p);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_psg01);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_psg02);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_psg13);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_psg14);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_psg15);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_psg16);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_psg17);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_psg18);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_pstm);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_psla);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d,", molstat.n_iso);
  strcat (out_buffer, tmp_buf);
  sprintf (tmp_buf, "%d", molstat.n_rad);
  strcat (out_buffer, tmp_buf);
}

static void
write_molstat_X_arr_dll (unsigned short *out_buffer)
{
  if (auto_ssr)			/* v0.3n */
    fix_ssr_ringcounts ();
  out_buffer[0] = n_heavyatoms;
  out_buffer[1] = n_heavybonds;
 
#ifdef REDUCED_SAR
  out_buffer[2] = n_countablerings;
#else
  out_buffer[2] = n_rings;
#endif

  out_buffer[3] = molstat.n_QA;
  
  out_buffer[4] = molstat.n_QB;
 
  out_buffer[5] = molstat.n_chg;
 
  out_buffer[6] = molstat.n_C1;
 
  out_buffer[7] = molstat.n_C2;
 
  out_buffer[8] = molstat.n_C;
 
  out_buffer[9] = molstat.n_CHB1p;
 
  out_buffer[10] = molstat.n_CHB2p;
 
  out_buffer[11] = molstat.n_CHB3p;
 
  out_buffer[12] = molstat.n_CHB4;
 
  out_buffer[13] = molstat.n_O2;
 
  out_buffer[14] = molstat.n_O3;
 
  out_buffer[15] = molstat.n_N1;
 
  out_buffer[16] = molstat.n_N2;
 
  out_buffer[17] = molstat.n_N3;
 
  out_buffer[18] = molstat.n_S;
 
  out_buffer[19] = molstat.n_SeTe;
 
  out_buffer[20] = molstat.n_F;
 
  out_buffer[21] = molstat.n_Cl;
 
  out_buffer[22] = molstat.n_Br;
 
  out_buffer[23] = molstat.n_I;
 
  out_buffer[24] = molstat.n_P;
 
  out_buffer[25] = molstat.n_B;
 
  out_buffer[26] = molstat.n_Met;
 
  out_buffer[27] = molstat.n_X;
 
  out_buffer[28] = molstat.n_b1;
 
  out_buffer[29] = molstat.n_b2;
 
  out_buffer[30] = molstat.n_b3;
 
  out_buffer[31] = molstat.n_bar;
 
  out_buffer[32] = molstat.n_C1O;
 
  out_buffer[33] = molstat.n_C2O;
 
  out_buffer[34] = molstat.n_CN;
 
  out_buffer[35] = molstat.n_XY;
 
  out_buffer[36] = molstat.n_r3;
 
  out_buffer[37] = molstat.n_r4;
 
  out_buffer[38] = molstat.n_r5;
 
  out_buffer[39] = molstat.n_r6;
 
  out_buffer[40] = molstat.n_r7;
 
  out_buffer[41] = molstat.n_r8;
 
  out_buffer[42] = molstat.n_r9;
 
  out_buffer[43] = molstat.n_r10;
 
  out_buffer[44] = molstat.n_r11;
 
  out_buffer[45] = molstat.n_r12;
 
  out_buffer[46] = molstat.n_r13p;
 
  out_buffer[47] = molstat.n_rN;
 
  out_buffer[48] = molstat.n_rN1;
 
  out_buffer[49] = molstat.n_rN2;
 
  out_buffer[50] = molstat.n_rN3p;
 
  out_buffer[51] = molstat.n_rO;
 
  out_buffer[52] = molstat.n_rO1;
 
  out_buffer[53] = molstat.n_rO2p;
 
  out_buffer[54] = molstat.n_rS;
 
  out_buffer[55] = molstat.n_rX;
 
  out_buffer[56] = molstat.n_rAr;
 
  out_buffer[57] = molstat.n_rBz;
 
  out_buffer[58] = molstat.n_br2p;
 
  out_buffer[59] = molstat.n_psg01;
 
  out_buffer[60] = molstat.n_psg02;
 
  out_buffer[61] = molstat.n_psg13;
 
  out_buffer[62] = molstat.n_psg14;
 
  out_buffer[63] = molstat.n_psg15;
 
  out_buffer[64] = molstat.n_psg16;
 
  out_buffer[65] = molstat.n_psg17;
 
  out_buffer[66] = molstat.n_psg18;
 
  out_buffer[67] = molstat.n_pstm;
 
  out_buffer[68] = molstat.n_psla;
 
  out_buffer[69] = molstat.n_iso;
 
  out_buffer[70] = molstat.n_rad;
 
}

static void
write_molstat_dll (char *out_buffer, int mode)
{
  char tmp_buf[256];
  char *sep1;
  char *sep2;
  switch (mode)
    {
    case 1:
      sep1 = "=";
      sep2 = " AND ";
      break;
    case 2:
      sep1 = "<=";
      sep2 = " AND ";
      break;
    default:
      sep1 = ":";
      sep2 = ";";
      break;
    }



  out_buffer[0] = '\0';

  if (auto_ssr)			/* v0.3n */
    fix_ssr_ringcounts ();
  sprintf (tmp_buf, "n_atoms%s%d%s", sep1, n_heavyatoms, sep2);
  strcat (out_buffer, tmp_buf);
  if (n_bonds > 0)
    {
      sprintf (tmp_buf, "n_bonds%s%d%s", sep1, n_heavybonds, sep2);
      strcat (out_buffer, tmp_buf);
    }

#ifdef REDUCED_SAR
  if (n_rings > 0)
    {
      sprintf (tmp_buf, "n_rings%s%d%s", sep1, n_countablerings, sep2);
      strcat (out_buffer, tmp_buf);
    }
#else
  if (n_rings > 0)
    {
      sprintf (tmp_buf, "n_rings%s%d%s", sep1, n_rings, sep2);
      strcat (out_buffer, tmp_buf);
    }
#endif

  if (opt_chg && molstat.n_chg > 0)	/* 0.3x */
    {
      sprintf (tmp_buf, "n_chg%s%d%s", sep1, molstat.n_chg, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_C1 > 0)
    {
      sprintf (tmp_buf, "n_C1%s%d%s", sep1, molstat.n_C1, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_C2 > 0)
    {
      sprintf (tmp_buf, "n_C2%s%d%s", sep1, molstat.n_C2, sep2);
      strcat (out_buffer, tmp_buf);
    }


  if (molstat.n_C > 0)
    {
      sprintf (tmp_buf, "n_C%s%d%s", sep1, molstat.n_C, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_CHB1p > 0)
    {
      sprintf (tmp_buf, "n_CHB1p%s%d%s", sep1, molstat.n_CHB1p, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_CHB2p > 0)
    {
      sprintf (tmp_buf, "n_CHB2p%s%d%s", sep1, molstat.n_CHB2p, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_CHB3p > 0)
    {
      sprintf (tmp_buf, "n_CHB3p%s%d%s", sep1, molstat.n_CHB3p, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_CHB4 > 0)
    {
      sprintf (tmp_buf, "n_CHB4%s%d%s", sep1, molstat.n_CHB4, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_O2 > 0)
    {
      sprintf (tmp_buf, "n_O2%s%d%s", sep1, molstat.n_O2, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_O3 > 0)
    {
      sprintf (tmp_buf, "n_O3%s%d%s", sep1, molstat.n_O3, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_N1 > 0)
    {
      sprintf (tmp_buf, "n_N1%s%d%s", sep1, molstat.n_N1, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_N2 > 0)
    {
      sprintf (tmp_buf, "n_N2%s%d%s", sep1, molstat.n_N2, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_N3 > 0)
    {
      sprintf (tmp_buf, "n_N3%s%d%s", sep1, molstat.n_N3, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_S > 0)
    {
      sprintf (tmp_buf, "n_S%s%d%s", sep1, molstat.n_S, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_SeTe > 0)
    {
      sprintf (tmp_buf, "n_SeTe%s%d%s", sep1, molstat.n_SeTe, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_F > 0)
    {
      sprintf (tmp_buf, "n_F%s%d%s", sep1, molstat.n_F, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_Cl > 0)
    {
      sprintf (tmp_buf, "n_Cl%s%d%s", sep1, molstat.n_Cl, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_Br > 0)
    {
      sprintf (tmp_buf, "n_Br%s%d%s", sep1, molstat.n_Br, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_I > 0)
    {
      sprintf (tmp_buf, "n_I%s%d%s", sep1, molstat.n_I, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_P > 0)
    {
      sprintf (tmp_buf, "n_P%s%d%s", sep1, molstat.n_P, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_B > 0)
    {
      sprintf (tmp_buf, "n_B%s%d%s", sep1, molstat.n_B, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_Met > 0)
    {
      sprintf (tmp_buf, "n_Met%s%d%s", sep1, molstat.n_Met, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_X > 0)
    {
      sprintf (tmp_buf, "n_X%s%d%s", sep1, molstat.n_X, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_b1 > 0)
    {
      sprintf (tmp_buf, "n_b1%s%d%s", sep1, molstat.n_b1, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_b2 > 0)
    {
      sprintf (tmp_buf, "n_b2%s%d%s", sep1, molstat.n_b2, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_b3 > 0)
    {
      sprintf (tmp_buf, "n_b3%s%d%s", sep1, molstat.n_b3, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_bar > 0)
    {
      sprintf (tmp_buf, "n_bar%s%d%s", sep1, molstat.n_bar, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_C1O > 0)
    {
      sprintf (tmp_buf, "n_C1O%s%d%s", sep1, molstat.n_C1O, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_C2O > 0)
    {
      sprintf (tmp_buf, "n_C2O%s%d%s", sep1, molstat.n_C2O, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_CN > 0)
    {
      sprintf (tmp_buf, "n_CN%s%d%s", sep1, molstat.n_CN, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_XY > 0)
    {
      sprintf (tmp_buf, "n_XY%s%d%s", sep1, molstat.n_XY, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_r3 > 0)
    {
      sprintf (tmp_buf, "n_r3%s%d%s", sep1, molstat.n_r3, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_r4 > 0)
    {
      sprintf (tmp_buf, "n_r4%s%d%s", sep1, molstat.n_r4, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_r5 > 0)
    {
      sprintf (tmp_buf, "n_r5%s%d%s", sep1, molstat.n_r5, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_r6 > 0)
    {
      sprintf (tmp_buf, "n_r6%s%d%s", sep1, molstat.n_r6, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_r7 > 0)
    {
      sprintf (tmp_buf, "n_r7%s%d%s", sep1, molstat.n_r7, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_r8 > 0)
    {
      sprintf (tmp_buf, "n_r8%s%d%s", sep1, molstat.n_r8, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_r9 > 0)
    {
      sprintf (tmp_buf, "n_r9%s%d%s", sep1, molstat.n_r9, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_r10 > 0)
    {
      sprintf (tmp_buf, "n_r10%s%d%s", sep1, molstat.n_r10, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_r11 > 0)
    {
      sprintf (tmp_buf, "n_r11%s%d%s", sep1, molstat.n_r11, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_r12 > 0)
    {
      sprintf (tmp_buf, "n_r12%s%d%s", sep1, molstat.n_r12, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_r13p > 0)
    {
      sprintf (tmp_buf, "n_r13p%s%d%s", sep1, molstat.n_r13p, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_rN > 0)
    {
      sprintf (tmp_buf, "n_rN%s%d%s", sep1, molstat.n_rN, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_rN1 > 0)
    {
      sprintf (tmp_buf, "n_rN1%s%d%s", sep1, molstat.n_rN1, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_rN2 > 0)
    {
      sprintf (tmp_buf, "n_rN2%s%d%s", sep1, molstat.n_rN2, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_rN3p > 0)
    {
      sprintf (tmp_buf, "n_rN3p%s%d%s", sep1, molstat.n_rN3p, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_rO > 0)
    {
      sprintf (tmp_buf, "n_rO%s%d%s", sep1, molstat.n_rO, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_rO1 > 0)
    {
      sprintf (tmp_buf, "n_rO1%s%d%s", sep1, molstat.n_rO1, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_rO2p > 0)
    {
      sprintf (tmp_buf, "n_rO2p%s%d%s", sep1, molstat.n_rO2p, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_rS > 0)
    {
      sprintf (tmp_buf, "n_rS%s%d%s", sep1, molstat.n_rS, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_rX > 0)
    {
      sprintf (tmp_buf, "n_rX%s%d%s", sep1, molstat.n_rX, sep2);
      strcat (out_buffer, tmp_buf);
    }
  if (molstat.n_rAr > 0)
    {
      sprintf (tmp_buf, "n_rar%s%d%s", sep1, molstat.n_rAr, sep2);
      strcat (out_buffer, tmp_buf);
    }

  if (molstat.n_rBz > 0)
    {
      sprintf (tmp_buf, "n_rbz%s%d%s", sep1, molstat.n_rBz, sep2);
      strcat (out_buffer, tmp_buf);
    }

  if (molstat.n_br2p > 0)
    {
      sprintf (tmp_buf, "n_br2p%s%d%s", sep1, molstat.n_br2p, sep2);
      strcat (out_buffer, tmp_buf);
    }

  if (molstat.n_psg01 > 0)
    {
      sprintf (tmp_buf, "n_psg01%s%d%s", sep1, molstat.n_psg01, sep2);
      strcat (out_buffer, tmp_buf);
    }

  if (molstat.n_psg02 > 0)
    {
      sprintf (tmp_buf, "n_psg02%s%d%s", sep1, molstat.n_psg02, sep2);
      strcat (out_buffer, tmp_buf);
    }

  if (molstat.n_psg13 > 0)
    {
      sprintf (tmp_buf, "n_psg13%s%d%s", sep1, molstat.n_psg13, sep2);
      strcat (out_buffer, tmp_buf);
    }

  if (molstat.n_psg14 > 0)
    {
      sprintf (tmp_buf, "n_psg14%s%d%s", sep1, molstat.n_psg14, sep2);
      strcat (out_buffer, tmp_buf);
    }

  if (molstat.n_psg15 > 0)
    {
      sprintf (tmp_buf, "n_psg15%s%d%s", sep1, molstat.n_psg15, sep2);
      strcat (out_buffer, tmp_buf);
    }

  if (molstat.n_psg16 > 0)
    {
      sprintf (tmp_buf, "n_psg16%s%d%s", sep1, molstat.n_psg16, sep2);
      strcat (out_buffer, tmp_buf);
    }

  if (molstat.n_psg17 > 0)
    {
      sprintf (tmp_buf, "n_psg17%s%d%s", sep1, molstat.n_psg17, sep2);
      strcat (out_buffer, tmp_buf);
    }

  if (molstat.n_psg18 > 0)
    {
      sprintf (tmp_buf, "n_psg18%s%d%s", sep1, molstat.n_psg18, sep2);
      strcat (out_buffer, tmp_buf);
    }

  if (molstat.n_pstm > 0)
    {
      sprintf (tmp_buf, "n_pstm%s%d%s", sep1, molstat.n_pstm, sep2);
      strcat (out_buffer, tmp_buf);
    }

  if (molstat.n_psla > 0)
    {
      sprintf (tmp_buf, "n_psla%s%d%s", sep1, molstat.n_psla, sep2);
      strcat (out_buffer, tmp_buf);
    }

  if (opt_iso && molstat.n_iso > 0)
    {
      sprintf (tmp_buf, "n_iso%s%d%s", sep1, molstat.n_iso, sep2);
      strcat (out_buffer, tmp_buf);
    }

  if (opt_rad && molstat.n_rad > 0)
    {
      sprintf (tmp_buf, "n_rad%s%d%s", sep1, molstat.n_rad, sep2);
      strcat (out_buffer, tmp_buf);
    }
}

static void
write_fg_code_dll (char *out_buffer)
{
  char tmp_buf[256];
  out_buffer[0] = '\0';
  if (fg[fg_cation - 1])
    {
      sprintf (tmp_buf, "000000T2;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_anion - 1])
    {
      sprintf (tmp_buf, "000000T1;");
      strcat (out_buffer, tmp_buf);
    }

  if (fg[fg_aldehyde - 1])
    {
      sprintf (tmp_buf, "C2O1H000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_ketone - 1])
    {
      sprintf (tmp_buf, "C2O1C000;");
      strcat (out_buffer, tmp_buf);
    }

  if (fg[fg_thioaldehyde - 1])
    {
      sprintf (tmp_buf, "C2S1H000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_thioketone - 1])
    {
      sprintf (tmp_buf, "C2S1C000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_imine - 1])
    {
      sprintf (tmp_buf, "C2N10000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_hydrazone - 1])
    {
      sprintf (tmp_buf, "C2N1N000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_semicarbazone - 1])
    {
      sprintf (tmp_buf, "C2NNC4ON;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_thiosemicarbazone - 1])
    {
      sprintf (tmp_buf, "C2NNC4SN;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_oxime - 1])
    {
      sprintf (tmp_buf, "C2N1OH00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_oxime_ether - 1])
    {
      sprintf (tmp_buf, "C2N1OC00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_ketene - 1])
    {
      sprintf (tmp_buf, "C3OC0000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_ketene_acetal_deriv - 1])
    {
      sprintf (tmp_buf, "C3OCC000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_carbonyl_hydrate - 1])
    {
      sprintf (tmp_buf, "C2O2H200;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_hemiacetal - 1])
    {
      sprintf (tmp_buf, "C2O2HC00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_acetal - 1])
    {
      sprintf (tmp_buf, "C2O2CC00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_hemiaminal - 1])
    {
      sprintf (tmp_buf, "C2NOHC10;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_aminal - 1])
    {
      sprintf (tmp_buf, "C2N2CC10;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_thiohemiaminal - 1])
    {
      sprintf (tmp_buf, "C2NSHC10;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_thioacetal - 1])
    {
      sprintf (tmp_buf, "C2S2CC00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_enamine - 1])
    {
      sprintf (tmp_buf, "C2CNH000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_enol - 1])
    {
      sprintf (tmp_buf, "C2COH000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_enolether - 1])
    {
      sprintf (tmp_buf, "C2COC000;");
      strcat (out_buffer, tmp_buf);
    }


  if (fg[fg_prim_alcohol - 1])
    {
      sprintf (tmp_buf, "O1H1C000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_sec_alcohol - 1])
    {
      sprintf (tmp_buf, "O1H2C000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_tert_alcohol - 1])
    {
      sprintf (tmp_buf, "O1H3C000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_1_2_diol - 1])
    {
      sprintf (tmp_buf, "O1H0CO1H;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_1_2_aminoalcohol - 1])
    {
      sprintf (tmp_buf, "O1H0CN1C;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_phenol - 1])
    {
      sprintf (tmp_buf, "O1H1A000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_1_2_diphenol - 1])
    {
      sprintf (tmp_buf, "O1H2A000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_enediol - 1])
    {
      sprintf (tmp_buf, "C2COH200;");
      strcat (out_buffer, tmp_buf);
    }

  if (fg[fg_dialkylether - 1])
    {
      sprintf (tmp_buf, "O1C0CC00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_alkylarylether - 1])
    {
      sprintf (tmp_buf, "O1C0CA00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_diarylether - 1])
    {
      sprintf (tmp_buf, "O1C0AA00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_thioether - 1])
    {
      sprintf (tmp_buf, "S1C00000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_disulfide - 1])
    {
      sprintf (tmp_buf, "S1S1C000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_peroxide - 1])
    {
      sprintf (tmp_buf, "O1O1C000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_hydroperoxide - 1])
    {
      sprintf (tmp_buf, "O1O1H000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_hydrazine - 1])
    {
      sprintf (tmp_buf, "N1N10000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_hydroxylamine - 1])
    {
      sprintf (tmp_buf, "N1O1H000;");
      strcat (out_buffer, tmp_buf);
    }


  if (fg[fg_prim_aliph_amine - 1])
    {
      sprintf (tmp_buf, "N1C1C000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_prim_arom_amine - 1])
    {
      sprintf (tmp_buf, "N1C1A000;");
      strcat (out_buffer, tmp_buf);
    }

  if (fg[fg_sec_aliph_amine - 1])
    {
      sprintf (tmp_buf, "N1C2CC00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_sec_mixed_amine - 1])
    {
      sprintf (tmp_buf, "N1C2AC00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_sec_arom_amine - 1])
    {
      sprintf (tmp_buf, "N1C2AA00;");
      strcat (out_buffer, tmp_buf);
    }

  if (fg[fg_tert_aliph_amine - 1])
    {
      sprintf (tmp_buf, "N1C3CC00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_tert_mixed_amine - 1])
    {
      sprintf (tmp_buf, "N1C3AC00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_tert_arom_amine - 1])
    {
      sprintf (tmp_buf, "N1C3AA00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_quart_ammonium - 1])
    {
      sprintf (tmp_buf, "N1C400T2;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_n_oxide - 1])
    {
      sprintf (tmp_buf, "N0O10000;");
      strcat (out_buffer, tmp_buf);
    }


  if (fg[fg_halogen_deriv - 1])
    {
      if (!fg[fg_alkyl_halide - 1] && !fg[fg_aryl_halide - 1] &&
	  !fg[fg_acyl_halide - 1])
	{
	  sprintf (tmp_buf, "XX000000;");
	  strcat (out_buffer, tmp_buf);
	}
    }

  if (fg[fg_alkyl_fluoride - 1])
    {
      sprintf (tmp_buf, "XF00C000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_alkyl_chloride - 1])
    {
      sprintf (tmp_buf, "XC00C000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_alkyl_bromide - 1])
    {
      sprintf (tmp_buf, "XB00C000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_alkyl_iodide - 1])
    {
      sprintf (tmp_buf, "XI00C000;");
      strcat (out_buffer, tmp_buf);
    }

  if (fg[fg_aryl_fluoride - 1])
    {
      sprintf (tmp_buf, "XF00A000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_aryl_chloride - 1])
    {
      sprintf (tmp_buf, "XC00A000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_aryl_bromide - 1])
    {
      sprintf (tmp_buf, "XB00A000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_aryl_iodide - 1])
    {
      sprintf (tmp_buf, "XI00A000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_organometallic - 1])
    {
      sprintf (tmp_buf, "000000MX;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_organolithium - 1])
    {
      sprintf (tmp_buf, "000000ML;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_organomagnesium - 1])
    {
      sprintf (tmp_buf, "000000MM;");
      strcat (out_buffer, tmp_buf);
    }

  if (fg[fg_carboxylic_acid - 1])
    {
      sprintf (tmp_buf, "C3O2H000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_carboxylic_acid_salt - 1])
    {
      sprintf (tmp_buf, "C3O200T1;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_carboxylic_acid_ester - 1])
    {
      sprintf (tmp_buf, "C3O2C000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_lactone - 1])
    {
      sprintf (tmp_buf, "C3O2CZ00;");
      strcat (out_buffer, tmp_buf);
    }

  if (fg[fg_carboxylic_acid_prim_amide - 1])
    {
      sprintf (tmp_buf, "C3ONC100;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_carboxylic_acid_sec_amide - 1])
    {
      sprintf (tmp_buf, "C3ONC200;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_carboxylic_acid_tert_amide - 1])
    {
      sprintf (tmp_buf, "C3ONC300;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_lactam - 1])
    {
      sprintf (tmp_buf, "C3ONCZ00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_carboxylic_acid_hydrazide - 1])
    {
      sprintf (tmp_buf, "C3ONN100;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_carboxylic_acid_azide - 1])
    {
      sprintf (tmp_buf, "C3ONN200;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_hydroxamic_acid - 1])
    {
      sprintf (tmp_buf, "C3ONOH00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_carboxylic_acid_amidine - 1])
    {
      sprintf (tmp_buf, "C3N2H000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_carboxylic_acid_amidrazone - 1])
    {
      sprintf (tmp_buf, "C3NNN100;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_nitrile - 1])
    {
      sprintf (tmp_buf, "C3N00000;");
      strcat (out_buffer, tmp_buf);
    }

  if (fg[fg_acyl_fluoride - 1])
    {
      sprintf (tmp_buf, "C3OXF000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_acyl_chloride - 1])
    {
      sprintf (tmp_buf, "C3OXC000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_acyl_bromide - 1])
    {
      sprintf (tmp_buf, "C3OXB000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_acyl_iodide - 1])
    {
      sprintf (tmp_buf, "C3OXI000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_acyl_cyanide - 1])
    {
      sprintf (tmp_buf, "C2OC3N00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_imido_ester - 1])
    {
      sprintf (tmp_buf, "C3NOC000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_imidoyl_halide - 1])
    {
      sprintf (tmp_buf, "C3NXX000;");
      strcat (out_buffer, tmp_buf);
    }

  if (fg[fg_thiocarboxylic_acid - 1])
    {
      sprintf (tmp_buf, "C3SOH000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_thiocarboxylic_acid_ester - 1])
    {
      sprintf (tmp_buf, "C3SOC000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_thiolactone - 1])
    {
      sprintf (tmp_buf, "C3SOCZ00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_thiocarboxylic_acid_amide - 1])
    {
      sprintf (tmp_buf, "C3SNH000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_thiolactam - 1])
    {
      sprintf (tmp_buf, "C3SNCZ00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_imido_thioester - 1])
    {
      sprintf (tmp_buf, "C3NSC000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_oxohetarene - 1])
    {
      sprintf (tmp_buf, "C3ONAZ00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_thioxohetarene - 1])
    {
      sprintf (tmp_buf, "C3SNAZ00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_iminohetarene - 1])
    {
      sprintf (tmp_buf, "C3NNAZ00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_orthocarboxylic_acid_deriv - 1])
    {
      sprintf (tmp_buf, "C3O30000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_carboxylic_acid_orthoester - 1])
    {
      sprintf (tmp_buf, "C3O3C000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_carboxylic_acid_amide_acetal - 1])
    {
      sprintf (tmp_buf, "C3O3NC00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_carboxylic_acid_anhydride - 1])
    {
      sprintf (tmp_buf, "C3O2C3O2;");
      strcat (out_buffer, tmp_buf);
    }

  if (fg[fg_carboxylic_acid_unsubst_imide - 1])
    {
      sprintf (tmp_buf, "C3ONCH10;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_carboxylic_acid_subst_imide - 1])
    {
      sprintf (tmp_buf, "C3ONCC10;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_co2_deriv - 1])
    {
      sprintf (tmp_buf, "C4000000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_carbonic_acid_deriv - 1])
    {
      sprintf (tmp_buf, "C4O30000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_carbonic_acid_monoester - 1])
    {
      sprintf (tmp_buf, "C4O3C100;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_carbonic_acid_diester - 1])
    {
      sprintf (tmp_buf, "C4O3C200;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_carbonic_acid_ester_halide - 1])
    {
      sprintf (tmp_buf, "C4O3CX00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_thiocarbonic_acid_deriv - 1])
    {
      sprintf (tmp_buf, "C4SO0000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_thiocarbonic_acid_monoester - 1])
    {
      sprintf (tmp_buf, "C4SOC100;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_thiocarbonic_acid_diester - 1])
    {
      sprintf (tmp_buf, "C4SOC200;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_thiocarbonic_acid_ester_halide - 1])
    {
      sprintf (tmp_buf, "C4SOX_00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_carbamic_acid_deriv - 1])
    {
      sprintf (tmp_buf, "C4O2N000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_carbamic_acid - 1])
    {
      sprintf (tmp_buf, "C4O2NH00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_carbamic_acid_ester - 1])
    {
      sprintf (tmp_buf, "C4O2NC00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_carbamic_acid_halide - 1])
    {
      sprintf (tmp_buf, "C4O2NX00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_thiocarbamic_acid_deriv - 1])
    {
      sprintf (tmp_buf, "C4SN0000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_thiocarbamic_acid - 1])
    {
      sprintf (tmp_buf, "C4SNOH00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_thiocarbamic_acid_ester - 1])
    {
      sprintf (tmp_buf, "C4SNOC00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_thiocarbamic_acid_halide - 1])
    {
      sprintf (tmp_buf, "C4SNXX00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_urea - 1])
    {
      sprintf (tmp_buf, "C4O1N200;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_isourea - 1])
    {
      sprintf (tmp_buf, "C4N2O100;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_thiourea - 1])
    {
      sprintf (tmp_buf, "C4S1N200;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_isothiourea - 1])
    {
      sprintf (tmp_buf, "C4N2S100;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_guanidine - 1])
    {
      sprintf (tmp_buf, "C4N30000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_semicarbazide - 1])
    {
      sprintf (tmp_buf, "C4ON2N00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_thiosemicarbazide - 1])
    {
      sprintf (tmp_buf, "C4SN2N00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_azide - 1])
    {
      sprintf (tmp_buf, "N4N20000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_azo_compound - 1])
    {
      sprintf (tmp_buf, "N2N10000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_diazonium_salt - 1])
    {
      sprintf (tmp_buf, "N3N100T2;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_isonitrile - 1])
    {
      sprintf (tmp_buf, "N3C10000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_cyanate - 1])
    {
      sprintf (tmp_buf, "C4NO1000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_isocyanate - 1])
    {
      sprintf (tmp_buf, "C4NO2000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_thiocyanate - 1])
    {
      sprintf (tmp_buf, "C4NS1000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_isothiocyanate - 1])
    {
      sprintf (tmp_buf, "C4NS2000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_carbodiimide - 1])
    {
      sprintf (tmp_buf, "C4N20000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_nitroso_compound - 1])
    {
      sprintf (tmp_buf, "N2O10000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_nitro_compound - 1])
    {
      sprintf (tmp_buf, "N4O20000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_nitrite - 1])
    {
      sprintf (tmp_buf, "N3O20000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_nitrate - 1])
    {
      sprintf (tmp_buf, "N4O30000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_sulfuric_acid_deriv - 1])
    {
      sprintf (tmp_buf, "S6O00000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_sulfuric_acid - 1])
    {
      sprintf (tmp_buf, "S6O4H000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_sulfuric_acid_monoester - 1])
    {
      sprintf (tmp_buf, "S6O4HC00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_sulfuric_acid_diester - 1])
    {
      sprintf (tmp_buf, "S6O4CC00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_sulfuric_acid_amide_ester - 1])
    {
      sprintf (tmp_buf, "S6O3NC00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_sulfuric_acid_amide - 1])
    {
      sprintf (tmp_buf, "S6O3N100;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_sulfuric_acid_diamide - 1])
    {
      sprintf (tmp_buf, "S6O2N200;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_sulfuryl_halide - 1])
    {
      sprintf (tmp_buf, "S6O3XX00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_sulfonic_acid_deriv - 1])
    {
      sprintf (tmp_buf, "S5O00000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_sulfonic_acid - 1])
    {
      sprintf (tmp_buf, "S5O3H000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_sulfonic_acid_ester - 1])
    {
      sprintf (tmp_buf, "S5O3C000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_sulfonamide - 1])
    {
      sprintf (tmp_buf, "S5O2N000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_sulfonyl_halide - 1])
    {
      sprintf (tmp_buf, "S5O2XX00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_sulfone - 1])
    {
      sprintf (tmp_buf, "S4O20000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_sulfoxide - 1])
    {
      sprintf (tmp_buf, "S2O10000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_sulfinic_acid_deriv - 1])
    {
      sprintf (tmp_buf, "S3O00000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_sulfinic_acid - 1])
    {
      sprintf (tmp_buf, "S3O2H000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_sulfinic_acid_ester - 1])
    {
      sprintf (tmp_buf, "S3O2C000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_sulfinic_acid_halide - 1])
    {
      sprintf (tmp_buf, "S3O1XX00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_sulfinic_acid_amide - 1])
    {
      sprintf (tmp_buf, "S3O1N000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_sulfenic_acid_deriv - 1])
    {
      sprintf (tmp_buf, "S1O00000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_sulfenic_acid - 1])
    {
      sprintf (tmp_buf, "S1O1H000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_sulfenic_acid_ester - 1])
    {
      sprintf (tmp_buf, "S1O1C000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_sulfenic_acid_halide - 1])
    {
      sprintf (tmp_buf, "S1O0XX00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_sulfenic_acid_amide - 1])
    {
      sprintf (tmp_buf, "S1O0N100;");
      strcat (out_buffer, tmp_buf);
    }

  if (fg[fg_alkylthiol - 1])
    {
      sprintf (tmp_buf, "S1H1C000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_arylthiol - 1])
    {
      sprintf (tmp_buf, "S1H1A000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_phosphoric_acid_deriv - 1])
    {
      sprintf (tmp_buf, "P5O0H000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_phosphoric_acid - 1])
    {
      sprintf (tmp_buf, "P5O4H200;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_phosphoric_acid_ester - 1])
    {
      sprintf (tmp_buf, "P5O4HC00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_phosphoric_acid_halide - 1])
    {
      sprintf (tmp_buf, "P5O3HX00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_phosphoric_acid_amide - 1])
    {
      sprintf (tmp_buf, "P5O3HN00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_thiophosphoric_acid_deriv - 1])
    {
      sprintf (tmp_buf, "P5O0S000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_thiophosphoric_acid - 1])
    {
      sprintf (tmp_buf, "P5O3SH00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_thiophosphoric_acid_ester - 1])
    {
      sprintf (tmp_buf, "P5O3SC00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_thiophosphoric_acid_halide - 1])
    {
      sprintf (tmp_buf, "P5O2SX00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_thiophosphoric_acid_amide - 1])
    {
      sprintf (tmp_buf, "P5O2SN00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_phosphonic_acid_deriv - 1])
    {
      sprintf (tmp_buf, "P4O30000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_phosphonic_acid - 1])
    {
      sprintf (tmp_buf, "P4O3H000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_phosphonic_acid_ester - 1])
    {
      sprintf (tmp_buf, "P4O3C000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_phosphine - 1])
    {
      sprintf (tmp_buf, "P3000000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_phosphinoxide - 1])
    {
      sprintf (tmp_buf, "P2O00000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_boronic_acid_deriv - 1])
    {
      sprintf (tmp_buf, "B2O20000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_boronic_acid - 1])
    {
      sprintf (tmp_buf, "B2O2H000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_boronic_acid_ester - 1])
    {
      sprintf (tmp_buf, "B2O2C000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_alkene - 1])
    {
      sprintf (tmp_buf, "000C2C00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_alkyne - 1])
    {
      sprintf (tmp_buf, "000C3C00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_aromatic - 1])
    {
      sprintf (tmp_buf, "0000A000;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_heterocycle - 1])
    {
      sprintf (tmp_buf, "0000CZ00;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_alpha_aminoacid - 1])
    {
      sprintf (tmp_buf, "C3O2HN1C;");
      strcat (out_buffer, tmp_buf);
    }
  if (fg[fg_alpha_hydroxyacid - 1])
    {
      sprintf (tmp_buf, "C3O2HO1H;");
      strcat (out_buffer, tmp_buf);
    }
}

/*static void mm_elab_molstat(void)
{
		count_neighbors();
//	init_molstat(&molstat);
	if (!found_arominfo) {
    chk_ringbonds();
    if (ringsearch_mode == rs_ssr)
      remove_redundant_rings();
    if (n_rings == max_rings) {
      if (opt_verbose)
	printf("warning: max. number of rings exceeded, reverting to SSR search\n");
      ringsearch_mode = rs_ssr;
      clear_rings();
      max_vringsize = 10;
      chk_ringbonds();
      remove_redundant_rings();
    }
    
    update_ringcount();
    update_atypes();
    update_Htotal();
    chk_arom();

    if (ringsearch_mode == rs_ssr) {
      do {
	prev_n_ar = count_aromatic_rings();
	chk_arom();
	n_ar = count_aromatic_rings();
      } while (prev_n_ar - n_ar != 0);
    }
	} else {
 update_atypes();
    update_Htotal();  
		
	}
 
 // get_molstat();
}*/

DLLEXPORT void
cm_molstat_X (char *buf)
{
  init_molstat (&molstat);
//mm_elab_molstat();
  get_molstat ();
  write_molstat_X_dll (buf);
}

DLLEXPORT void
cm_molstat_X_arr (unsigned short *buf)
{
  init_molstat (&molstat);
//mm_elab_molstat();
  get_molstat ();
  write_molstat_X_arr_dll (buf);
}

DLLEXPORT void
cm_molstat (char *buf)
{
  init_molstat (&molstat);
//mm_elab_molstat();
  get_molstat ();
  write_molstat_dll (buf, 0);
}

DLLEXPORT void
cm_molstat_sql_exact (char *buf)
{
  init_molstat (&molstat);
//mm_elab_molstat();
  get_molstat ();
  write_molstat_dll (buf, 1);
}

DLLEXPORT void
cm_molstat_sql_substruct (char *buf)
{
  init_molstat (&molstat);
//mm_elab_molstat();
  get_molstat ();
  write_molstat_dll (buf, 2);
}

DLLEXPORT void
cm_fg_codes (char *buf)
{
//mm_elab_molstat();
  chk_functionalgroups ();
  write_fg_code_dll (buf);
}

static void
write_MDLmolfile_dll (char *out_buffer)
{
  int i;
  char tmpstr[256];
  char wline[256];
  int a_chg;
  int a_iso;
  int a_rad;
  char tmflabel[256];		/* v0.3m */
  char STR1[256], STR7[256];
  int FORLIM;
  *out_buffer = '\0';
  *tmpstr = '\0';
  *wline = '\0';
  sprintf (tmflabel, "%i", tweaklevel);	/* v0.3m */
  while (strlen (tmflabel) < 2)	/* v0.3m */
    sprintf (tmflabel, "0%s", strcpy (STR1, tmflabel));
  sprintf (tmflabel, "TMF%s", strcpy (STR1, tmflabel));	/* v0.3m */
  if (strlen (molname) > 80)
    sprintf (molname, "%.80s", strcpy (STR1, molname));
  strncat (out_buffer, molname, 80);
  sprintf (wline, "\n  CheckMol                        %s", tmflabel);	/* v0.3m */
  if (ringsearch_mode == rs_sar)	/* v0.3m */
    strcat (wline, ":r0");
  if (ringsearch_mode == rs_ssr)	/* v0.3m */
    strcat (wline, ":r1");
  if (opt_metalrings)
    strcat (wline, ":m1");
  else
    strcat (wline, ":m0");
  /* v0.3m */
  sprintf (tmpstr, "\n%s\n", molcomment);
  strcat (wline, tmpstr);
  sprintf (tmpstr, "%d", n_atoms);
  lblank (3L, tmpstr);
  strcat (wline, tmpstr);
  /* first 3 digits: number of atoms */
  sprintf (tmpstr, "%d", n_bonds);
  lblank (3L, tmpstr);
  strcat (wline, tmpstr);
  *tmpstr = '\0';		/* next 3 digits: number of bonds */
  strcpy (tmpstr, "  0");
  strcat (wline, tmpstr);
  /* next 3 digits: number of atom lists (not used by us) */
/* p2c: checkmol.pas, line 2388:
 * Note: Turbo Pascal conditional compilation directive was ignored [218] */
#ifdef REDUCED_SAR
  sprintf (tmpstr, "%d", n_countablerings);
  /* v0.3n; changed n_rings into n_countablerings */
#else
  sprintf (tmpstr, "%d", n_rings);
#endif
  lblank (3L, tmpstr);
  strcat (wline, tmpstr);
  /* officially "obsolete", we use it for the number of rings */
  strcat (wline, "  ");		/* v0.3n: obey chiral flag */
  if (chir_flag)
    strcat (wline, "1");
  else
    strcat (wline, "0");
  /* v0.3n */
  strcat (wline, "               999 V2000\n");
  /* v0.3n (adjust string length) */
  strcat (out_buffer, wline);
  FORLIM = n_atoms;
  for (i = 0; i < FORLIM; i++)
    {
      *wline = '\0';
      *tmpstr = '\0';
      sprintf (tmpstr, "%1.4f", atom[i].x);
      lblank (10L, tmpstr);
      strcat (wline, tmpstr);
      sprintf (tmpstr, "%1.4f", atom[i].y);
      lblank (10L, tmpstr);
      strcat (wline, tmpstr);
      sprintf (tmpstr, "%1.4f", atom[i].z);
      lblank (10L, tmpstr);
      strcat (wline, tmpstr);
      strcpy (tmpstr, atom[i].element);
      /* tmpstr := lowercase(tmpstr); REPLACE!!! */
      //tmpstr[0] = toupper (tmpstr[0]);
      all_lowercase (tmpstr);
      tmpstr[0] = toupper (tmpstr[0]);
      /*wline := wline + ' '+atom^[i].element+' '; */
      sprintf (wline + strlen (wline), " %s ", tmpstr);
      strcat (wline, " 0");	/* mass difference (isotopes) */
      /* now we code aromaticity into the old-style charge column (charges are now in the M  CHG line) */
      if (atom[i].arom)
	strcpy (tmpstr, " 00");
      else
	strcpy (tmpstr, "  0");
      strcat (wline, tmpstr);
      strcat (wline, "  0  0  0  0  0  0  0  0  0  0\n");
      strcat (out_buffer, wline);
    }
  FORLIM = n_bonds;
  for (i = 0; i < FORLIM; i++)
    {
      *wline = '\0';
      *tmpstr = '\0';
      sprintf (tmpstr, "%d", bond[i].a1);
      lblank (3L, tmpstr);
      strcat (wline, tmpstr);
      sprintf (tmpstr, "%d", bond[i].a2);
      lblank (3L, tmpstr);
      strcat (wline, tmpstr);
      if (bond[i].btype == 'S')
	strcpy (tmpstr, "  1");
      if (bond[i].btype == 'D')
	strcpy (tmpstr, "  2");
      if (bond[i].btype == 'T')
	strcpy (tmpstr, "  3");
      if (bond[i].btype == 'A')
	strcpy (tmpstr, "  4");
      if (bond[i].btype == 'l')
	strcpy (tmpstr, "  5");
      if (bond[i].btype == 's')
	strcpy (tmpstr, "  6");
      if (bond[i].btype == 'd')
	strcpy (tmpstr, "  7");
      if (bond[i].btype == 'a')
	strcpy (tmpstr, "  8");
      /* now encode our own aromaticity information */
      if (bond[i].arom)
	tmpstr[1] = '0';
      strcat (wline, tmpstr);	/* next, encode bond stereo property (v0.3f) */
      /*if (bond^[i].stereo = bstereo_up) then wline := wline + '  1' else */
      /*  if (bond^[i].stereo = bstereo_down) then wline := wline + '  6' else */
      /*    wline := wline + '  0'; */
      /* restore original value from MDL molfile (v0.3n) */
      /* wline := wline + '  ' + inttostr(bond^[i].mdl_stereo);    REPLACE!!! */
      *tmpstr = '\0';
      sprintf (tmpstr, "%i", bond[i].mdl_stereo);
      strcat (wline, "  ");
      strcat (wline, tmpstr);
      *tmpstr = '\0';
      /* now encode the ring_count of this bond (using a field which officially is "not used") */
      /* tmpstr := inttostr(bond^[i].ring_count); REPLACE!!! */
      sprintf (tmpstr, "%i", bond[i].ring_count);
      while (strlen (tmpstr) < 3)
	sprintf (tmpstr, " %s", strcpy (STR1, tmpstr));
      sprintf (wline + strlen (wline), "%s  0  0\n", tmpstr);
      strcat (out_buffer, wline);
    }
  FORLIM = n_atoms;
  for (i = 1; i <= FORLIM; i++)
    {
      a_chg = atom[i - 1].formal_charge;
      if (a_chg != 0)
	{
	  strcpy (wline, "M  CHG  1 ");
	  sprintf (tmpstr, "%d", i);
	  lblank (3L, tmpstr);
	  sprintf (wline + strlen (wline), "%s ", tmpstr);
	  sprintf (tmpstr, "%d", a_chg);
	  lblank (3L, tmpstr);
	  strcat (wline, tmpstr);
	  strcat (out_buffer, wline);
	  strcat (out_buffer, "\n");
	}
    }
  for (i = 1; i <= FORLIM; i++)	/* 0.3x */
    {
      a_iso = atom[i - 1].nucleon_number;
      if (a_iso != 0)
	{
	  strcpy (wline, "M  ISO  1 ");
	  sprintf (tmpstr, "%d", i);
	  lblank (3L, tmpstr);
	  sprintf (wline + strlen (wline), "%s ", tmpstr);
	  sprintf (tmpstr, "%d", a_iso);
	  lblank (3L, tmpstr);
	  strcat (wline, tmpstr);
	  strcat (out_buffer, wline);
	  strcat (out_buffer, "\n");
	}
    }
  for (i = 1; i <= FORLIM; i++)	/* 0.3x */
    {
      a_rad = atom[i - 1].radical_type;
      if (a_rad != 0)
	{
	  strcpy (wline, "M  RAD  1 ");
	  sprintf (tmpstr, "%d", i);
	  lblank (3L, tmpstr);
	  sprintf (wline + strlen (wline), "%s ", tmpstr);
	  sprintf (tmpstr, "%d", a_rad);
	  lblank (3L, tmpstr);
	  strcat (wline, tmpstr);
	  strcat (out_buffer, wline);
	  strcat (out_buffer, "\n");
	}
    }
  strcat (out_buffer, "M  END\n");
}

DLLEXPORT void
cm_tweak_molfile (char *buf)
{
//chk_functionalgroups();
  write_MDLmolfile_dll (buf);
}
#endif
