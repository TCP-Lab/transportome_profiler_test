
???
Blood_cancer
Bone_marrow_cancer
???

trasformare 'Colorectal_cancer' in 'Colon cancer' o aggiungere 'Rectum' di TCGA ai case  ???

sample_type... verificare (tenere i metastatici ??)


non ci sono campioni Vagina in TCGA (manca il case)... eliminare il type 'Vaginal_cancer'





TCGA
9                                            Soft tissue,Bone
22                                           Lymphatic tissue
29                                           White blood cell


GTEX
1                                                 Bone Marrow
8                                                       Blood
22                                                     Spleen




https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html#:~:text=To%20perform%20the%20median%20of,which%20we%20will%20see%20later.


** fold_change
** norm_fold_change

deseq_shrinkage

** cohen_d
** norm_cohen_d

** s2n_ratio
** norm_s2n_ratio

-- bws_test
-- norm_bws_test

---------------
Rb+ nella heatmap finale non va bene... dai!

da dove è stata presa l'informazione sul pore loops e pore forming (i.e., la tabella 'structure' dell' MTP-DB)?



perché porre la limitazione di 40 geni per la ricorsione sulle colonne successivedella tabella?
previene in generale la generazione di liste valide (geni > 10)...



 "^whole_transportome///pores///channels///gating_mechanism::voltage///carried_solute::Rb\\+///pore_forming::1.0$"

ENSG00000115474 compare a sproposito perché ha
gating_mechanism == NULL

"^whole_transportome///pores///channels///gating_mechanism::voltage///carried_solute::cation///pore_forming::0.0$"

ha molti canali non corretti

--------------

----------------------
make_genesets.py -> make_genesets_mod.py
----------------------

Misspellings:
metasample.py @ line 170
	raise ValueError("Row name column not fonud in metadata")
	raise ValueError("Row name column not found in metadata")

Update comment:
select_and_run.py @ line 70
	# Now we can run run_deseq.R

Update help:
ranking_methods.py (or generanker)
`s2n_ratio` description (desc = "Use the signal to noise ratio (diff of means divided by variance)") or standard deviation?

riscrivere il JSON per
lavorare a livello di studio TCGA e non a livello di primary site

riflettere sul comportamento delle metriche di ranking per geni assenti (righe nulle)



nel database 
amino phospolipid
credo sia un refuso



Add this:
```python
import sys

def exit_program(string):
    print(string)
    sys.exit(0)
```



# 1. Generate large tables
uses `make_large_tables()` to make 9 large tables based on hardcoded `basic_gene_lists` found in
`./data/in/basic_gene_lists.json`
```
whole_transportome
|___pores
|	|___channels
|	|___aquaporins
|___transporters
	|___solute_carriers
	|___atp_driven
		|___ABC
		|___pumps
```
```python
	# Check the size of tables
	for name, table in large_tables.items():
	    print(name + ":  " + str(table.shape))
	exit_program("\nExit_1\n")

	# 2. Generate lists from large tables
```
# 2. Generate lists from large tables
For each large table, use `bonsai` to generate tree structures of all the possible gene sets, based on the 3 parameters of the function `generate_gene_list_trees()`, with the following default values:
```
min_pop_score: float = 0.5,
min_set_size: int = 10,
min_recurse_set_size: int = 40,
```

```python
	# Manually check these two tables
	trees["solute_carriers"].to_representation(sys.stdout)
	trees["pumps"].to_representation(sys.stdout)
	exit_program("\nExit_2\n")

	# 3. Make the union of the genesets following the structure
```

# 3. Make the union of the genesets following the structure
Paste all the trees togheter
```python
# Visualize the assembled tree and check the total size
large_tree.to_representation(sys.stdout)
print("\n")
for name, table in large_tables.items():
    print(trees[name])

print("--- TOTAL ---")
print(large_tree)
exit_program("\nExit_3.1\n")

if not args.no_prune:
```
