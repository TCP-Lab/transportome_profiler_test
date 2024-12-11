# How to run the ___Transportome Profiler___ locally

## Clone the repo
For any Linux system, including WSL:
```bash
git clone git@github.com:TCP-Lab/transportome_profiler.git
cd ./transportome_profiler
```

## Install R requirements
```bash
Rscript ./src/helper_scripts/install_r_pkgs.R
```

## (Create and) Activate a Python virtual environment and install requirements
```bash
python -m venv env
source ./env/bin/activate
pip install -r ./src/requirements.txt

# To update internal packages (i.e., bonsai, gene_ranker, metasplit, and panid)
# from the respective git repos  
pip install --force-reinstall -r ./src/requirements.txt

# When finished:
deactivate
```
More on Virtual Environments [here](https://docs.python.org/3/tutorial/venv.html).

## Install `xsv`
This is a program required by `metasplit` for the fast reshaping of very large
CSV files.
```bash
sudo pacman -Syu xsv
```

## Install `fast-cohen`
This is a program written in __Rust__ that performs a fast computation of the
Cohen's _d_ statistics.
```bash
cargo install --git https://github.com/MrHedmad/fast-cohen.git
```

## Install `Kerblam!`
This is our workflow manager.
```bash
# Install a prebuilt binary
curl --proto '=https' --tlsv1.2 -LsSf https://github.com/MrHedmad/kerblam/releases/latest/download/kerblam-installer.sh | sh

# Or, alternatively, install it from source with cargo
cargo install kerblam

# Check it with
kerblam --version
```

## Fetch MTP-DB, TCGA and GTEx datasets from remote (through XENA)
```bash
kerblam data fetch
```
will download in `/data/in` the following objects:
1. `expression_matrix.tsv.gz`, an archive containing all the log2 raw counts for
	transcript abundance quantification, as provided by TCGA and GTEx projects;
1. `expression_matrix_metadata.tsv.gz`, an archive containing the related
	metadata for each sample (i.e., patient or, better, specimen);
1. `MTPDB.sqlite.gz`, an archive containing our
	[Membrane Transport Protein Database](https://github.com/TCP-Lab/MTP-DB)
	used for gene set definition;
1. `ensg_data.csv`, giving for each ENSG ID the corresponding HGNC gene symbol; 
1. `geo`, a folder containing the tables of raw counts and metadata for
	individual small studies retrieved from GEO for validation purposes.

### Raw counts
The `expression_matrix.tsv.gz` archive is our Zenodo copy of the
_gene expression RNAseq - RSEM expected_count_ data set by
[XENA](https://xenabrowser.net/datapages/?dataset=TcgaTargetGtex_gene_expected_count&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443),
containing both TCGA and GTEx _harmonized_ transcriptomics data.
- author: _UCSC TOIL RNA-seq recompute_
- unit: __log2(expected_count+1)__
- size: 60,499 identifiers (genes) x 19,109 samples
- download: https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_gene_expected_count.gz

#### Sample check
```bash
zcat expression_matrix.tsv.gz | head -n1 | grep -oE 'TCGA-|GTEX-|TARGET-|\sK-' | wc -l
zcat expression_matrix.tsv.gz | head -n1 | grep -oE 'TCGA-' | wc -l
zcat expression_matrix.tsv.gz | head -n1 | grep -oE 'GTEX-' | wc -l
zcat expression_matrix.tsv.gz | head -n1 | grep -oE 'TARGET-' | wc -l
zcat expression_matrix.tsv.gz | head -n1 | grep -oE '\sK-' | wc -l
```
returned:
| Source        | Sample size   |
| ------------- |:-------------:|
|**TOT**        |19,109         |
|TCGA           |10,530         |
|GTEX           |7,775          |
|TARGET         |734            |
|K              |70             |

#### Gene check
```bash
zcat expression_matrix.tsv.gz | wc -l
```
returned: `60499`.



# Run `heatmaps` workflow on test data set
## Generate the reduced data set for testing
This pipeline makes use of `metasample.py`
```bash
kerblam run gen_test_data
```
> [!NOTE]  
> The accuracy of this step is verified by the `metasample` section of the
> `profiler_tests.R` script.

## Run the analysis pipeline
Edit the `./data/in/config/heatmaps_runtime_options.json` JSON file based on the desired options
```json
{
    "rank_method": "norm_fold_change",
    "threads": 3,
    "save_extra_plots": false,
    "prune_similarity": 0.9,
    "prune_direction": "bottomup",
    "run_unweighted": false,
    "alpha_threshold": 0.20,
    "cluster_heatmap_cols": false
}
```
In particular, use `generanker --list-methods` within the Python virtual environment, to see the available ranking metrics implemented by _Gene Ranker_. Then
```bash
kerblam run heatmaps -l --profile test
```
this will run the following modules (from `./src/modules/`) along with the
related dependencies:
1. `ranking/select_and_run.py`
	- metasplit
	- gene_ranker
		- fast-cohen
1. `make_genesets.py`
	- bonsai
1. `run_gsea.R`
1. `plotting/plot_large_heatmap.R`

### Rank genes
`metasplit` is the program used by `select_and_run.py` to parse the JSON query
file (`./data/in/config/DEA_queries/dea_queries.json`) and extract within-group (i.e., for each cancer type) case and control submatrices from the global expression matrix.
Then, for each cancer type, `gene_ranker` is run to calculate the ranking metric selected through the `rank_method` JSON property.
Final ranks are saved in the `./data/deas/` directory as two-coulmn (sample,ranking) CSV tables.

### Make the gene sets
#### 1. Generate large tables
`make_genesets.py` uses `make_large_tables()` Ariadne's function to make 9
fundamental "large tables" based on the queries hardcoded in
`./data/in/config/gene_lists/basic.json`
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
#### 2. Generate lists from large tables
For each large table, `bonsai` is used to generate tree structures representing
all the possible gene sets, based on the the 3 parameters of the function
`generate_gene_list_trees()`, with the following default values:
```python
min_pop_score: float = 0.5,
min_set_size: int = 10,
min_recurse_set_size: int = 40,
```
with the following meaning:
- `min_pop_score`: minimum portion of non-NA values in a column to be considered for gene lists.
- `min_set_size`: Minimum number of genes to produce a valid gene set.
- `min_recurse_set_size`: minimum parent-geneset size to have before running
recursion on it (effective if `recurse` boolean is `True`).

These parameters cannot be assigned by editing the `heatmaps_runtime_options.json` JSON file because they are considered lower-level parameters, however they can be passed as arguments to the `make_genesets.py` script within the `heatmaps.makefile` workflow. 
```make
# E.g.,
python $(mods)/make_genesets.py ./data/MTPDB.sqlite ./data/in/config/gene_lists/basic.json \
	./data/genesets.json ./data/genesets_repr.txt \
	--min_pop_score 0.7 \
	--min_set_size 15 \
	--min_recurse_set_size 20 \
	--prune_direction $(PRUNE_DIRECTION) \
	--prune_similarity $(PRUNE_SIMILARITY) \
	--verbose
```

#### 3. Make the union of the genesets following the structure
All gene sets are merged together into a large tree structure, then the
`peune()` function is used to remove redundancy. The two parameters of `prune()`
function (`similarity` and `direction`) are set in the
`./data/in/config/heatmaps_runtime_options.json` file, with the following
defaults:
```json
"prune_similarity": 0.5,
"prune_direction": "bottomup",
```




Notes on DESeq2

1. You have to set the option `minReplicatesForReplace` to `Inf` in `DESeq` in
	order to never replace outliers (and so have the baseMeans exactly equal to
	the mean of the MOR-normalized counts across **all** samples)
	`dds2 <- DESeq2::DESeq(dds, minReplicatesForReplace = Inf)`
1. **log2 Fold Change in DESeq2 is not identical to FC calculated from
	normalized count** (https://support.bioconductor.org/p/p134193/)
	_...turning off fold change shrinkage should make log2foldchange from
	DESeq2 be simply equal to (mean of normalized counts group B) / (mean of
	normalized counts group A). However, it seems that some degree of fold
	change moderation is done even when betaPrior is False._
	**It's not always equal to the ratio of the mean of normalized counts
	depending on the fit of the GLM, but close (when no other factors are
	present in the design).** Michael Love
