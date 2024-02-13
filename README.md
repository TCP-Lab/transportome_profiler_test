# transportome_profiler_test
R tests for the `transportome_profiler` project


# How to locally run the ___Transportome Profiler___

## Clone the repo
```bash
git clone git@github.com:TCP-Lab/transportome_profiler.git
cd ./transportome_profiler
```

## Install R requirements
```bash
sudo Rscript ./src/helper_scripts/install_r_pkgs.R
```

## (Create and) Activate a Python virtual environment and install requirements
```bash
python -m venv env
source ./env/bin/activate
pip install -r ./src/requirements.txt

# When finished:
deactivate
```

## Install `xsv`
This is a program required by `metasplit` for the fast reshaping of very large
CSV files
```bash
sudo pacman -Syu xsv
```

## Fetch MTP-DB, TCGA and GTEx datasets from remote (through XENA)
```bash
kerblam data fetch
```
will download in `/data/in` the following three archives
1. `expression_matrix.tsv.gz`
1. `expression_matrix_metadata.tsv.gz`
1. `MTPDB.sqlite.gz`

In particular, `expression_matrix.tsv.gz` is the
`gene expression RNAseq - RSEM expected_count` data set from
[XENA](https://xenabrowser.net/datapages/?dataset=TcgaTargetGtex_gene_expected_count&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443);
- author: _UCSC TOIL RNA-seq recompute_
- unit: __log2(expected_count+1)__
- size: 60,499 identifiers (genes) x 19,109 samples
- download: https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_gene_expected_count.gz

### Sample check
```bash
head -n1 expression_matrix.tsv | grep -E -o 'TCGA-|GTEX-|TARGET-|\sK-' | wc -l
```
- TOT 		19109
- TCGA		10530
- GTEX 		7775
- TARGET 	734
- K 		70

### Gene check
```bash
wc -l expression_matrix.tsv
```
60499

## Generate the reduced data set for testing purposes by `metasample.py`
```bash
kerblam run gen_test_data
```

## Run the analysis pipeline on test dataset
```bash
kerblam run heatmaps -l --profile test
```
this will run
1. `select_and_run.py`
	- metasplit
	- generanker
		- fast-cohen




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
