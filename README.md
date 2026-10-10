# ISpy

Snakemake pipeline that finds insertion sequence (IS) elements in Illumina WGS data, keeps a curated local IS database up to date, and shows which genes each insertion hits.

Read preprocessing is the same as [ISScreener](https://github.com/kylegontjes/ISScreener) (Kyle Gontjes). The ISfinder web search is replaced by a local BLAST against ISFinder, ISOSDB or both.

## Workflows

Run them in this order. Each one reads the previous one's results.

| Snakefile | What it does | Main outputs |
|---|---|---|
| `ISpy.smk` | Trimmomatic → bwa → Picard → panISa. Then local BLAST of the IS ends against `blast_dbs`, genes near each site (bedtools), plots, and ISEScan on each assembly | `results/<prefix>/<sample>/{blast,annotation,plots}/`, `results/<prefix>/<prefix>_<db>_master.tsv` |
| `update_db.smk` | Merges the `base` databases with every ISEScan element and removes duplicates. Elements that match nothing are renamed `IS_novel_NNNN_<family>`. The result is saved as a new database | `IS_db/<name>/<name>.fna` and reports |
| `reblast.smk` | Re-BLASTs every panISa site against the new database, then re-annotates and re-plots | same as `ISpy.smk` for `<name>`, plus `<prefix>_<name>_no_identity.tsv` |

## Setup

1. Edit `config/config.yaml` and `config/samples.tsv` (one `sample_id` column).
2. Put the input data in place (not tracked in git):
   - reads: `fastq/<sample>_R1.fastq.gz`, `fastq/<sample>_R2.fastq.gz`
   - assemblies: `assembly/<sample>/<sample>.fna`
   - reference FASTA and GFF: `ref/`
   - IS databases: `IS_db/ISFinder/IS.fna` (from [thanhleviet/ISfinder-sequences](https://github.com/thanhleviet/ISfinder-sequences)), `IS_db/ISOSDB/ISOSDB.V3.fna`, and the ISOSDB family table `IS_db/ISOSDB/IS_fam_annot.txt` (columns `IS`, `IS_fam`)
3. Install Snakemake (≥ 7.32) and conda/mamba. Every rule's tools come from `envs/*.yaml`.

## Run

```bash
snakemake -s ISpy.smk     --use-conda --conda-prefix ~/.ispy_conda -j 8
snakemake -s update_db.smk --use-conda --conda-prefix ~/.ispy_conda -j 8
snakemake -s reblast.smk   --use-conda --conda-prefix ~/.ispy_conda -j 8
```

Add `-n -p` for a dry run. On Great Lakes, `sbatch ISpy.sbat` runs all three (edit the email and account first).

Use `--conda-prefix` with a path that has no spaces. Conda tools break when their environment path contains spaces, such as a Dropbox folder.

## Key options (`config/config.yaml`)

- `blast_dbs`: which databases `ISpy.smk` searches: `ISFinder`, `ISOSDB` and/or `combined`. `combined: auto` builds ISFinder + ISOSDB with exact duplicates removed.
- `annotation.method`: how genes are picked for each site.
  - `window`: every gene within `distance` bp (default 1000).
  - `closest`: the `n_closest` nearest genes (default 3).
- `update_db`: which `base` databases to merge, the new database's `name` and `outdir`, the identity/coverage cut-offs for calling an element known, and `panisa_rescue`. The rescue cuts elements behind "No identity" sites out of the assemblies.
- `trim: false` with `trimmed_reads`: start from reads that are already trimmed.

## Scripts

| Script | Purpose |
|---|---|
| `scripts/local_db_blast.py` | Local version of panISa's `ISFinder_search.py`; same output table |
| `scripts/update_IS_db.py` | Merges databases, clusters ISEScan elements (cd-hit-est), names novel IS, runs the panISa rescue |
| `scripts/gff_to_bed.py`, `make_sites_bed.py`, `summarize_sites.py` | Inputs for bedtools and the per-site gene table |
| `scripts/plot_IS_sites.R`, `compare_dbs.R` | Gene neighbourhood, genome overview and compressed plots; ISFinder vs ISOSDB Venn |
| `scripts/bwa.sh` | bwa mem + read group + sort (from ISScreener) |

## Credits

- Preprocessing and `bwa.sh` are adapted from [ISScreener](https://github.com/kylegontjes/ISScreener) © 2024 Kyle Gontjes (MIT License).
- `local_db_blast.py` reimplements `ISFinder_search.py` from [panISa](https://github.com/bvalot/panISa) (Treepong et al., GPL-3.0).
