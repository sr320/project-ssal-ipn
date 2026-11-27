# project-ssal-ipn

DNA methylation response of Atlantic salmon (*Salmo salar*) to infectious pancreatic necrosis virus (IPNV), by whole-genome bisulfite sequencing (WGBS).

**Dashboard:** https://sr320.github.io/project-ssal-ipn/ (progress and major issues; served from `docs/`)

## Data

- Macrogen order HN00292532, delivered 2026-11-26. Report in `20261126_HN00292532_WGB_Report/`.
- 18 libraries: groups G1, G2, G6 with six libraries each (HK1–HK6, head kidney). Design in `data/samples.tsv` and below.
- Paired-end 2 × 151 bp, Accel-NGS Methyl-Seq library kit, EZ DNA Methylation-Gold conversion.
- 36 FASTQ files, 85 GB compressed, 218 Gb of sequence. 34.7–49.2 M read pairs per library, about 3.5–5.0× raw depth of the ~3 Gb genome.
- Download links and md5sums: `20261126_HN00292532_WGB_Report/assets/spgs/HN00292532_18samples_md5sum_DownloadLink.txt`. **Links expire two weeks after delivery (about 2026-12-10).**

## Experimental design

Fish were acclimated in fresh water (FW), moved to seawater (SW) at 40 days post infection (dpi, counted from the FW challenge on day 1), and all sampled for head kidney 14 days after the SW challenge (day 68).

| Group | Label | IPNV in FW (day 1) | IPNV in SW (54 dpi) | What it represents |
|---|---|---|---|---|
| G1 | CTRL | no | no | Uninfected control |
| G2 | SW-I | no | yes | First infection, in seawater |
| G6 | FW+SW-I | yes | yes | Infected in fresh water, reinfected in seawater |

Contrasts for differential methylation:

- **G2 vs G1**: response to a primary IPNV infection in SW, 14 days in.
- **G6 vs G1**: response to reinfection in fish with a prior FW infection.
- **G6 vs G2**: effect of the earlier FW infection on the response to the same SW challenge (priming or lasting imprint). This is the main question the design adds.

Still to confirm: one fish per library, fish IDs, tank per group (if each group was one tank, tank and group are confounded), whether G1 got a mock challenge, IPNV strain and dose, and whether G6 fish carried virus before the SW challenge. G3–G5 are presumably other trial groups that were not sequenced.

## Workflow

Gannet holds the data; Hyak (klone) only computes.

- Raw reads are downloaded from Macrogen straight to gannet. Hyak never holds the full set.
- Each Hyak job pulls the one library it needs from gannet over https, works on `/gscratch/scrubbed`, and deletes the reads when it finishes.
- Finished files are pushed back to gannet with rsync. Files no longer needed on Hyak are removed there once the transfer succeeds.
- Jobs go to the checkpoint partition (`ckpt-all`, account `srlab`) first. The one exception is the genome index build.
- The dashboard links to the MultiQC reports and per-library files on gannet.

| # | Step | Script | Where | Ends up on gannet |
|---|---|---|---|---|
| 0 | Check container, gannet addresses, ssh | `code/00-preflight.sh` | klone login node | |
| 1 | Download from Macrogen, md5 check | `code/01-download-to-gannet.sh` | **gannet** | `raw/` (36 FASTQ) |
| 2 | Raw read QC, flowcell/lane per library | `code/02-fastqc-raw.sh` | ckpt array (18) | `output/qc/fastqc-raw/` |
| 3 | Ssal_v3.1 genome and Bismark index | `code/03-genome-index.sh` | lab partition | (kept in `/gscratch/srlab`) |
| 4 | Trim, align in chunks, deduplicate, M-bias | `code/04-trim-align-dedup.sh` | ckpt array (18) | logs, reports, M-bias |
| 5 | Methylation extraction, merged CpG coverage | `code/05-extract.sh` | ckpt array (18) | `output/methylation/`, `output/dedup/*.bam` |
| 2b | MultiQC across libraries | `code/02b-multiqc.sh` | ckpt job | `output/qc/multiqc-*/` |
| 6 | Numbers and links for the dashboard | `code/06-collect-metrics.sh` | klone login node | |
| 7 | Catch-up push of anything jobs could not send | `code/07-archive-to-gannet.sh` | klone login node | |
| 8 | Sample QC: coverage, correlation, PCA | to write | Hyak job (R) | |
| 9 | Differential methylation (DSS): G2 vs G1, G6 vs G1, G6 vs G2 | to write | Hyak job (R) | |
| 10 | Annotation and enrichment | to write | Hyak job / local | |

### Set up once

1. The gannet project directory is `/volume1/v1_web/project-ssal-ipn`, read by Hyak at `https://gannet.fish.washington.edu/v1_web/project-ssal-ipn` (`GANNET_SSH` and `GANNET_WEB` at the top of `code/00-config.sh`). `GANNET_SSH` logs in with your klone username; change it there if your gannet username differs.
2. On gannet (or a machine with it mounted), inside tmux: `bash code/01-download-to-gannet.sh <that directory>/raw`.
3. For jobs to push results themselves, klone needs a passwordless ssh key for gannet. Without one, jobs leave results on scrubbed and `code/07-archive-to-gannet.sh` (which asks for the password) moves them.

### Running on Hyak

```bash
ssh klone.hyak.uw.edu
cd /gscratch/scrubbed/$USER
git clone https://github.com/sr320/project-ssal-ipn.git && cd project-ssal-ipn
mkdir -p logs
bash code/00-install-cutadapt.sh       # once, on a compute node (salloc); the container lacks cutadapt
bash code/00-preflight.sh              # everything must say ok
sbatch code/03-genome-index.sh         # can start before the reads are staged
sbatch code/02-fastqc-raw.sh
sbatch code/02b-multiqc.sh             # MultiQC, pushed to gannet; review before trimming
sbatch --array=1 --export=ALL,SUBSET=200000,MAXINS=500  code/04-trim-align-dedup.sh   # pilot
sbatch --array=1 --export=ALL,SUBSET=200000,MAXINS=1000 code/04-trim-align-dedup.sh   # pilot
sbatch code/04-trim-align-dedup.sh     # full run (defaults CLIP=15, MAXINS=1000)
sbatch code/05-extract.sh              # ignores 3 bp of R1 and 9 bp of R2 (from the pilot's M-bias)
bash code/07-archive-to-gannet.sh      # only if jobs could not push
bash code/06-collect-metrics.sh        # then commit and push docs/ and output/
```

### Pilot: one library first

Depth is low (about 3.5-5.0x raw per library), so one library goes through the whole pipeline before the rest, to measure what survives. It also tests the scripts. The pilot fetches `G1_HK1` straight from Macrogen (`RAW_SOURCE=macrogen`), so it needs neither gannet nor the other 34 files; the genome index must be built first.

```bash
sbatch code/03-genome-index.sh
# after the index is built:
sbatch --array=1 --export=ALL,RAW_SOURCE=macrogen,SUBSET=200000,MAXINS=500  code/04-trim-align-dedup.sh
sbatch --array=1 --export=ALL,RAW_SOURCE=macrogen,SUBSET=200000,MAXINS=1000 code/04-trim-align-dedup.sh
grep "Mapping efficiency" logs/04-trim-align_*    # pick MAXINS
sbatch --array=1 --export=ALL,RAW_SOURCE=macrogen code/04-trim-align-dedup.sh
sbatch --array=1 code/05-extract.sh               # after 04 finishes
cat /gscratch/scrubbed/$USER/project-ssal-ipn/methylation/G1_HK1.coverage.tsv
bash code/06-collect-metrics.sh                   # then commit and push
```

`G1_HK1.coverage.tsv` holds mean CpG depth (strands merged) and the share of CpGs with at least 1, 3, 5 and 10 reads; `G1_HK1.coverage-hist.tsv` has the full distribution. With the duplication rate from `dedup/G1_HK1.deduplication_report.txt`, these decide between analysing regions at the current depth and asking for deeper sequencing. `G1_HK1` has median yield; the lowest-yield libraries (`G1_HK4`, `G2_HK2`, `G6_HK4`) will be about 12% lower.

If gannet is not set up yet, the pushes at the end of each job fail harmlessly and the files stay on scrubbed.

### Checkpoint partition

Hyak stops and requeues checkpoint jobs every 4–5 hours and can preempt them at any time. The scripts are written for that:

- Every stage leaves a marker file, so a requeued task skips what is done.
- Alignment runs in chunks of 10 million pairs; a preempted task repeats at most one chunk. Duplicates are then removed across all chunks of a library together.
- The genome index build cannot resume, so it is sent to `cpu-g2-mem2x`.
- Methylation extraction cannot resume part-way either. It should fit in the window at this depth; this is untested.
- If ckpt is not enough (a library keeps restarting, or the queue is slow), resubmit on an allocation that is not preempted. Command-line options override the `#SBATCH` lines:

| Account / partition | Limit | Example |
|---|---|---|
| `srlab` / `cpu-g2-mem2x` | all of it | `sbatch --partition=cpu-g2-mem2x --array=7 code/04-trim-align-dedup.sh` |
| `coenv` / `compute` | | `sbatch --account=coenv --partition=compute --array=7 code/04-trim-align-dedup.sh` |
| `coenv` / `cpu-g2` | 28 CPUs | `sbatch --account=coenv --partition=cpu-g2 --cpus-per-task=28 --array=7 code/04-trim-align-dedup.sh` |

  On `coenv`/`cpu-g2`, 28 CPUs is one alignment task at a time (the script asks for 32 by default; Bismark's `PARALLEL=6` fits 28).

### What is kept where

| Files | Hyak scrubbed | Gannet |
|---|---|---|
| Raw FASTQ | only while a job uses that library | kept |
| Trimmed FASTQ | deleted after alignment | not kept (`KEEP_TRIMMED=1` moves them instead) |
| Per-chunk BAMs | deleted after deduplication | not kept |
| Deduplicated BAM | until extraction finishes | moved |
| CpG report, bedGraph, coverage, per-read context files | until pushed | moved |
| Merged CpG coverage (`*.merged_CpG_evidence.cov.gz`) | kept for differential methylation | copied |
| FastQC, MultiQC, logs, Bismark reports, M-bias | kept (small) | copied |
| Genome index | `/gscratch/srlab/$USER/project-ssal-ipn/genome` (group storage, not scrubbed) | not copied |

Scrubbed deletes files 30 days after creation, so the merged coverage files may need pulling back from gannet for later steps. The trim/align and extraction arrays run six libraries at a time (`%6`) to limit how much sits on scrubbed at once; my rough estimate is 30-40 GB per library in flight.

### Decision points

1. **After raw FastQC:** length of the adaptase tail to clip (`CLIP`, default 15 bp); whether flowcell/lane lines up with group.
2. **After the alignment pilot:** maximum fragment length (`MAXINS` 500 or 1000); expected mapping rate and usable depth.
3. **After M-bias:** bases to ignore at read ends in extraction.
4. **After sample QC:** coverage filter and any libraries to drop.
5. **Before differential methylation:** group definitions and contrasts.

### Status of the scripts

The trimming and Bismark settings come from the tested scripts in [project-salmon-vaccine-wgbs](https://github.com/sr320/project-salmon-vaccine-wgbs) (same kit, same genome). The SLURM wrapping, gannet streaming and chunked alignment are new. They pass `bash -n`, and the gannet push and metrics logic were tested against mock files, but **nothing has been run on Hyak**. To confirm on the first run:

- compute nodes can reach gannet and NCBI over https, and gannet over ssh;
- the output file names of `deduplicate_bismark --multiple` and the extractor match what `05-extract.sh` and the archive step expect (run one library through first: `--array=1`);
- the cores and memory requested suit ckpt.

## Updating the dashboard

- Steps and issues: edit `docs/data.js`.
- Per-library numbers: `bash code/06-collect-metrics.sh` on klone rewrites `docs/metrics.js`, including the gannet links.
- Commit and push; GitHub Pages serves `docs/` from `main` (Settings → Pages → Deploy from a branch → `main` / `docs`).
