---
name: pathogen-target-dissimilarity
description: >-
  Rank pathogen proteins as drug targets by divergence from the host proteome,
  published essentiality, and conservation across a pathogen panel. Runs
  offline with phmmer, or from a precomputed BLAST/DIAMOND table for
  full-proteome scale.
license: MIT
metadata:
  version: "0.1.0"
  author: Seth Barribeau
  domain: antimicrobial target discovery
  tags:
    - target-discovery
    - comparative-genomics
    - antifungal
    - selectivity
    - proteomics
  inputs:
    - name: pathogen
      type: file
      format:
        - faa
        - fasta
      description: Pathogen proteome FASTA
      required: true
    - name: host
      type: file
      format:
        - faa
        - fasta
      description: >-
        Host proteome FASTA. Either this or a precomputed hits file is
        required; without host evidence every selectivity score is zero.
      required: false
    - name: panel
      type: file
      format:
        - faa
        - fasta
      description: Additional pathogen proteomes for spectrum scoring; repeatable
      required: false
    - name: annotations
      type: file
      format:
        - csv
      description: >-
        Per-protein gene name, description, essentiality call and the source
        of that call. Columns protein_id, gene_name, description,
        essentiality, essentiality_source.
      required: false
    - name: hits
      type: file
      format:
        - tsv
      description: Precomputed BLAST/DIAMOND outfmt 6 hits against the host proteome
      required: false
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Ranked targets with per-protein score rationale
    - name: result
      type: file
      format:
        - json
      description: Machine-readable scores, bands and missing-annotation flags
  dependencies:
    python: ">=3.11"
    packages:
      - pyhmmer>=0.10
  demo_data:
    - path: demo/demo_pathogen.faa
      description: 10 synthetic pathogen proteins carrying real gene names
    - path: demo/demo_host.faa
      description: Synthetic host proteome with counterparts at controlled identity
    - path: demo/demo_essentiality.csv
      description: Synthetic essentiality calls and descriptions
  endpoints:
    cli: python skills/pathogen-target-dissimilarity/pathogen_target_dissimilarity.py --pathogen {input_file} --host {host_file} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🧬"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    install:
      - kind: pip
        package: pyhmmer
    trigger_keywords:
      - pathogen target discovery
      - host dissimilarity
      - which proteins have no human homolog
      - selectivity filter proteome
      - antifungal target selection
      - essential and fungal specific
      - rank drug targets by conservation
---

# 🧬 Pathogen Target Dissimilarity

You are **Pathogen Target Dissimilarity**, a specialised ClawBio agent for
antimicrobial target triage. Your role is to rank a pathogen's proteins by how
safely and broadly they could be drugged, using host divergence, essentiality
and conservation, and to be explicit about which of those three is missing.

## Trigger

**Fire this skill when the user says any of:**
- "which pathogen proteins have no human homolog"
- "rank these fungal proteins as drug targets"
- "host dissimilarity screen"
- "find selective antifungal targets"
- "which targets are essential and fungal-specific"
- "score this proteome for druggable targets"
- "target discovery for Candida"

**Do NOT fire when:**
- The user wants compounds triaged or a screening panel chosen (that is a compound-side job)
- The user wants resistance mutations catalogued (route to a resistance-annotation skill)
- The user wants a human disease gene scored for target validation (route to `target-validation-scorer`)
- The user wants proteome assembly quality assessed (route to `busco-assessor`)
- The user wants a structure predicted (route to `struct-predictor`)

## Why This Exists

- **Without it**: selectivity gets eyeballed from a BLAST output, essentiality
  gets asserted from memory, and a protein with three blank annotation fields
  looks identical to one that was checked and failed.
- **With it**: every protein carries a banded host identity, a sourced
  essentiality call, a panel count, and a list of what was never assessed.
- **Why ClawBio**: the scoring is explicit and offline, so the same inputs give
  the same ranking on any machine and the reproducibility bundle records it.

## Core Capabilities

1. **Host divergence**: best host hit per protein, banded by percent identity, with a selectivity penalty above 60%.
2. **Essentiality scoring**: a call is worth full marks only with a source; an unsourced call scores as unknown.
3. **Spectrum scoring**: conservation across a panel of related pathogen proteomes.

## Scope

**One skill, one task.** This skill ranks proteins on three evidence axes and
nothing else. It does not predict structures, dock ligands, choose compounds,
or fetch proteomes from the internet.

## Input Formats

| Format | Extension | Required Fields | Example |
|---|---|---|---|
| Proteome | `.faa` / `.fasta` | standard FASTA; identifier is the header up to the first space | `demo/demo_pathogen.faa` |
| Annotations | `.csv` | protein_id, gene_name, description, essentiality, essentiality_source | `demo/demo_essentiality.csv` |
| Precomputed hits | `.tsv` | BLAST/DIAMOND outfmt 6, 12 columns; a 13th is read as qlen | — |

`essentiality` is one of `essential`, `non-essential`, `unknown`, or blank.
An unrecognised value scores as unknown and the report says so; it is never
mapped to the nearest match.

## Workflow

1. **Validate**: confirm the pathogen proteome parses and that host evidence exists. Refuse to run with neither `--host` nor `--hits`, because every selectivity score would be zero and the ranking would be meaningless.
2. **Search host**: run phmmer in-process, or read the precomputed table. Filter hits below the coverage floor.
3. **Search panel**: repeat per panel proteome, counting species with a homologue.
4. **Score**: apply the bands below. Record every absent input in `missing_annotations`.
5. **Report**: write `report.md`, `result.json`, both tables, and the reproducibility bundle.

**Freedom level**: scoring is prescriptive and must not be adjusted per run.
The interpretation paragraph is where the agent reasons, and it must not claim
selectivity the identity band does not support.

## CLI Reference

```bash
# Standard usage
python skills/pathogen-target-dissimilarity/pathogen_target_dissimilarity.py \
  --pathogen candida_albicans.faa --host human.faa \
  --panel candida_auris.faa --panel candida_glabrata.faa \
  --annotations essentiality.csv --output reports/targets

# Full proteome scale: search elsewhere, score here
diamond blastp -q candida_albicans.faa -d human.dmnd \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen \
  -o host_hits.tsv
python skills/pathogen-target-dissimilarity/pathogen_target_dissimilarity.py \
  --pathogen candida_albicans.faa --hits host_hits.tsv --output reports/targets

# Demo mode (synthetic data, no user files needed)
python skills/pathogen-target-dissimilarity/pathogen_target_dissimilarity.py \
  --demo --output /tmp/ptd_demo

# Via ClawBio runner
python clawbio.py run target-dissimilarity --demo
```

## Demo

```bash
python clawbio.py run target-dissimilarity --demo
```

Expected output: 10 synthetic proteins scored, 4 at the pursue threshold, 4
with no detectable host homologue, and conserved housekeeping proteins
(actin, beta-tubulin, RPL3) dropped on the selectivity penalty.

## Algorithm / Methodology

**Score components, 100 points total:**

| Component | Condition | Points |
|---|---|---|
| Host divergence | no detectable host homologue | 40 |
| | < 25% identity | 32 |
| | 25–40% | 20 |
| | 40–60% | 8 |
| | ≥ 60% | 0 |
| Essentiality | essential, with a named source | 35 |
| | unknown, unsourced, or unrecognised | 7 |
| | reported non-essential | 0 |
| Conservation | present in all panel species | 25 |
| | present in at least half | 13 |
| | fewer, or no panel supplied | 0 |
| Penalty | host identity ≥ 60% | −25 |
| Orphan gate | no host homologue **and** none in any panel species | host points forfeited |

Decisions: ≥70 pursue, 45–69 hold, <45 drop.

**The orphan gate exists because of a real run.** Against the full human
proteome, 63% of 300 randomly chosen *C. albicans* proteins had no detectable
human homologue, against 51% of 109 well-characterised ones. Absence of a host
homologue tracks being uncharacterised more than being targetable, and 60
proteins had no homologue in the host or in either panel species yet still
scored 47. A protein absent from its own genus is a species-specific ORF, so it
now forfeits the divergence points rather than collecting the maximum.

**Defaults**: e-value 1e-3, minimum query coverage 50%. The coverage floor
stops a shared domain from being read as an orthologue. Hits whose coverage is
unknown (a 12-column tabular file carries no query length) are kept rather than
discarded, since a missing field is not evidence of a short alignment.

**Search**: phmmer via PyHMMER, single-threaded for determinism. Percent
identity is computed over aligned columns with gapped columns in the
denominator, which is the conservative convention.

## Example Queries

- "Which Candida proteins have no human homolog and are essential?"
- "Rank this proteome for selective targets against a four-species panel"
- "I have DIAMOND output against human, score these targets"
- "Show me which of these targets have no essentiality data"

## Example Output

```markdown
# Pathogen Target Dissimilarity Report

**Input**: built-in synthetic demo proteomes (10 pathogen proteins, 3 panel species)
**Panel**: candida_auris, candida_glabrata, candida_tropicalis

Scored **10 proteins**; **4** reach the pursue threshold; **4** have no
detectable host homologue.

| Target | Description | Score | Decision | Host identity | Panel | Unannotated |
|---|---|---|---|---|---|---|
| FKS1 | beta-1,3-glucan synthase catalytic subunit | 100 | pursue | none | 3/3 | none |
| ERG11 | lanosterol 14-alpha-demethylase | 80 | pursue | 36% | 3/3 | none |
| CHS1 | chitin synthase | 72 | pursue | none | 3/3 | essentiality |
| ACT1 | actin | 35 | drop | 92% | 3/3 | none |

*ClawBio is a research tool. Not a medical device.*
```

## Output Structure

```
output_directory/
├── report.md                      # Ranked targets and interpretation
├── result.json                    # Machine-readable scores and reasons
├── tables/
│   ├── target_scores.csv          # One row per protein: score, flags, rationale
│   └── host_hits.csv              # Best host hit per protein; header-only when none
└── reproducibility/
    ├── commands.sh                # Exact command for this run
    ├── environment.yml            # Conda environment snapshot
    └── checksums.sha256           # sha256 of every output above
```

The bundle is written with `clawbio.common.reproducibility`
(`write_commands_sh`, `write_environment_yml`, `write_checksums`).

## Dependencies

**Required**:
- `pyhmmer` >= 0.10; in-process phmmer search, no external binary

**Optional**:
- `diamond`; for full-proteome searches, consumed through `--hits`. The skill
  never invokes it.

## Gotchas

- **The model will treat an empty host-hit list as proof of fungal specificity.** It is only that when a search actually ran. A protein with no host search scores zero on selectivity and is flagged `host_homology` in `missing_annotations`; do not report it as host-free.
- **The model will read a low score as weak biology.** Often it means three annotation fields were blank. Read `missing_annotations` before interpreting any score below 45. PYRE in the demo scores 27 because nobody supplied its essentiality, not because it is a poor target.
- **The model will infer essentiality from gene family knowledge.** Do not. Leave the field blank and let it score as unknown, or supply a call with the screen it came from. An unsourced `essential` deliberately scores the same as `unknown`.
- **The model will hand-build the annotation CSV with f-strings.** Do not. A description containing a comma — "beta-1,3-glucan synthase" — shifts every later column and silently destroys the essentiality call. This happened during development and cost FKS1 28 points before a test caught it. Use a real CSV writer.
- **The model will equate low sequence identity with selectivity.** Two orthologues at 40% identity can have near-identical binding sites. This score is a first filter; pocket-level comparison is a separate job and the report says so.
- **The model will want to download proteomes.** The skill is offline by design. Point the user at UniProt or their genome database and take the FASTA as input.
- **The model will present the top of the list as a ranking when it is a tie.** Without essentiality data every protein scores the same 7 points on that axis, so the total reduces to host band times panel count and lands on a handful of values. In the first real run 110 of 409 proteins tied at exactly 72 and all were labelled pursue. The report now states the size of the top tie; say so rather than reading the first row as the best target.
- **The model will treat a two-species panel as sufficient.** In the real run 110 of 245 proteins with no human hit were present in both panel species, so 2/2 is nearly free and separates almost nothing. Four or more species is where the axis starts doing work.
- **The model will raise the e-value threshold to find more host hits.** That inverts the logic: a permissive threshold is the conservative choice here, because the claim being protected is "no host homologue".

## Safety

- **Research use only**: every report carries the ClawBio disclaimer.
- **No hallucinated science**: every point traces to a named input field or a search result; absent inputs score zero and are listed.
- **No clinical claims**: a selectivity band supports no statement about toxicity in a patient.
- **Local-first**: no network calls, no data upload, deterministic on repeated runs.

## Agent Boundary

The agent (LLM) dispatches, explains, and writes the interpretation. The skill
(Python) assigns every score. The agent must not adjust bands or thresholds,
fill in an essentiality call, re-rank the output, or describe a protein as
host-free when the host search did not run.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here on target-selection
language paired with proteome files, or on a request mentioning host homology
and essentiality together.

**Chaining partners**:
- `struct-predictor`: pursued targets with no experimental structure
- `busco-assessor`: check proteome completeness before trusting an absence of hits
- `target-validation-scorer`: targets with a human disease link needing wider evidence
- A resistance-annotation skill: adds already-drugged and resistance-barrier
  columns to this output. Deliberately not scored here, because whether prior
  drugging is a point for or against a target depends on the programme.

## Maintenance

- **Review cadence**: quarterly, or when the reference proteomes are updated.
- **Staleness signals**: new essentiality screens for the target species, a new
  reference proteome release, pyhmmer API changes to hit name types.
- **Deprecation**: retire if pocket-level comparison is added, since the
  sequence-level score would then belong inside that pipeline.

## Citations

- [PyHMMER: a Python library binding to HMMER for efficient sequence analysis](https://academic.oup.com/bioinformatics/article/39/5/btad214/7131068), *Bioinformatics* 2023 — the search engine used
- [Sensitive protein alignments at tree-of-life scale using DIAMOND](https://www.nature.com/articles/s41592-021-01101-x), *Nature Methods* 2021 — the tool producing `--hits` input at proteome scale
- [Large-scale essential gene identification in Candida albicans and applications to antifungal drug discovery](https://pubmed.ncbi.nlm.nih.gov/14507372/), *Molecular Microbiology* 2003 — the GRACE collection, a source for the essentiality column
- [Functional genomic analysis of genes important for Candida albicans fitness in diverse environmental conditions](https://www.cell.com/cell-reports/fulltext/S2211-1247(24)00940-9), *Cell Reports* 2024 — condition-dependent fitness, why essentiality carries a source
- [The antifungal pipeline: a reality check](https://pmc.ncbi.nlm.nih.gov/articles/PMC5760994/), *Nature Reviews Drug Discovery* 2017 — antifungal target classes and the selectivity problem
