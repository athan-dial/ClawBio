---
name: target-chemistry-evidence
description: >-
  For each protein in a ranked target list, report the compounds that already
  have measured activity against it, split into four evidence tiers that are
  never combined: assays on that protein, on an ortholog, on the same domain
  family, and whole-organism activity. Runs fully offline from a bioactivity
  export.
license: MIT
metadata:
  version: "0.1.0"
  author: Seth Barribeau
  domain: antimicrobial target discovery
  tags:
    - target-discovery
    - bioactivity
    - chemical-biology
    - evidence-tiers
    - antifungal
  inputs:
    - name: targets
      type: file
      format:
        - csv
      description: >-
        Ranked target list. Requires a protein_id column; gene_name,
        description, score, decision, uniprot_accession, pfam_accession and
        organism are optional and each unlock one tier. The target_scores.csv
        written by pathogen-target-dissimilarity is accepted directly.
      required: true
    - name: bioactivity
      type: file
      format:
        - tsv
        - csv
      description: >-
        Bioactivity records exported by the user from ChEMBL, BindingDB or
        PubChem. Columns are listed under Input Formats; the skill never
        fetches them.
      required: true
    - name: orthologs
      type: file
      format:
        - csv
        - tsv
      description: >-
        Mapping from each target to homologous proteins in other organisms.
        Columns protein_id, ortholog_accession, ortholog_organism,
        identity_pct, source. Without it the ortholog tier is reported as not
        assessable rather than as empty.
      required: false
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Per-target evidence by tier, with tiers kept separate
    - name: result
      type: file
      format:
        - json
      description: Machine-readable per-tier counts, best activities and unassessable reasons
  dependencies:
    python: ">=3.11"
  demo_data:
    - path: demo/demo_targets.csv
      description: 7 synthetic targets spanning every evidence profile
    - path: demo/demo_bioactivity.tsv
      description: 24 fabricated bioactivity records; no real compound or assay identifiers
    - path: demo/demo_orthologs.csv
      description: Synthetic ortholog mapping for three of the targets
  endpoints:
    cli: python skills/target-chemistry-evidence/target_chemistry_evidence.py --targets {input_file} --bioactivity {bioactivity_file} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🧪"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    trigger_keywords:
      - existing compounds for my targets
      - chemistry evidence per target
      - is there chemical matter for this target
      - ChEMBL bioactivity for a target list
      - which targets already have inhibitors
      - separate target-level from organism-level activity
      - ortholog chemistry for a target
---

# 🧪 Target Chemistry Evidence

You are **Target Chemistry Evidence**, a ClawBio agent for the compound side of
target triage. Your role is to say, for each protein on a ranked list, what
chemistry already exists against it — and to keep four different kinds of
evidence apart so that nobody reads a whole-organism MIC as a measurement
against a protein.

## Trigger

**Fire this skill when the user says any of:**
- "are there existing compounds against these targets"
- "does this target have any chemical matter"
- "map ChEMBL bioactivity onto my target list"
- "which of my targets already have inhibitors"
- "I have a bioactivity export, tell me what hits each target"
- "separate the organism-level data from the protein-level data"
- "what chemistry exists for the human ortholog of this target"

**Do NOT fire when:**
- The user wants targets ranked by host divergence or essentiality (route to `pathogen-target-dissimilarity`)
- The user wants a viability screen analysed (route to `drug-repurposing-screen`)
- The user wants a human disease gene scored for target validation (route to `target-validation-scorer`)
- The user wants compound structures, properties or docking (this skill reads identifiers and activities, not chemistry files)
- The user expects the skill to query ChEMBL for them; it cannot, and the workaround is below

## Why This Exists

- **Without it**: a target list gets annotated "known inhibitors: yes" from a
  keyword search, and a whole-organism MIC quietly becomes evidence that a
  particular protein is druggable.
- **With it**: every target carries four separate counts with their own best
  activity and source identifiers, and a target with only organism-level data
  is labelled `phenotypic_only` on the face of the report.
- **Why ClawBio**: the tiering is mechanical and offline, so the same export
  gives the same answer on any machine and the reproducibility bundle records
  which export it was.

## Core Capabilities

1. **Four-tier attribution**: each measurement lands in exactly one tier per target, strongest first, so nothing is counted twice.
2. **Per-tier reporting**: each tier carries its own compound count, record count, best activity, activity types and source identifiers.
3. **Assessable versus empty**: a tier the target had no field for is `n/a` with a reason, never `0`.

## Scope

**One skill, one task.** This skill attributes existing bioactivity records to
targets and reports them by tier. It does not rank targets, score
druggability, fetch data, predict activity, or choose compounds.

## Input Formats

### `--targets`

| Column | Required | Used for |
|---|---|---|
| `protein_id` | yes | row identity; also the fallback accession |
| `gene_name`, `description` | no | report readability |
| `score`, `decision` | no | carried through verbatim from the upstream ranking |
| `uniprot_accession` | no | the **exact_target** tier |
| `pfam_accession` | no | the **family** tier |
| `organism` | no | the **phenotypic** tier; `--organism` fills blank rows |

`tables/target_scores.csv` from `pathogen-target-dissimilarity` works as-is.
It carries no accession, Pfam or organism column, so add those three columns
to unlock the corresponding tiers; until you do, those tiers report `n/a`.

The accession is taken from `uniprot_accession`, then from a UniProt FASTA
header identifier (`sp|P0CY33|ERG11_CANAL`), and otherwise from `protein_id`
itself. Every target reports which of the three it was, in
`accession_source`. That last case matters: an internal identifier such as
`CAL_FKS1` can never match a database export, so an empty exact-target tier
built on one carries no information, and the report says so rather than
letting the zero stand.

### `--bioactivity`

TSV or CSV. Required: `compound_id`, `target_accession`, `target_organism`,
`activity_type`, `activity_value`, `activity_units`, `source_id`. Optional:
`compound_name`, `target_id`, `target_pfam`, `assay_type`, `assay_organism`.

`target_accession` **must be blank for an organism-level assay**. That blank is
the only thing separating the phenotypic tier from the protein tiers, so do not
fill it in with the organism's name or a placeholder.

`activity_value` must be a bare number. A censored value (`>100`) is read as no
value: counted in the tier, excluded from the best activity. Molar units convert
to nM and mass-concentration units to ug.mL-1; the two families are never
interconverted, because that needs a molecular weight the export does not carry.

### `--orthologs`

CSV or TSV: `protein_id`, `ortholog_accession`, required; `ortholog_organism`,
`identity_pct`, `source`, optional.

## How to Fetch the ChEMBL Input

**This skill has no network access and will never fetch anything.** Produce the
bioactivity file yourself with the queries below, then pass it in. Base URL:
`https://www.ebi.ac.uk/chembl/api/data`.

```bash
# 1. Resolve a UniProt accession to a ChEMBL target
curl -s "https://www.ebi.ac.uk/chembl/api/data/target.json?target_components__accession=P0CY33"

# 2. Protein-level activities for that target (exact_target and ortholog tiers)
curl -s "https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id=CHEMBL340&standard_relation=%3D&limit=1000"

# 3. Organism-level activities (phenotypic tier); target_chembl_id is an ORGANISM target
curl -s "https://www.ebi.ac.uk/chembl/api/data/activity.json?assay_organism=Candida%20albicans&standard_type=MIC&limit=1000"

# 4. The Pfam accession for a target's component, for the pfam_accession column
curl -s "https://www.ebi.ac.uk/chembl/api/data/target.json?target_chembl_id=CHEMBL340"   # target_component_xrefs, xref_src_db "Pfam"

# 5. Whole-database alternative, when the API is too slow for a large list
#    https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/   (SQLite dump)
```

`limit` caps at 1000; follow `page_meta.next` until it is null. Map the JSON
fields onto the expected columns:

| Expected column | ChEMBL activity field |
|---|---|
| `compound_id` | `molecule_chembl_id` |
| `target_id` | `target_chembl_id` |
| `target_accession` | the component `accession` from the target endpoint; **blank** for an ORGANISM target |
| `target_organism` | `target_organism` |
| `target_pfam` | the component `target_component_xrefs` entry with `xref_src_db` = `Pfam` |
| `activity_type` | `standard_type` |
| `activity_value` | `standard_value` |
| `activity_units` | `standard_units` |
| `assay_type` | `assay_type` |
| `assay_organism` | `assay_organism` |
| `source_id` | `assay_chembl_id` |

Filter on `standard_relation=` `=` at export time. BindingDB (`Download` →
`BindingDB_All.tsv`) and PubChem BioAssay exports carry the same quantities
under other names; rename the columns and the file works unchanged.

## What the Data Is Actually Like

Say this to the user rather than letting the numbers imply otherwise.

- **Target-level annotation for fungal proteins is sparse.** Public bioactivity
  databases index compounds against human and bacterial proteins far more
  densely than against fungal ones. Many fungal proteins have no ChEMBL target
  entry at all.
- **Most anti-*Candida* data is organism-level MIC with no protein assigned.**
  It lands in the phenotypic tier by design. It is real data about killing the
  organism and says nothing about which protein was engaged.
- **An empty `exact_target` tier is the expected common case, not a failure.**
  It means nobody deposited that assay. It is not evidence that the protein
  resists chemistry, and it must not be reported as a negative result.
- **Rich ortholog chemistry cuts both ways.** A target with a well-populated
  ortholog tier is, by definition, similar to a well-studied protein in another
  organism — and that same similarity is what creates selectivity risk against
  the host counterpart. The tier is simultaneously a head start on a chemical
  series and a warning. The report says so; do not report only the head start.

## Workflow

1. **Validate**: confirm the target list carries `protein_id` and the bioactivity export carries all seven required columns. Refuse a missing column rather than defaulting it.
2. **Resolve identity**: take each target's accession, Pfam accession and organism from its columns, `--organism`, or the `protein_id`, and record which are absent.
3. **Attribute**: assign every record to at most one tier per target, testing exact_target, then ortholog, then family, then phenotypic.
4. **Summarise per tier**: count distinct compounds, record the best activity per unit family, and list activity types and source identifiers. Never add tiers together.
5. **Report**: write `report.md`, `result.json`, three tables and the reproducibility bundle. Mark every unassessable tier with its reason.

**Freedom level**: tiering is prescriptive and must not be adjusted per run.
The interpretation paragraph is where the agent reasons, and it must not
promote a phenotypic count into a claim about a protein.

## CLI Reference

```bash
# Standard usage
python skills/target-chemistry-evidence/target_chemistry_evidence.py \
  --targets target_scores.csv --bioactivity chembl_export.tsv \
  --orthologs orthologs.csv --organism "Candida albicans" \
  --output reports/chemistry

# Minimal: target list and export only; ortholog and family tiers report n/a
python skills/target-chemistry-evidence/target_chemistry_evidence.py \
  --targets target_scores.csv --bioactivity chembl_export.tsv \
  --output reports/chemistry

# Demo mode (synthetic data, no user files needed)
python skills/target-chemistry-evidence/target_chemistry_evidence.py \
  --demo --output /tmp/tce_demo

# Via ClawBio runner
python clawbio.py run target-chemistry --demo
```

## Demo

```bash
python clawbio.py run target-chemistry --demo
```

Expected output: 7 synthetic targets; 1 with exact-target evidence, 2
`phenotypic_only`, 1 `no_evidence`, and the ortholog tier `n/a` for the 4
targets absent from the mapping.

## Algorithm / Methodology

**Tier precedence** — a record is attributed to the first tier it qualifies for,
so no measurement is counted twice for the same target:

| Order | Tier | Condition |
|---|---|---|
| 1 | `exact_target` | record accession equals the target's accession |
| 2 | `ortholog` | record accession is in the target's ortholog mapping |
| 3 | `family` | record Pfam accession equals the target's Pfam accession |
| 4 | `phenotypic` | record has **no** accession and its organism matches the target's |

**Assessability.** Each tier needs one field: an accession, an ortholog
mapping, a Pfam accession, an organism. Missing the field makes the tier
`n/a` with a reason, which is not the same as the tier being evaluated and
found empty.

**Per tier, reported separately**: distinct compound count, record count, best
activity per unit family, activity types seen, count of measurements whose
units did not convert, and the source identifiers. **There is no total.** Each
target carries a `profile` label instead — `exact_target_only`,
`ortholog+family+phenotypic`, `no_evidence` — which names the tiers rather than
ranking them.

**Normalisation**: accessions are upper-cased and stripped of isoform and
version suffixes. Organisms match case-insensitively and tolerate a strain
suffix at a word boundary, so `Candida albicans SC5314` matches `Candida
albicans` and `Candida glabrata` does not. Molar units become nM; g/L, mg/L,
ug/mL and ng/mL become ug.mL-1. Unit strings are matched exactly and never
case-folded, because folding would map `mM` onto `M`.

**Ordering**: targets stay in input order. The input is already a ranking, and
re-sorting by chemistry here would silently replace the upstream priority with
a chemistry-driven one.

## Example Queries

- "Which of these 40 targets already have compounds with measured activity?"
- "I exported ChEMBL activities for Candida albicans — map them onto my target list"
- "Show me targets where the only evidence is a whole-organism MIC"
- "Which targets have human ortholog chemistry I should worry about for selectivity?"

## Example Output

```markdown
# Target Chemistry Evidence Report

**Input**: built-in synthetic demo data (7 targets, 24 fabricated bioactivity records)
**Organism**: Candida albicans

Checked **7 targets** against **24 bioactivity records**; **6** carry evidence in at least one tier.

| Target | Description | Exact | Ortholog | Family | Phenotypic | Profile |
|---|---|---|---|---|---|---|
| ERG11 | lanosterol 14-alpha-demethylase | 3 | 1 | 1 | 7 | exact_target+ortholog+family+phenotypic |
| FKS1 | beta-1,3-glucan synthase catalytic subunit | 0 | n/a | 0 | 7 | phenotypic_only |
| HSP90 | molecular chaperone | 0 | 2 | 1 | 7 | ortholog+family+phenotypic |
| ERG3 | C-5 sterol desaturase | 0 | n/a | 0 | 0 | no_evidence |

**Exact-target evidence**: 1 of 7 targets. Target-level annotation for fungal
proteins is sparse in public bioactivity databases, so an empty exact-target
tier is the common case and is not a failure of the search.

*ClawBio is a research tool. Not a medical device.*
```

## Output Structure

```
output_directory/
├── report.md                            # Per-target tiers and interpretation
├── result.json                          # Per-tier counts, best activities, reasons
├── tables/
│   ├── target_evidence_summary.csv      # One row per target: counts, profile, accession source
│   ├── evidence_by_tier.csv             # One row per target and tier
│   └── compound_evidence.csv            # One row per attributed measurement
└── reproducibility/
    ├── commands.sh                      # Exact command for this run
    ├── environment.yml                  # Conda environment snapshot
    └── checksums.sha256                 # sha256 of every output above
```

The bundle is written with `clawbio.common.reproducibility`
(`write_commands_sh`, `write_environment_yml`, `write_checksums`).

`evidence_by_tier.csv` is long rather than wide on purpose: every row names its
tier, so nothing downstream can add two tiers together without noticing.

## Dependencies

**Required**:
- Python 3.11+ standard library only. No third-party package, no external binary.

**Optional**:
- `curl` or any HTTP client, to produce the `--bioactivity` file. The skill
  never invokes one.

## Gotchas

- **The model will report a phenotypic count as chemistry against the protein.** Do not. A whole-organism MIC says a compound kills the organism; it says nothing about which of its thousands of proteins was hit. Quote the tier name every time a count is mentioned.
- **The model will read the phenotypic tier as target-specific.** It is not. Every target in the same organism receives the *same* organism-level records, because those records were never attributed to a protein. Seven targets each showing "7 phenotypic compounds" is one set of seven compounds, not forty-nine.
- **The model will treat an empty `exact_target` tier as a negative result.** It is an absence of deposited data. Fungal target-level annotation is sparse, and the common outcome of a correct run is zero. Say "no assay has been deposited", never "no compounds are active".
- **The model will want to combine the tiers into one druggability number.** Do not. That is the one thing this skill refuses to do, and adding it downstream destroys the distinction the whole output exists to preserve. A target with four exact-target IC50s and one with four whole-organism MICs are not interchangeable at any weighting.
- **The model will present ortholog chemistry as pure upside.** It is also the selectivity problem. A target whose ortholog tier is full is similar to a well-studied protein elsewhere, which is exactly the similarity that makes a selective compound hard. Report both directions or neither.
- **The model will fill in a blank `target_accession` to make a row "complete".** That silently promotes organism-level data into the exact-target tier. A blank accession is load-bearing data.
- **The model will compare a MIC in ug/mL to an IC50 in nM.** They are not comparable without a molecular weight the export does not carry. The skill keeps the two unit families apart and so must the narrative.
- **The model will try to fetch ChEMBL when the export is missing.** There is no network here. Give the user the URL patterns above and take the file as input.
- **The model will read `exact_target: 0` on an internal identifier as a real search.** It is not one. When `accession_source` is `protein_id fallback`, the accession was never a database accession, so nothing could have matched. The report flags it; do not summarise past the flag.
- **The model will read `n/a` as zero.** `n/a` means the target had no accession, no ortholog mapping, no Pfam accession or no organism, so the tier was never evaluated. The fix is a column, not a wider export.
- **The model will raise a potency cutoff to make a target look better.** The skill applies none; every record that attributes is counted. Filtering belongs in the export, where the user can see it.

## Safety

- **Research use only**: every report carries the ClawBio disclaimer.
- **No hallucinated science**: every count traces to a row in the supplied export; no activity is predicted, inferred or imputed.
- **No clinical claims**: an activity value supports no statement about efficacy or safety in a patient.
- **Local-first**: no network calls, no data upload, deterministic on repeated runs.

## Agent Boundary

The agent (LLM) dispatches, explains, and writes the interpretation. The skill
(Python) does every attribution. The agent must not merge tiers, invent an
overall score, fill in a missing accession, re-order the targets, or describe
phenotypic activity as evidence against a protein.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here on compound-evidence
language paired with a target list, or on a request mentioning ChEMBL and a
set of proteins together.

**Chaining partners**:
- `pathogen-target-dissimilarity`: produces the `--targets` file this skill consumes
- `struct-predictor`: targets with exact-target chemistry and no experimental structure
- `drug-repurposing-screen`: takes the compounds surfaced here into a viability screen
- `target-validation-scorer`: targets whose ortholog chemistry raises a host-protein question

## Maintenance

- **Review cadence**: quarterly, or on a ChEMBL release.
- **Staleness signals**: a ChEMBL API schema change to the activity or target
  endpoints, new unit strings appearing in exports, a new bulk download layout.
- **Deprecation**: retire if a maintained offline mirror of the chemistry
  databases becomes part of ClawBio, since attribution would then belong inside
  that component.

## Citations

- [The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods](https://academic.oup.com/nar/article/52/D1/D1180/7337608), *Nucleic Acids Research* 2024 — the source of the `--bioactivity` input
- [BindingDB in 2024: a FAIR knowledgebase of protein-small molecule binding data](https://academic.oup.com/nar/article/53/D1/D1633/7909300), *Nucleic Acids Research* 2025 — an alternative export producing the same columns
- [PubChem 2023 update](https://academic.oup.com/nar/article/51/D1/D1373/6777787), *Nucleic Acids Research* 2023 — BioAssay exports as a third input source
- [Pfam: The protein families database in 2021](https://academic.oup.com/nar/article/49/D1/D412/5943818), *Nucleic Acids Research* 2021 — the family tier's level of grouping
- [The antifungal pipeline: a reality check](https://pmc.ncbi.nlm.nih.gov/articles/PMC5760994/), *Nature Reviews Drug Discovery* 2017 — why most antifungal activity data is organism-level
- [Antifungal drug resistance: an emergent health threat](https://pmc.ncbi.nlm.nih.gov/articles/PMC8395843/), *Journal of Fungi* 2021 — target-level versus whole-cell evidence in antifungal discovery
