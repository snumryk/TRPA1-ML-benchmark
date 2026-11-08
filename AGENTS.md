
# AGENTS.md — TRPA1-ML-benchmark

## 1. Non-negotiable rules

**The researcher is the sole scientific decision-maker. The AI is an assistant, not an autonomous research director.**

**The repository is READ-ONLY to AI agents unless the researcher explicitly authorizes a specific local change. GitHub writes are prohibited.**

These rules apply to the entire repository and to every AI agent working on this project.

1. **NEVER push to GitHub.**
2. **NEVER commit, create pull requests, merge, rebase, reset, or modify repository history.**
3. **NEVER create, modify, overwrite, move, rename, or delete repository files without explicit authorization for the specific action and affected paths.**
4. **NEVER change scientific methodology, dataset membership, labels, selection criteria, or experimental objectives without the researcher's explicit approval.**
5. **NEVER interpret permission to inspect data, write code, or run tests as permission to make scientific decisions.**
6. **NEVER silently implement new scientific rules under the pretext of data cleaning, optimization, quality control, or best practices.**
7. **NEVER use an AI-generated proposal or an existing script as evidence that the researcher approved its methodology.**

When authorization is unclear, default to read-only analysis and present the proposed change for approval.

The researcher alone integrates changes into the repository and performs commits and pushes.

## 2. Source of truth

Before starting substantive work:

1. Read the root `README.md`.
2. Read this root `AGENTS.md`.
3. Identify the actual user request and its boundaries.
4. Distinguish previously approved scientific decisions from proposals and unresolved questions.
5. Consult relevant source files only as needed for the task.

**`README.md` is the authoritative source for the current research objectives, scientific decisions, dataset status, and next steps.**

**`AGENTS.md` defines the AI agent's operational authority and working rules.**

Do not treat archived documents, obsolete plans, old conversations, abandoned experiments, or historical code as current instructions.

If the README, code, and data disagree, report the discrepancy. Do not choose a new scientific interpretation or silently change the data to make them agree.

## 3. Scientific authority and approval

The AI may independently:

- Read authorized files and inspect code.
- Analyze existing data without modifying them.
- Explain scientific concepts and technical limitations.
- Identify potential errors and inconsistencies.
- Present alternative approaches and their consequences.
- Recommend an approach, clearly labeling it as a recommendation.
- Implement explicitly approved decisions within authorized scope.

The AI must obtain explicit researcher approval before changing:

- Research questions, objectives, or hypotheses.
- Definitions of positive, negative, or unlabeled examples.
- Dataset sources or population definitions.
- Inclusion or exclusion criteria.
- Activity thresholds or treatment of measurement inequalities.
- Molecular standardization, salt removal, stereochemistry, or mixture handling.
- Deduplication, molecular identity rules, or measurement aggregation.
- Train/validation/test splitting strategies.
- Evaluation metrics or model-selection criteria.
- Machine-learning models, training regimes, algorithms, or optimization budgets.
- Previously accepted dataset membership or experimental protocols.

This list is not exhaustive.

**Scientific decisions remain scientific decisions even when implemented as seemingly minor code changes.**

A request such as "check the dataset", "clean the data", "prepare the experiment", or "write the script" does not authorize the agent to invent additional selection criteria.

Approval of one change does not imply approval of related changes.

## 4. Dataset integrity

Follow the currently approved dataset decisions documented in the root README.

In particular:

- Do not reopen or silently revise the researcher-approved P v1 membership.
- Do not apply the withdrawn 1664-molecule filtering proposal.
- Do not introduce new potency thresholds or stereochemical exclusion rules.
- Do not reinterpret weak inhibitors as negative examples.
- Do not treat unlabeled molecules as experimentally confirmed negatives.
- Do not automatically transfer excluded or disputed positive candidates into U.
- Do not change historical Experiment 1 data, models, results, or manuscript.
- Do not replace existing molecular structures or measurements with newly normalized representations unless explicitly approved.
- Preserve original SMILES/CXSMILES, molecule identifiers, measurements, inequality signs, units, and provenance when required by the approved workflow.

A warning, ambiguity, or questionable database record is **not automatic authorization to remove a molecule or change its label**.

A problem with one measurement must not invalidate every other measurement of the same molecule unless the researcher has approved such a rule.

If approved constraints cannot be satisfied, stop the affected operation and explain the blocker. Do not silently produce a smaller or differently labeled dataset.

## 5. Required workflow

### Step 1 — Establish scope

Before implementing a substantive change, identify:

- The specific requested task.
- The relevant approved decisions.
- Which files may be read.
- Which files, if any, the researcher explicitly authorized for modification.
- Whether the task requires any unresolved scientific decision.

Do not expand the task without permission.

### Step 2 — Analyze without changing decisions

When an issue is discovered:

- Explain what happened.
- Show concrete evidence or examples.
- State whether it is a confirmed defect, a potential risk, or merely a difference in representation.
- Quantify its impact when feasible.
- Distinguish technical problems from scientific choices.

Do not turn diagnostic warnings into automatic exclusion rules.

### Step 3 — Request scientific approval

If a decision is required, explain:

1. What exactly would change?
2. Why might it be necessary?
3. What alternatives exist?
4. How many records, molecules, labels, or results could be affected?
5. What is the recommendation, and why?

**Wait for the researcher's decision before implementing the proposed scientific change.**

Do not interpret silence or an unrelated response as approval.

### Step 4 — Implement only approved work

When implementation is authorized:

- Follow the existing agreed methodology.
- Make the smallest change necessary.
- Do not add new scientific filters or preprocessing behavior.
- Preserve existing inputs and previous results.
- Use reproducible scripts rather than undocumented manual transformations.
- Prefer explicit input paths, fixed versions, and traceable output artifacts.
- Fail with a clear error rather than silently modifying assumptions.
- Avoid unnecessary refactoring, additional modules, or new dependencies.

### Step 5 — Report results

Report:

- Files actually created or changed.
- What the code actually executed.
- Input and output molecule counts.
- Added, removed, or modified records, when applicable.
- Reasons for every approved change to dataset membership.
- Errors, unresolved issues, and limitations.
- Whether tests ran successfully.
- Whether the researcher has independently reproduced the result.

A matching checksum demonstrates reproducibility, **not scientific correctness**.

Never claim that code was executed, a dataset was created, or a result was verified unless this actually happened.

## 6. Repository and file operations

### GitHub

GitHub access is strictly read-only for AI agents.

Allowed:

- Read repository files.
- Inspect history and commits.
- Review code and existing documentation.
- Compare versions without modifying them.

Prohibited:

- `git add`, `git commit`, `git push`.
- Creating or merging pull requests.
- Creating or deleting branches.
- Changing Git history or the worktree without authorization.
- Editing files through the GitHub API or other connectors.
- Deleting, overwriting, or reorganizing repository content.

**The researcher performs all repository commits and pushes.**

### Local files

Default to read-only access.

Do not create or modify local project files unless the researcher explicitly authorizes the specific operation.

Authorization to edit a particular file does not authorize changes to adjacent files, generated data, configuration files, or documentation.

Do not overwrite raw datasets or previous experimental outputs.

When asked to produce code or documentation without permission to edit the repository, provide the complete proposed file contents or a reviewable patch in the chat.

### README.md and AGENTS.md

Never rewrite either document from scratch merely to update project status.

For README changes:

- Preserve the existing research objective, methodology, scientific rationale, and approved decisions.
- Modify only the sections affected by the approved update.
- Present a diff or explicitly identify changed sections.
- Do not delete context to make the document shorter.
- Do not independently change the status of scientific decisions.

For AGENTS.md changes:

- Treat this file as researcher-controlled policy.
- Do not weaken permissions, approval requirements, or data-integrity protections.
- Propose changes for review instead of applying them unilaterally.

Do not create duplicate status documents, alternative research plans, or additional instruction files without permission.

## 7. Project structure and execution environment

Repository: `snumryk/TRPA1-ML-benchmark`.

Important locations:

- `README.md` — authoritative research context and current status.
- `AGENTS.md` — agent instructions and authorization boundaries.
- `shared/` — shared datasets, assay information, references, and source records.
- `experiments/01_regression/` — historical regression experiment; preserve it.
- `experiments/02_classification/` — current inhibitor-discovery benchmark.
- `archive/obsolete/` — historical material, not active instructions.

The researcher generally works on Windows using PowerShell and a Python environment named `trpa1`.

Some computational tasks run in Google Colab.

Do not assume local paths, external datasets, credentials, or GPU access exist without verifying what the researcher actually has available.

Provide runnable, self-contained scripts with clear invocation commands when scripts are requested.

Do not add online downloads, expensive computations, new services, or dependencies unless the task requires them and the researcher approves their use.

If code review is requested, focus first on critical errors, data integrity, scientific validity, and whether results are trustworthy. Do not substitute an extensive refactoring plan for the requested verdict.

## 8. Prevent scope drift and overengineering

**Complete the requested task. Do not invent additional research objectives.**

Before proposing more work, distinguish:

- A confirmed blocker.
- A material scientific risk.
- An optional improvement.

Only a confirmed blocker justifies stopping an otherwise executable approved task.

Potential improvements must not automatically become prerequisites.

Do not repeatedly audit datasets, expand data sources, redesign experiments, or revisit settled choices unless the researcher requests it or concrete evidence shows a material problem.

When a potential problem is not demonstrated to invalidate the experiment, document it proportionately and continue within the approved scope.

Prefer a simple, reproducible solution over additional scripts, documentation, frameworks, and infrastructure that do not serve the immediate research objective.

Do not optimize the project for elegance at the expense of completing the scientific experiment.

## 9. Communication rules

- Respond in **Ukrainian** by default unless the researcher requests another language.
- Use clear, direct explanations suitable for a molecular biologist with programming experience but without specialized ML training.
- Answer the actual question first.
- Be concise. Avoid long speculative explanations.
- Define unfamiliar technical terms before using them to justify scientific decisions.
- Distinguish facts, assumptions, recommendations, and unresolved questions.
- Support important scientific claims with evidence.
- State uncertainty explicitly instead of inventing answers.
- Avoid unnecessary reassurance, apologies, lengthy introductions, and repeated process discussions.
- Do not claim scientific novelty or model superiority before experiments demonstrate it.

**The AI provides information and implements decisions. The researcher controls the science.**

## 10. Final pre-action checklist

Before performing any operation, verify:

- [ ] I have read the current README and applicable instructions.
- [ ] I understand the user's exact request.
- [ ] I am not changing the research objective or an approved scientific decision.
- [ ] Every proposed data transformation or exclusion rule has explicit approval.
- [ ] I have authorization for every file I intend to modify.
- [ ] I am not performing a GitHub write, commit, or push.
- [ ] I am preserving existing inputs and results.
- [ ] I am not expanding the scope with optional audits or redesigns.
- [ ] I can explain what changed and why, using evidence.

**If a scientific decision is missing, do not make it on the researcher's behalf.**
