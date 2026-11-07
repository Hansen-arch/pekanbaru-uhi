# Conventions for this repository

Rules for anyone, human or agent, working here. They exist because earlier work
in this project broke them.

## Prose

- **No em dashes.** Use a comma, colon, full stop or parentheses instead.
- Write as the analyst who chose the method. Do not describe the work as an
  assignment: no "the brief asked for", "the requirement was", "as requested",
  "for the deliverable", "the task". State a constraint as a property of the
  data ("Level-2A begins in December 2018"), not as a workaround for a spec.
- No filler or hype: no "delve", "leverage", "comprehensive", "seamless", "unlock".
- Every number in the README must come from a file in `outputs/` and be
  reproducible by a script in `src/` or the repo root. If a number cannot be
  reproduced, say so in the text instead of quoting it.
- Report limitations honestly, including results that contradict the expected
  sign or the earlier version of the analysis.

## Git

- **Human authorship only.** Commits are authored by Hansen-arch. No tool credit
  lines, no assistant co-author trailers, no robot emoji in messages.
- Commit messages say why the change matters to a reader. Never write a message
  about a cosmetic edit ("remove em dashes", "fix typo", "tidy formatting");
  fold such a change into the commit it belongs to.
- **One repository per project.** Before creating a new repo, check
  `~/Projects` for an existing analysis of the same subject and extend that
  instead of building a parallel version.

## Data and secrets

- Never print a token, paste one into a command, or store one in a git remote
  URL. Read credentials from the environment or a credential helper.
- Do not commit regenerable large rasters; `.gitignore` covers them. Commit the
  small tables and the figures the README references.

## Verification

- Rerun the script that produces any figure or table you change, then open the
  PNG and confirm it renders. Do not claim a figure works because the script
  exited zero.
- `.githooks/pre-commit` rejects em dashes and `.githooks/commit-msg` rejects
  non-human authorship lines. Hooks are enabled with
  `git config core.hooksPath .githooks`.

## Hooks installed here

`.githooks/pre-commit` blocks em dashes and credential-shaped strings in staged
text files. `.githooks/commit-msg` blocks non-human authorship lines. Enable with
`git config core.hooksPath .githooks`.
