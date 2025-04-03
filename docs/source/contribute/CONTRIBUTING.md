# Contributing guide

Thank you for your interest in contributing to this project!
Please follow the guidelines below to ensure smooth collaboration and maintain a clean and structured codebase.

---

## Branching model

- `main`: The **stable production branch**. It contains only **squashed commits** that correspond to each official release.
- `dev`: The **development branch**. It contains **integrated features** and **bug fixes**, pending final validation through manual testing on real data before going into production.

---

## Workflow

### New branch

**Always create your working branch from `dev`**:

```bash
git checkout dev
git pull
git checkout -b <type>/<short_task_name>
```

```{note}
Branch names should follow this pattern:
`<type>/<short_task_name>`

Where `<type>` can be:
- `feat`: New feature
- `fix`: Bug fix
- `doc`: Documentation update
- `test`: Adding missing tests
- `chore`: Technical tasks with no impact on functionality (e.g., updating dependencies, CI config, build scripts, etc.)
```

### Commit messages

Use clear and conventional commit messages. We recommend following [Conventional Commits](https://www.conventionalcommits.org/):

Template:
```
<type>(<scope>): Summary title #<issue_number>

- List of task or bugfix done in this commit
```

Example:
```
doc(contributor): Write file for commit good practices #167

- Content reminder
- Message template
- This example
```

```{note}
During the Pull Request to `dev`, the "squash & merge" method will be used, so this clear message is only needed for the squash commit.
```

### Pull Request (PR)

Once your development is complete:

```bash
git push origin <your-branch-name>
```

Then:
- Open a **pull request to `dev`**
- A **PR template is provided** with a checklist to support the validation process
- Select **"Squash & merge"** when merging the PR
- Use a **clear and descriptive PR title** — it will become the commit message on `dev`

---

### Continuing development (while PR is pending)

```{warning}
Avoid working again on the same branch that is already under review.
```

While waiting for your PR to be reviewed and merged, you can continue development by creating a new branch from `dev` (recommended) or from your previous branch if you need the latest update to continue working.

```bash
git checkout dev
git pull
git checkout -b <type>/<new_task>
```

---

## Resolving merge conflicts

Merge conflicts may happen when:
- Two PRs modify the same part of the codebase.
- A long time has passed between when your branch was created and when it's merged.

To resolve conflicts:
1. Pull the latest `dev` into your branch:
```bash
git checkout <your-branch>
git pull origin dev
```

2. Git will notify you of conflicts — open the conflicting files, resolve manually, then:
```bash
git add <resolved-files>
git commit
```

3. Push the updated branch:
```bash
git push origin <your-branch>
```

If the PR is not closed, it will automatically integrate your new conflict resolution commit.



---

## Release

- This is done **via a pull request** from `dev` to `main`, and must use **squash & merge**.
- It's a **significant product update**.
- The commit message on `main` should clearly describe the release content, e.g.:
  `release: v1.3.0 - performance improvements & dashboard redesign`

After the merge:
- Tag the new release commit on `main` with a semantic version, e.g.: `v1.3.0`.

---

## Testing & Validation

- Make sure all tests pass before submitting a PR.
- Add tests for any new feature or fix when possible.
- Manual tests on real (big) data are expected before merging into `main`.

---

## Communication

- Use GitHub Issues for reporting bugs or proposing features.
- Use discussions or team communication tools (Zulip) for technical coordination.

---

**Thank you for contributing!**
