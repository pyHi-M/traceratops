# 🚀 Release guide

*This guide describes the release process used in this project.*

| Component   | Purpose                                                      |
| ----------- | ------------------------------------------------------------ |
| `dev`       | Clean, linear development history – suitable for debugging.|
| `release/*` | Temporary branch used to stabilize and validate a release. |
| `main`      |Each release is merged via a **merge commit** for traceability. |
| **Tagging** | The tag is placed on the final squash commit in `dev`.|

------

## 1. Create a release branch from `dev`

```bash
git checkout dev
git pull origin dev
git checkout -b release/1.2.3
git push -u origin release/1.2.3
```

```{note}
This branch is used to validate the upcoming release. You can make fixes freely here.
```

------

## 2. QA & Final fixes on `release/1.2.3`

Process all [checks required to merge into `main`](https://traceratops.readthedocs.io/en/main/contribute/checklist/pr_template_COPY.html#main-checklist).

Make any necessary commits (fixes, version bump, etc.) directly on the `release` branch.

```bash
# example fix
git add .
git commit -m "Fix: final QA issue"
git push
```

------

## 3. Open a Pull Request: `release/1.2.3` → `dev`

- Use **"Squash and merge"** on GitHub.
- The result is a **single clean commit** on `dev` representing the release.
- Use a clear commit message like: `"Release 1.2.3"`

This keeps `dev` linear and avoids duplicate commits.

------

## 4. Reset `release/1.2.3` to `dev`

Now that `dev` contains the release, we align `release` with it to avoid SHA divergence:

```bash
git checkout dev
git pull origin dev

git checkout release/1.2.3
git reset --hard dev
git push --force
```

This ensures that `release` and `dev` point to the **exact same commit**.

------

## 4.(bis) Revalidate the release branch

```{note}
After resetting release/* to match dev, it’s important to revalidate the branch before merging into main — even if you already validated it earlier.
```

### Why is this necessary?

Between the time you:
1. created release/1.2.3
2. and squash-merged it into dev

👉 another developer may have merged a PR into dev, adding changes not present in the original release branch.

After the squash merge, when you reset --hard release ← dev, the release now includes those other commits, even though they were never validated in the original release context.

------

## 5. Tag the Final Release Commit (on `dev`)

Before merging into `main`, **tag the last squash commit on `dev`**, which represents the true released code:

```bash
git tag -a v1.2.3 -m "Release v1.2.3"
git push origin v1.2.3
```

> ✅ This makes the tag visible in both `dev` and `main` (since `main` merges the same commit from the release branch).

------

## 6. Open a Pull Request: `release/1.2.3` → `main`

- On GitHub, use **"Create a merge commit"**
- This gives a clear entry point in `main` for the release.
- Merge message:
   ➤ `"Release 1.2.3 – merged from release/1.2.3"`

This keeps `main` non-linear but readable and semantically meaningful.


------

## 7. Cleanup

You may delete the temporary release branch:

```bash
git branch -d release/1.2.3
git push origin --delete release/1.2.3
```
