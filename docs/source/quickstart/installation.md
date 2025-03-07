# Installation

|OS|Linux|Windows|Mac|
|:-:|:-:|:-:|:-:|
|**compatibility**|Yes|?|?|

## Install conda (virtual environment manager)

*We use conda environment to avoid version problem with other application dependencies.*

We recommend to download the lighter version `miniconda`.

- [Installing miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions)

- [Installing anaconda distribution](https://www.anaconda.com/docs/getting-started/anaconda/install#basic-install-instructions)

## Create conda environment

Open a **terminal** (for Windows user: from the Start menu, open the **Anaconda Prompt**). Create a conda environment and activate it:
```
conda create -n traceratops python=3.11
conda activate traceratops
```

## Install package

```bash
cd $HOME/Repositories/traceratops
pip install -e .
```

## Developer mode installation

```bash
cd $HOME/Repositories/traceratops
pip install -e ".[dev]"
```

### [Recommended] Setup pre-commit in local

- check if it's well installed

  `pre-commit --version`

- install command of the file ".pre-commit-config.yaml" inside ".git/hooks/pre-commit"

  `pre-commit install`

- test pre-commit without any commit

  `pre-commit run --all-files`

```{note}
- First time, it can take few minutes because it need to install all dependancies.
- Then, you can re-run this command when you want to simulate action that are execute, by hidden way, when you try to make a commit.
- For info, if one step of the pre-commit fail, your commit fail and you need to fix pre-commit error to be allow to commit.
```

- [Optional] fix strange issue or warning

  `pre-commit autoupdate --repo https://github.com/pre-commit/pre-commit-hooks`

- [Optional] Update pre-commit file
  - `pre-commit clean`
  - `pre-commit autoupdate`
  - `pre-commit install`
