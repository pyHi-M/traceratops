
# Installation [DEV]

## Create conda environment

```bash
conda create -n traceratops python=3.11
conda activate traceratops
```

## Download & install the source code

```bash
cd $HOME/Repositories
git clone git@github.com:pyHi-M/traceratops.git

cd $HOME/Repositories/traceratops
pip install -e ".[dev]"
```

## [Recommended] Setup pre-commit in local

- check if it's well installed

  `pre-commit --version`

- install command of the file ".pre-commit-config.yaml" inside ".git/hooks/pre-commit"

  `pre-commit install`

- test pre-commit without any commit

  `pre-commit run --all-files`

```{note}
- First time, it can take few minutes because it need to install all dependencies.
- Then, you can re-run this command when you want to simulate action that are execute, by hidden way, when you try to make a commit.
- For info, if one step of the pre-commit fail, your commit fail and you need to fix pre-commit error to be allow to commit.
```

### Tips

- [Optional] fix strange issue or warning

  `pre-commit autoupdate --repo https://github.com/pre-commit/pre-commit-hooks`

- [Optional] Update pre-commit file
  - `pre-commit clean`
  - `pre-commit autoupdate`
  - `pre-commit install`

### Codespell

- [Very optional] Activate locally codespell inside .pre-commit-config.yaml
  - Uncomment codespell section
  - Ignore locally modification on the `.pre-commit-config.yaml` file:
    ```git update-index --assume-unchanged .pre-commit-config.yaml```
  - To undo: ```git update-index --no-assume-unchanged .pre-commit-config.yaml```
  - to see 'assume unchanged' files ("-v" option ==> use lowercase letters)
    ```git ls-files -v```

## Strange error

- If you get things like this during the installation:

```rest
ERROR: pip's dependency resolver does not currently take into account all the packages that are installed. This behaviour is the source of the following dependency conflicts.
pylint 2.12.2 requires isort<6,>=4.2.5, which is not installed.
pylint 2.12.2 requires mccabe<0.7,>=0.6, which is not installed.
pylint 2.12.2 requires toml>=0.9.2, which is not installed.
```

**OR**

- If a script name seems to refer to another package:

  - Identify the error package name (here it's pylint).
  - If it's not in the dependancies of your project, it's a "ghost".
  - To exorcise it:
    - Identify exe path : `which pylint` (give me `/home/user/.local/bin/pylint`)
    - Remove it : `rm ~/.local/bin/pylint`
  - Re-install your package to be sure `pip install -e ".[dev]"`
