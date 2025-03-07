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
