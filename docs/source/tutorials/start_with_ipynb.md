# Start with Jupyter Notebook

*A Jupyter Notebook is an interactive file where you can find both markdown text and executable code with its outputs displayed.*

*The main tutorials are in this format.*

*Locally, you can manipulate a Jupyter Notebook file in various ways: in your web browser using the basic ‘Jupyter Notebook’ or the modern ‘JupyterLab’, or in an IDE such as Spyder, VSCode or PyCharm. We recommend using JupyterLab to get started.*


## Installation of JupyterLab

1. Activate your [conda environment]. Hereafter we will assume your environment is called `my_env`. Replace this below if this is not the case.

```sh
conda activate my_env
```


```{note}
Replace `my_env` by your `environment name` if you called it something else.
```

2. Install a tool to manage jupyter notebook like [JupyterLab: <img src="../_static/Download-Icon.png" width="50"/>](https://jupyter.org/install#jupyterlab)
```sh
conda install jupyterlab
```

3. Create a specific kernel to ensure that the correct environment is used during code execution:

```
python -m ipykernel install --user --name my-kernel
```

```{note}
Replace `my-kernel` by a distinctive name.
```

## Open tutorial with JupyterLab

1. Find the path to your tutorial file, it's ending with `.ipynb` extension.

2. Open a terminal inside your tutorial folder and activate your [conda environment]):
```sh
conda activate my_env
```

3. Open your tutorial with JupiterLab:
```sh
jupyter-lab my_tutorial.ipynb
```

4. Once you spin up a jupyter lab, select the `my-kernel` (the kernel that you defined at step 3.) by clicking on panel "Kernel > Change Kernel..." to be able to run your specific functions.

5. Now you can follow the tutorial by running each cell with the `run` icon (or `Shift+Enter` on keyboard).
