# [How to] Import a script from pyHiM

## 1<sup>st</sup> - On traceratops
`cd $HOME/Repositories/traceratops`

1. Create a new branch:
	```bash
	git checkout dev
	git pull
	git checkout -b import/<script_name>
	```

2. Update duplicated code files (located inside `core/` folder). Example for `chromatin_trace_table.py`:

   1. ```bash
      cd $HOME/Repositories/pyHiM
      git checkout development
      git pull
      ```

   2. ```bash
      cp $HOME/Repositories/pyHiM/src/matrixOperations/chromatin_trace_table.py $HOME/Repositories/traceratops/traceratops/core/chromatin_trace_table.py
      ```

   3. Fix directly the pseudo merge conflicts via the "source control" of VSCode

   4. Run the all tests present on traceratops to ensure that you break nothing:

      ```bash
      cd $HOME/Repositories/traceratops
      conda activate traceratops
      pytest test -vv
      ```

   5. Make a commit with ONLY the update changes (should be about only one file).

3. Import the script from pyHiM (development branch)

4. **Tests** each option case (run args) & **Document** & **Refactor** (only the code that are tested)

5. Rebase on `dev` branch

6. Wait to check if the continuous integration succeed (pytest & pre-commit)

7. [Recommended] Make manual tests with real (big) data

8. Update MINOR number of the version (ex: 0.2.2 ==> 0.3.0) inside the `pyproject.toml`

9. Make a Pull Request on `main`

## 2<sup>nd</sup> - On pyHiM (development / branch tempo)
`cd $HOME/Repositories/pyHiM`

1. Create a new branch:

   ```bash
   git checkout development
   git pull
   git checkout -b export/<script_name>
   ```

2. Check one last time that this script has never been imported from another pyHiM file

3. [Recommended] Check that the script did not call external functions that it was the only one to call ==> in which case you need to delete them to clean up the code.

4. Delete exported script

5. Delete all references inside the documentation

6. Delete reference inside the `pyproject.toml`

7. Update version inside the `pyproject.toml` with the the "dev" keyword (ex: 1.0.0.dev1 ==> 1.0.0.dev2)

8. Make a Pull Request on `development`
