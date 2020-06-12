# Machine Learning Powered Riemann Solver

## Installation

Using `setyp.py`:

```bash
git clone git@gitlab.com:rhyme-org/ml_riemann_solver.git
cd MLRiemannSolver
python3 setup.py install # --user in case you want to install it locally
```

Using `pip`:

```bash
git clone git@gitlab.com:rhyme-org/ml_riemann_solver.git
cd MLRiemannSolver
pip install . # --user in case you want to install it locally
pip3 install . # --user in case you want to install it locally
```

To install the package locally in editable mode, run:

```bash
$ pip install -e . --user
$ pip3 install -e . --user
```

## Usage

```python
from ml_riemann_solver import MLRiemannSolver as MLRS
```

## Running tests

To run the tests, execute the following command:

```bash
$ python3 ./setup.py test
```

or

```bash
$ py.test # -s to show the stdout
```
