# tbabspy

## tbabs X-ray Absorption in Python

This is a very light and unofficial wrapper over Jörn Wilms' `tbabs` XSPEC model.

### Install

Clone the repo, `cd` into the directory, and just run `pip`:

```
$ git clone https://github.com/geoffryan/tbabspy.git
$ cd tbabspy
$ pip install -e .
```

### Run

`tbabspy` reproduces the `tbabs()`, `ztbabs()`, `tbfeo()`, and `tbgas()` models.  Each are top-level functions which take an array-like `Energy` (in keV) as their first argument and the standard parameters as their following arguments.

```python
>>> import numpy as np
>>> import tbabspy
>>> Ebins = np.linspace(0.3, 10, 1001)
>>> NH = 1.0   # in units of 1e22 cm^-2
>>> photabs = tbabspy.tbabs(Ebins, NH)
```

### Examples

Examples will be found in the `examples/` directory.
