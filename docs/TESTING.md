[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Gallery](./GALLERY.md) | [Recipes](./RECIPES.md) | [FAQ](./FAQ.md) | [ABOUT](./ABOUT.md)

# Testing

This page defines the test suites used for routine checks and releases.

## Test Suites

### Fast (skip slow tests)
```bash
pytest tests/ -v -m "not slow"
```

### Full
```bash
pytest tests/ -v
```

### Regression (reference SVGs)
```bash
pytest tests/ -v -m "regression"
```

### Circular-only
```bash
pytest tests/ -v -m "circular"
```

### Linear-only
```bash
pytest tests/ -v -m "linear"
```

### Output Comparison
```bash
pytest tests/test_output_comparison.py::TestOutputComparison -v
```

### Regenerate Reference Outputs (intentional visual changes)
```bash
pytest tests/test_output_comparison.py::TestGenerateReferences -v
```

### Lint
```bash
ruff check gbdraw/ --select=E,F,W --ignore=E501,W503
```

[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Gallery](./GALLERY.md) | [Recipes](./RECIPES.md) | [FAQ](./FAQ.md) | [ABOUT](./ABOUT.md)
