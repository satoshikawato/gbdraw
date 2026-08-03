[Documentation home](../../DOCS.md) | [Python API](../../REFERENCE/python-api.md) | [Python how-to guides](README.md)

# How to draw Linear diagrams and comparisons from Python

Draw a whole-genome Lambda/DE3 comparison from fixed LOSATN evidence. The two
complete records occupy separate Linear rows.

## Prerequisites

- Install gbdraw in a Python 3.10 or newer environment.
- Start in an empty working directory.
- Copy these public fixtures into that directory:
  - [`NC_001416.gb`](../../../gbdraw/web/tutorial-data/lambda/NC_001416.gb),
    SHA-256 `4b76b8bacc8026aac3f19a4a915f4ac772ad61e7ec18f0e2cc859229f95a66e7`;
  - [`NC_042057.1.gb`](../../../gbdraw/web/tutorial-data/de3/NC_042057.1.gb),
    SHA-256 `288eb87480f8fe6eab6246fe1fc4af78c85a6cb7e591c47fd7d7f0170c932e09`;
  - [`lambda-de3.losatn.tsv`](../../../gbdraw/web/tutorial-data/lambda-de3-comparison/lambda-de3.losatn.tsv),
    SHA-256 `703b0ac749669152a8e6d5fa6fb246cf5973cb1a3e0ca9db2f346ee33628317c`.

The [fixture manifest](../../../gbdraw/web/tutorial-data/manifest.json) records
their provenance. Lambda (`NC_001416.1`, 48,502 bp) and DE3 (`NC_042057.1`,
42,925 bp) are complete genomes. `lambda-de3.losatn.tsv` is the fixed six-row
LOSATN result for that pair.

## Draw the comparison

Save this program as `draw_linear_examples.py` beside the three fixture files:

<!-- executable:H-PY-02:start -->
```python
from pathlib import Path

from Bio.SeqRecord import SeqRecord

from gbdraw import (
    Diagram,
    FeatureOptions,
    LinearComparisonOptions,
    LinearOptions,
    Thresholds,
    TitleOptions,
    draw_linear,
    read_genbank,
)


comparison_records = read_genbank(
    [Path("NC_001416.gb"), Path("NC_042057.1.gb")]
)
assert [(record.id, len(record)) for record in comparison_records] == [
    ("NC_001416.1", 48_502),
    ("NC_042057.1", 42_925),
]
assert all(isinstance(record, SeqRecord) for record in comparison_records)

comparison_options = LinearComparisonOptions(
    blast_files=("lambda-de3.losatn.tsv",),
)
comparison_diagram = draw_linear(
    comparison_records,
    options=LinearOptions(
        features=FeatureOptions(types=("CDS",)),
        comparisons=comparison_options,
        thresholds=Thresholds(evalue=1e-5),
        title=TitleOptions(text="Whole-genome Lambda and DE3 LOSATN comparison", position="top"),
        legend="none",
    ),
)
comparison_svg = comparison_diagram.to_svg()
comparison_bytes = comparison_diagram.to_bytes("svg")
comparison_path = comparison_diagram.save(Path("python_linear_comparison.svg"))

assert isinstance(comparison_diagram, Diagram)
assert comparison_diagram.mode == "linear"
assert comparison_svg.encode("utf-8") == comparison_bytes
assert comparison_path.read_bytes() == comparison_bytes

print("Saved python_linear_comparison.svg")
```
<!-- executable:H-PY-02:end -->

Run it from the working directory:

```bash
python draw_linear_examples.py
```

## Verification

The comparison has two complete record rows and six ribbons. Each ribbon is read
from the fixed LOSATN table; the program does not run a search.

![Whole-genome Lambda and DE3 Linear comparison from six LOSATN rows](../../images/h-py-02/python_linear_comparison.svg)

The recipe runner executes the literal code in a clean directory. It verifies the
two record identities and lengths, 130 displayed features, the two-row topology,
and all six LOSATN endpoint pairs.

## Troubleshooting

- `Could not find comparison file`: keep `lambda-de3.losatn.tsv` beside the
  program. It is a precomputed table, not an input genome.
- The comparison has fewer than six ribbons: retain the documented
  `Thresholds(evalue=1e-5)` setting. It includes every fixed row.
- The two record rows overlap: use the default Linear layout and keep both whole
  GenBank inputs in the documented order.

See the [Python API reference](../../REFERENCE/python-api.md) for Linear option
types and the [comparison result reference](../../REFERENCE/comparison-programs-thresholds-and-results.md)
for table fields and threshold meaning.
