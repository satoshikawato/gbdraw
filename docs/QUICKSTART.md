[Documentation home](./DOCS.md) | [Tutorials](./TUTORIALS/README.md) | [Technical documentation](./REFERENCE/README.md) | [FAQ](./FAQ.md) | [Gallery](./GALLERY.md) | [Installation](./INSTALL.md)

# Quickstart

To draw without installing anything, open the web app at
[https://gbdraw.app/](https://gbdraw.app/).

With gbdraw installed and a source checkout
(`git clone https://github.com/satoshikawato/gbdraw.git`), draw two figures
from the inputs in its `gbdraw/web/tutorial-data/` folder: a circular map of
the human mitochondrial genome (NCBI `NC_012920.1`), and the phage Lambda
(`NC_001416.1`) and DE3 (`NC_042057.1`) genomes joined by their LOSATN
comparison. Run in the directory that holds the `gbdraw` checkout:

```bash
data=gbdraw/gbdraw/web/tutorial-data
gbdraw circular \
  --gbk "$data/human-mitochondrion/HmmtDNA.gbk" \
  --separate_strands \
  --track_type middle \
  --labels out \
  --species "<i>Homo sapiens</i>" \
  -o human_mitochondrion \
  -f svg
gbdraw linear \
  --gbk "$data/lambda/NC_001416.gb" "$data/de3/NC_042057.1.gb" \
  --blast "$data/lambda-de3-comparison/lambda-de3.losatn.tsv" \
  --separate_strands \
  --show_labels first \
  -o lambda_de3 \
  -f svg
```

gbdraw writes `human_mitochondrion.svg` and `lambda_de3.svg`.

To build a figure step by step, follow one of the two first Tutorials. Each
page lists the files to download, then gives the steps for the web app, the
command line, and Python.

- [Draw a labeled circular map of the human mitochondrial genome](./TUTORIALS/first-circular-genome-diagram.md):
  [web app](./TUTORIALS/first-circular-genome-diagram.md#in-the-web-app),
  [command line](./TUTORIALS/first-circular-genome-diagram.md#on-the-command-line),
  or [Python](./TUTORIALS/first-circular-genome-diagram.md#in-python).
- [Draw a labeled linear map of the Lambda phage genome](./TUTORIALS/first-linear-genome-diagram.md):
  [web app](./TUTORIALS/first-linear-genome-diagram.md#in-the-web-app),
  [command line](./TUTORIALS/first-linear-genome-diagram.md#on-the-command-line),
  or [Python](./TUTORIALS/first-linear-genome-diagram.md#in-python).

If this is your first Tutorial, read [Get the tutorial inputs](./GETTING_TUTORIAL_DATA.md)
before downloading the inputs. You do not need a source checkout.

After your first diagram, choose another [Tutorial](./TUTORIALS/README.md). Use
the [FAQ](./FAQ.md) to choose an approach and [Technical
documentation](./REFERENCE/README.md) for exact controls, options, schemas,
APIs, and output behavior.

[Documentation home](./DOCS.md) | [Tutorials](./TUTORIALS/README.md) | [Technical documentation](./REFERENCE/README.md) | [FAQ](./FAQ.md) | [Gallery](./GALLERY.md) | [Installation](./INSTALL.md)
