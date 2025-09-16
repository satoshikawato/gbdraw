# Tutorial 3: Advanced Customization

[< Back to Tutorials](./)

**Goal**: Go beyond the basics to gain fine-grained control over plot aesthetics using configuration files for colors and labels.

---

## Part 1: Advanced Color Control

While `--palette` is great for general styling, you often need to highlight specific features. `gbdraw` offers two ways to do this with simple tab-separated value (TSV) files.

### Method 1: Override Default Colors (`-d`)

This method replaces the default color for an entire feature type (e.g., make all `CDS` features gray).

1.  **Create a default-override TSV file.** Let's call it `modified_defaults.tsv`. This file has two columns: `FeatureType` and `Color`.

    ```tsv
    # modified_defaults.tsv
    CDS	#d3d3d3
    rRNA	#a83232
    ```

2.  **Use it in your command.**
    ```bash
    gbdraw circular \
      --gbk your_genome.gbk \
      -d modified_defaults.tsv \
      -o plot_custom_defaults.svg
    ```
    In the output, all CDS features will be gray and all rRNA features will be dark red.

### Method 2: Feature-Specific Colors (`-t`)

This method is more powerful. It colors individual features that match a specific rule, such as a gene's product name. This is perfect for highlighting genes of interest.

1.  **Create a feature-specific color TSV file.** Let's call it `feature_colors.tsv`. The file has 5 columns: `FeatureType`, `Qualifier`, `RegexPattern`, `Color`, and `LegendLabel`.

    ```tsv
    # feature_colors.tsv
    # This colors all tyrosine recombinases red and all BIRP proteins yellow.
    CDS	product	tyrosine recombinase	red	Tyrosine Recombinase
    CDS	product	baculoviral IAP repeat-containing protein	yellow	BIRP
    ```
2.  **Combine both methods for maximum control.** Let's make all CDS features gray *except* for the specific ones we want to highlight.

    ```bash
    gbdraw circular \
      --gbk your_genome.gbk \
      -d modified_defaults.tsv \
      -t feature_colors.tsv \
      --legend right \
      -o plot_highlighted.svg
    ```
The result is a plot where most genes are gray, but specific genes of interest are colored and labeled in the legend.

![Plot with custom colors](https://github.com/satoshikawato/gbdraw/blob/main/examples/LC738868_middle_separate_strands.png)

---

## Part 2: Advanced Label Control

When using `--show_labels`, you can prevent visual clutter using whitelists and blacklists.

### Blacklisting Labels (`--label_blacklist`)

This is the most common use case: hiding uninformative labels like "hypothetical protein".

* **As a command-line argument:**
    ```bash
    gbdraw circular \
      --gbk your_genome.gbk \
      --show_labels \
      --label_blacklist "hypothetical protein,uncharacterized protein" \
      -o plot_filtered_labels.svg
    ```
* **As a file:** Create a file `blacklist.txt` with one term per line, then use `--label_blacklist blacklist.txt`.

### Whitelisting Labels (`--label_whitelist`)

This is the opposite: you specify exactly which genes should be labeled, and all others will be hidden. This is useful for showing only a few key genes in a large genome.

1.  **Create a whitelist TSV file.** The format is `FeatureType`, `Qualifier`, `RegexPattern`.
    ```tsv
    # whitelist.tsv
    # Only show labels for these two genes
    CDS	locus_tag	ABC_00123
    CDS	gene	dnaA
    ```
2.  **Use it in your command.**
    ```bash
    gbdraw circular \
      --gbk your_genome.gbk \
      --show_labels \
      --label_whitelist whitelist.tsv \
      -o plot_whitelisted.svg
    ```

### Changing Label Content (`--qualifier_priority`)

By default, `gbdraw` uses the `product` qualifier for labels. If you prefer to use the `gene` name or `locus_tag`, you can specify a new priority.

1.  **Create a priority file.**
    ```tsv
    # priority.tsv
    # For CDS features, use the 'gene' qualifier first.
    CDS	gene
    ```
2.  **Use it in your command.**
    ```bash
    gbdraw circular \
      --gbk your_genome.gbk \
      --show_labels \
      --qualifier_priority priority.tsv \
      -o plot_gene_labels.svg
    ```
    