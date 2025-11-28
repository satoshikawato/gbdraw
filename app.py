#!/usr/bin/env python
# coding: utf-8

import subprocess
import os
import io
import re 
import logging
import shutil
import uuid
import tomllib
import time
import tempfile
import gbdraw.version
import streamlit as st
from pathlib import Path
from importlib import resources
from contextlib import redirect_stdout, redirect_stderr
from streamlit.runtime.uploaded_file_manager import UploadedFile
import streamlit.components.v1 as components
from gbdraw.circular import circular_main
from gbdraw.linear import linear_main


# --- Google Analytics settings ---
GA_ID = "G-GG6JMKM02Y"
GA_JS = f"""
<script async src="https://www.googletagmanager.com/gtag/js?id={GA_ID}"></script>
<script>
    window.dataLayer = window.dataLayer || [];
    function gtag(){{dataLayer.push(arguments);}}
    gtag('js', new Date());
    gtag('config', '{GA_ID}');
</script>
"""

components.html(GA_JS, height=0, width=0)

# --- Basic Application Settings ---
st.set_page_config(
    page_title="gbdraw Web App",
    page_icon="🧬",
    layout="wide", menu_items={
        'Get help': 'https://github.com/satoshikawato/gbdraw/blob/main/docs/DOCS.md',
        'Report a bug': "https://github.com/satoshikawato/gbdraw/issues",
        'About': "# 🧬 gbdraw Web App\nA genome diagram generator for microbes and organelles.\nhttps://github.com/satoshikawato/gbdraw/"
    })

# --- Cookie Consent Banner ---
if 'cookie_consent' not in st.session_state:
    cookie_container = st.container()
    
    st.markdown("""
        <style>
        .cookie-banner {
            position: fixed;
            bottom: 0;
            left: 0;
            width: 100%;
            background-color: #f0f2f6;
            padding: 1rem;
            z-index: 100000;
            border-top: 1px solid #ccc;
            text-align: center;
            box-shadow: 0px -2px 10px rgba(0,0,0,0.1);
        }
        </style>
    """, unsafe_allow_html=True)

    with cookie_container:
        st.markdown('<div class="cookie-banner">', unsafe_allow_html=True)
        
        cols = st.columns([0.8, 0.2])
        with cols[0]:
            st.write("🍪 **Privacy Notice:** This website uses Google Analytics to collect anonymous usage data (e.g., visitor country) to improve the tool. No personal genomic data is stored.")
        with cols[1]:
            if st.button("Got it!", type="primary", key="cookie_btn"):
                st.session_state['cookie_consent'] = True
                st.rerun()
        
        st.markdown('</div>', unsafe_allow_html=True)

def get_version_info():
    """Retrieves the gbdraw version and commit ID."""
    try:
        version = gbdraw.version.__version__
    except AttributeError:
        version = "N/A"
    
    try:
        commit_id = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"]
        ).strip().decode("utf-8")
    except (subprocess.CalledProcessError, FileNotFoundError):
        commit_id = "N/A"
        
    return version, commit_id

VERSION, COMMIT_ID = get_version_info()

st.title("🧬 gbdraw Web App")
st.caption("A genome diagram generator for microbes and organelles")

st.markdown(
    "📖 **Documentation:** See the [**Official Documentation**](https://github.com/satoshikawato/gbdraw/blob/main/docs/DOCS.md) for detailed usage and examples."
)

# --- Define the list of selectable feature keys ---
FEATURE_KEYS = [
    "assembly_gap", "C_region", "CDS", "centromere", "D-loop", "D_segment",
    "exon", "gap", "intron", "J_segment", "mat_peptide", "misc_binding",
    "misc_difference", "misc_feature", "misc_RNA", "misc_structure",
    "mobile_element", "modified_base", "mRNA", "ncRNA", "operon", "oriT",
    "precursor_RNA", "primer_bind", "propeptide", "protein_bind", "regulatory",
    "repeat_region", "rep_origin", "rRNA", "sig_peptide", "stem_loop",
    "telomere", "tmRNA", "transit_peptide", "tRNA", "unsure", "V_region",
    "V_segment", "variation", "3'UTR", "5'UTR"
]
QUALIFIER_KEYS = ["allele", "anticodon", "artificial_location", "bound_moiety", "codon_start", "direction", "EC_number", "estimated_length", "exception", "experiment", "frequency", "function", "gap_type", "gene", "gene_synonym", "inference", "linkage_evidence", "locus_tag", "mobile_element_type", "mod_base", "ncRNA_class", "note", "number", "operon", "PCR_conditions", "product", "pseudo", "pseudogene", "regulatory_class", "replace", "ribosomal_slippage", "rpt_family", "rpt_type", "rpt_unit_seq", "satellite", "tag_peptide", "translation", "transl_except", "transl_table", "trans_splicing"]

# --- Helper functions and Session State for Dynamic Priority Input ---

def sanitize_filename(filename: str) -> str:
    """
    Strips directory traversal characters and removes characters that are unsafe
    for filenames.
    """
    # 1. Strip path traversal characters
    sanitized = os.path.basename(filename)
    
    # 2. Remove illegal characters for most filesystems
    #    (allows letters, numbers, underscore, hyphen, and dot)
    sanitized = re.sub(r'[^a-zA-Z0-9_.-]', '', sanitized)
    
    # 3. If the name is empty after sanitizing, provide a default
    if not sanitized:
        return ""
        
    return sanitized

def add_priority_row():
    """Appends a new empty row to the list in session_state."""
    if st.session_state.manual_priorities:
        new_id = max(row['id'] for row in st.session_state.manual_priorities) + 1
    else:
        new_id = 0
    st.session_state.manual_priorities.append({'id': new_id, 'feature': '', 'qualifiers': ''})

def remove_priority_row(row_id):
    """Removes a specific row by its ID."""
    st.session_state.manual_priorities = [
        row for row in st.session_state.manual_priorities if row['id'] != row_id
    ]

# --- Whitelist Rows Functions ---
def add_whitelist_row():
    """Appends a new empty row to the whitelist in session_state."""
    if st.session_state.manual_whitelist:
        new_id = max(row['id'] for row in st.session_state.manual_whitelist) + 1
    else:
        new_id = 0
    st.session_state.manual_whitelist.append({'id': new_id, 'feature': 'CDS', 'qualifier': 'product', 'keyword': ''})

def remove_whitelist_row(row_id):
    """Removes a specific whitelist row by its ID."""
    st.session_state.manual_whitelist = [
        row for row in st.session_state.manual_whitelist if row['id'] != row_id
    ]

# Initialize session state for manual priorities if it doesn't exist.
# This runs only once per session.
if 'manual_priorities' not in st.session_state:
    st.session_state.manual_priorities = [{'id': 0, 'feature': 'CDS', 'qualifiers': 'product,gene'}]
if 'manual_whitelist' not in st.session_state:
    st.session_state.manual_whitelist = []


st.markdown(
    """
    <p>
    <a href="https://anaconda.org/bioconda/gbdraw"><img src="https://anaconda.org/bioconda/gbdraw/badges/version.svg" alt="version"></a>
    <a href="https://anaconda.org/bioconda/gbdraw"><img src="https://anaconda.org/bioconda/gbdraw/badges/platforms.svg" alt="platforms"></a>
    <a href="http://bioconda.github.io/recipes/gbdraw/README.html"><img src="https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat" alt="install with bioconda"></a>
    <a href="https://anaconda.org/bioconda/gbdraw"> <img src="https://anaconda.org/bioconda/gbdraw/badges/downloads.svg" /> </a>
    <a href="https://anaconda.org/bioconda/gbdraw"><img src="https://anaconda.org/bioconda/gbdraw/badges/license.svg" alt="license"></a>
    <a href="https://deepwiki.com/satoshikawato/gbdraw"><img src="https://deepwiki.com/badge.svg" alt="Ask DeepWiki"></a>
    </p>
    """,
    unsafe_allow_html=True
)

# --- Temporary Directory and Session State Initialization ---
# Upload destination for input files
TEMP_DIR = Path("gbdraw_temp")
UPLOAD_DIR = TEMP_DIR / "uploads"
UPLOAD_DIR.mkdir(exist_ok=True, parents=True)

# Initialize session state
if 'uploaded_files' not in st.session_state:
    st.session_state.uploaded_files = {}
if 'circular_result' not in st.session_state:
    st.session_state.circular_result = None
if 'linear_result' not in st.session_state:
    st.session_state.linear_result = None
if 'linear_seq_count' not in st.session_state:
    st.session_state.linear_seq_count = 1
if 'color_picker_seed' not in st.session_state:
    st.session_state.color_picker_seed = 0

# --- Helper Functions ---
@st.cache_data
def get_palettes():
    """Dynamically get the list of color palettes from gbdraw's internal files."""
    try:
        with resources.files("gbdraw").joinpath("data").joinpath("color_palettes.toml").open("rb") as fh:
            doc = tomllib.load(fh)
        return [""] + sorted(k for k in doc if k != "title")
    except (FileNotFoundError, ModuleNotFoundError, AttributeError):
        st.warning("Could not dynamically load palettes from gbdraw. Using a default list.")
        return ["default"]

@st.cache_data
def get_palette_colors(palette_name: str) -> dict:
    """Loads the colors for a specific palette from the internal TOML file."""
    try:
        with resources.files("gbdraw").joinpath("data").joinpath("color_palettes.toml").open("rb") as fh:
            all_palettes = tomllib.load(fh)
        return all_palettes.get(palette_name, all_palettes.get("default", {}))
    except (FileNotFoundError, ModuleNotFoundError, AttributeError) as e:
        st.warning(f"Could not dynamically load palette colors: {e}. Using a default list.")
        return {
            "CDS": "#89d1fa", "rRNA": "#71ee7d", "tRNA": "#e8b441",
            "tmRNA": "#ded44e", "ncRNA": "#c4fac3", "repeat_region": "#d3d3d3",
            "misc_feature": "#d3d3d3", "default": "#d3d3d3", "skew_high": "#6dded3",
            "skew_low": "#ad72e3", "gc_content": "#a1a1a1", "pairwise_match": "#d3d3d3"
        }

def manual_select_update(key):
    """Callback to update session state from a widget's value."""
    manual_key = f"{key}_manual"
    if key in st.session_state and manual_key in st.session_state:
        if st.session_state[key] != st.session_state[manual_key]:
            st.session_state[manual_key] = st.session_state[key]

def create_manual_selectbox(label, options, key, help=None):
    """Creates a selectbox with robust, manual state management."""
    manual_key = f"{key}_manual"
    if manual_key not in st.session_state:
        st.session_state[manual_key] = ""

    try:
        current_index = options.index(st.session_state[manual_key])
    except ValueError:
        current_index = 0
    
    st.selectbox(
        label,
        options=options,
        key=key,
        index=current_index,
        on_change=manual_select_update,
        help=help,
        args=(key,)
    )

PALETTES = get_palettes()

# --- Sidebar (File Management) ---
with st.sidebar:
    st.header("📂 File Management")
    uploaded_files_list = st.file_uploader(
        "Upload GenBank, GFF3, FASTA, Color, or BLAST files",
        accept_multiple_files=True
    )
    if uploaded_files_list:
        for uploaded_file in uploaded_files_list:
            safe_name = uploaded_file.name.replace(" ", "_").replace("(", "").replace(")", "")
            if safe_name not in st.session_state.uploaded_files:
                save_path = UPLOAD_DIR / f"{uuid.uuid4().hex[:8]}_{safe_name}"
                with open(save_path, "wb") as f:
                    f.write(uploaded_file.getbuffer())
                st.session_state.uploaded_files[safe_name] = str(save_path)
        st.success(f"{len(uploaded_files_list)} file(s) uploaded successfully!")

    st.subheader("Uploaded Files")
    if not st.session_state.uploaded_files:
        st.info("No files uploaded yet.")
    else:
        for fname in st.session_state.uploaded_files.keys():
            st.write(f"- `{fname}`")

    if st.button("Clear All Uploaded Files"):
        for key in list(st.session_state.keys()):
            del st.session_state[key]
        if UPLOAD_DIR.exists():
            shutil.rmtree(UPLOAD_DIR)
            UPLOAD_DIR.mkdir(exist_ok=True, parents=True)
        st.success("All uploaded files and states have been cleared.")
        st.rerun()
    st.markdown("---")
    with st.expander("📜 Privacy & Terms"):
        st.markdown("""
        **Privacy Policy**
        Uploaded data is processed temporarily for visualization and is automatically deleted immediately after processing. We do not store or view your data.

        **Analytics**
        This tool uses **Google Analytics** to collect anonymous usage statistics (e.g., visitor country) to improve the service. No personal data is linked to these statistics.
        
        **Disclaimer**
        This tool is provided "as is" without warranty of any kind. The author is not liable for any damages arising from its use.
        
        **License**
        Source code is available under the [MIT License](https://github.com/satoshikawato/gbdraw).
        """)
    with st.expander("📬 Contact"):
            st.markdown("""
            **Bug Reports & Feature Requests:**
            Please open an issue on [GitHub Issues](https://github.com/satoshikawato/gbdraw/issues).
            
            **Other Inquiries:**
            Satoshi Kawato (Tokyo Univ. of Marine Science and Tech.)
            Email: `kawato [at] kaiyodai.ac.jp`
            """)
    with st.expander("📚Citation"):
            st.markdown("""
            If you use gbdraw in your research, I would appreciate it if you could cite the official GitHub repository:
            [https://github.com/satoshikawato/gbdraw/](https://github.com/satoshikawato/gbdraw/)
            """)
# --- Main Content (Mode Selection) ---
file_options = [""] + sorted(st.session_state.uploaded_files.keys())

selected_mode = st.radio(
    "Select Mode",
    ["🔵 Circular", "📏 Linear"],
    horizontal=True,
    label_visibility="collapsed"
)

# --- CIRCULAR MODE ---
if selected_mode == "🔵 Circular":
    st.header("Circular Mode")
    st.markdown("---")
    st.subheader("Input Genome File")

    c_input_type = st.radio(
        "Input file type",
        ("GenBank", "GFF3 + FASTA"),
        key="c_input_type",
        horizontal=True,
    )

    if c_input_type == "GenBank":
        create_manual_selectbox("GenBank file:", file_options, "c_gb")
    else:
        col1, col2 = st.columns(2)
        with col1:
            create_manual_selectbox("GFF3 file:", file_options, "c_gff")
        with col2:
            create_manual_selectbox("FASTA file:", file_options, "c_fasta")

    st.markdown("---")
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("🎨 Color Customization")
        create_manual_selectbox("Default-override color file (optional):", file_options, "c_d_color", help="Tab-separated value (TSV) file that overrides the color palette. Feature types not specified will use the selected palette colors. See [here](https://github.com/satoshikawato/gbdraw/blob/main/docs/TUTORIALS/3_Advanced_Customization.md#method-1-override-default-colors--d) for details.")
        create_manual_selectbox("Feature-specific color file (optional):", file_options, "c_t_color", help="TSV file that overrides the color palette for specific features. Features not specified will use the selected palette colors or the default-override colors if provided. See [here](https://github.com/satoshikawato/gbdraw/blob/main/docs/TUTORIALS/3_Advanced_Customization.md#method-2-feature-specific-colors--t) for details.")

        # Color Palette Selection and Customization
        def circular_palette_changed():
            new_palette = st.session_state.c_palette_selector
            st.session_state.custom_circular_colors = get_palette_colors(new_palette).copy()
            st.session_state.color_picker_seed += 1

        if 'custom_circular_colors' not in st.session_state:
            default_palette = PALETTES[0] if PALETTES and PALETTES[0] != "" else "default"
            st.session_state.custom_circular_colors = get_palette_colors(default_palette).copy()

        st.selectbox(
            "Base color palette:", 
            PALETTES, 
            key="c_palette_selector", index=11,
            on_change=circular_palette_changed,
            help="Select a base color palette. You can further customize individual feature colors below. See [here](https://github.com/satoshikawato/gbdraw/blob/main/examples/color_palette_examples.md) for palette examples.",
        )

        st.write("Click on the color boxes below to customize the default colors for each feature.")
        cols = st.columns(5)
        color_keys = sorted(st.session_state.custom_circular_colors.keys())
        for i, feature in enumerate(color_keys):
            col = cols[i % 5]
            with col:
                st.session_state.custom_circular_colors[feature] = st.color_picker(
                    label=feature,
                    value=st.session_state.custom_circular_colors[feature],
                    key=f"c_color_picker_{feature}_{st.session_state.color_picker_seed}"
                )
    with col2:
    # --- Qualifier Priority & Label Filtering Section (OUTSIDE the form) ---
        st.subheader("✏️ Label Content Customization")
        
        # Qualifier Priority
        st.markdown("##### Qualifier Priority")
        st.markdown("Select a TSV file from the sidebar OR define priorities manually below. If a file is selected, manual entries are ignored.")
        create_manual_selectbox(
            "Qualifier Priority File (optional)",
            file_options,
            "c_qualifier_priority_file",
            help="Tab-separated value (TSV) file defining qualifier priorities for feature labels. See [here](https://github.com/satoshikawato/gbdraw/blob/main/docs/TUTORIALS/3_Advanced_Customization.md#changing-label-content---qualifier_priority) for details."
        )
        with st.expander("Or, define priorities manually:", expanded=False):
            col1, col2, _ = st.columns([3, 5, 1])
            col1.markdown("**Feature Type**")
            col2.markdown("**Qualifier Priority (comma-separated)**")

            for row in st.session_state.manual_priorities:
                col1, col2, col3 = st.columns([3, 5, 1])
                col1.text_input("Feature", value=row['feature'], key=f"feature_{row['id']}", label_visibility="collapsed", help="Feature type (e.g., CDS, tRNA, rRNA, misc_feature). See [here](https://docs.google.com/spreadsheets/d/1qosakEKo-y9JjwUO_OFcmGCUfssxhbFAm5NXUAnT3eM/edit?gid=0#gid=0) for a list of common feature types and asssociated qualifiers.")
                col2.text_input("Qualifiers", value=row['qualifiers'], key=f"qualifiers_{row['id']}", label_visibility="collapsed", help="Type of qualifiers to use for label content, in order of priority (e.g., product,gene,locus_tag).See [here](https://docs.google.com/spreadsheets/d/1qosakEKo-y9JjwUO_OFcmGCUfssxhbFAm5NXUAnT3eM/edit?gid=0#gid=0) for a list of common feature types and asssociated qualifiers.")
                col3.button("➖", key=f"remove_{row['id']}", on_click=remove_priority_row, args=(row['id'],))

            st.button("➕ Add Row", on_click=add_priority_row, use_container_width=True)

        # Label Content Filtering
        st.markdown("##### Label Content Filtering")
        c_filter_mode = st.radio(
            "Select label filtering mode:",
            ("None", "Blacklist (exclude keywords)", "Whitelist (include keywords)"),
            key="c_filter_mode",
            horizontal=True, help="Choose how to filter feature labels based on keywords. See [here](https://github.com/satoshikawato/gbdraw/blob/main/docs/TUTORIALS/3_Advanced_Customization.md#part-2-advanced-label-control) for details.")


        if c_filter_mode == "Blacklist (exclude keywords)":
            st.info("Select a file with blacklist keywords (one per line) OR enter them manually below. If a file is selected, manual entry is ignored. See [here](https://github.com/satoshikawato/gbdraw/blob/main/docs/TUTORIALS/3_Advanced_Customization.md#part-2-advanced-label-control) for details.")
            create_manual_selectbox(
                "Blacklist File (optional)",
                file_options,
                "c_blacklist_file"
            )
            st.text_area(
                "Or, enter keywords manually (comma-separated):",
                value="hypothetical, uncharacterized, putative, unknown",
                help="Features with these keywords in their labels will be hidden.",
                key="c_blacklist_manual"
            )
        elif c_filter_mode == "Whitelist (include keywords)":
            st.info("Select a file with label whitelist (one per line) OR enter them manually below. If a file is selected, manual entry is ignored. See [here](https://github.com/satoshikawato/gbdraw/blob/main/docs/TUTORIALS/3_Advanced_Customization.md#whitelisting-labels---label_whitelist) for details.")        
            create_manual_selectbox(
                "Whitelist File (optional)",
                file_options,
                "c_whitelist_file"
            )
            st.info("Define rules to ONLY show labels containing specific keywords. For example, show 'CDS' features where the 'product' contains 'DNA polymerase'.")
            with st.container():
                wl_col1, wl_col2, wl_col3, _ = st.columns([3, 3, 4, 1])
                wl_col1.markdown("**Feature Type**")
                wl_col2.markdown("**Qualifier**")
                wl_col3.markdown("**Keyword to Include**")

                for row in st.session_state.manual_whitelist:
                    wl_col1, wl_col2, wl_col3, wl_col4 = st.columns([3, 3, 4, 1])
                    row['feature'] = wl_col1.selectbox("Feature", options=FEATURE_KEYS, index=FEATURE_KEYS.index(row['feature']) if row['feature'] in FEATURE_KEYS else 0, key=f"wl_feature_{row['id']}", label_visibility="collapsed")
                    row['qualifier'] = wl_col2.selectbox("Qualifier", options=QUALIFIER_KEYS, index=QUALIFIER_KEYS.index(row['qualifier']) if row['qualifier'] in QUALIFIER_KEYS else 0, key=f"wl_qualifier_{row['id']}", label_visibility="collapsed")
                    row['keyword'] = wl_col3.text_input("Keyword", value=row['keyword'], key=f"wl_keyword_{row['id']}", label_visibility="collapsed")
                    wl_col4.button("➖", key=f"wl_remove_{row['id']}", on_click=remove_whitelist_row, args=(row['id'],))

            st.button("➕ Add Whitelist Row", on_click=add_whitelist_row, use_container_width=True)
    st.markdown("---")

    # --- Main Drawing Form ---
    with st.form("circular_form"):
        st.header("Drawing Options")
        col1, col2 = st.columns(2)
        with col1:
            st.subheader("Basic Settings")
            st.markdown("##### Output Settings")
            c_prefix = st.text_input("Output prefix (optional):", help="Default is the basename of the input file name.")
            c_fmt = st.selectbox("Output format:", ["svg", "png", "pdf", "eps", "ps"], index=0, key="c_fmt", help="Output file format. Default is SVG, which is the fastest and most flexible for web display.")
            st.markdown("##### Definition Text")
            c_species = st.text_input("Species name (optional):", help='e.g., "<i>Escherichia coli</i>". Combining italic and block elements like "<i>Ca.</i> Tyloplasma litorale" cannot be reliably converted from SVG to PDF/PNG/EPS/PS. As a workaround, export to SVG format and convert to other formats using external tools like [Inkscape](https://inkscape.org/).')
            c_strain = st.text_input("Strain/isolate name (optional):", help='e.g., "K-12"')
            st.markdown("##### Legend")
            c_legend = st.selectbox("Legend:", ["right", "left", "upper_left", "upper_right", "lower_left", "lower_right", "none"], index=0, key="c_legend", help="Position of the legend. 'none' hides the legend.")
        with col2:
            st.subheader("Display Options")
            st.markdown("##### Track layout")
            c_track_type = st.selectbox("Track type:", ["tuckin", "middle", "spreadout"], index=0, key="c_track", help="Choose how features are displayed in the circular track. 'tuckin' is the default and most compact, 'middle' places features along the middle of the circle, and 'spreadout' spreads them around the circle. See [here](https://github.com/satoshikawato/gbdraw/blob/main/docs/TUTORIALS/1_Customizing_Plots.md#track-layout-style---track_type) for examples.")
            c_separate_strands = st.checkbox("Separate strands", value=True, key="c_strands", help="Display features on separate strands for better distinction of forward and reverse strands. See [here](https://github.com/satoshikawato/gbdraw/blob/main/docs/TUTORIALS/1_Customizing_Plots.md#strand-separation---separate_strands) for examples.")
            st.markdown("##### Label layout")
            c_show_labels = st.checkbox("Show labels", value=False, key="c_labels", help="Display feature labels on the circular map.")
            c_allow_inner_labels = st.checkbox("Allow inner labels", value=False, key="c_inner_labels", help="Enable inner labels as well as outer labels. This can help avoid label overlap in some casees but suppresses GC content and GC skew tracks.")
            st.markdown("##### Dinucleotide Tracks")        
            if c_allow_inner_labels:
                st.info("💡 Inner labels are enabled. GC content and GC skew tracks will be automatically suppressed to avoid overlap.")
                c_suppress_gc = True
                c_suppress_skew = True
            else:
                c_suppress_gc = st.checkbox("Suppress GC content track", value=False, key="c_gc_suppress", help="Suppress the GC content track.")
                c_suppress_skew = st.checkbox("Suppress GC skew track", value=False, key="c_skew_suppress", help="Suppress the GC skew track.")

        with st.expander("🔧 Advanced Options"):
            st.subheader("Advanced Drawing")
            adv_cols1, adv_cols2 = st.columns(2)
            with adv_cols1:
                st.markdown("##### Feature Selection")
                c_adv_feat = st.multiselect("Features (-k):", options=FEATURE_KEYS, default=["CDS","rRNA","tRNA","tmRNA","ncRNA","misc_RNA","repeat_region"], key="c_feat", help="Select which features to include in the circular map. Default includes CDS, tRNA, rRNA, and repeat regions.")
                st.markdown("##### Dinucleotide and Window/Step Size")
                c_adv_nt = st.text_input("Dinucleotide (--nt):", value="GC", key="c_nt", help="Dinucleotide to use for GC content and skew calculations. Default is 'GC'. See [here](https://github.com/satoshikawato/gbdraw/blob/main/docs/FAQ.md#q-can-i-plot-the-at-content-instead-of-gc-content) for details.")
                c_adv_win = st.number_input("Window size (-w):", key="c_win", step=1, help="Window size for GC content and skew calculations. Default: 1kb for genomes < 1Mb, 10kb for genomes <10Mb, 100kb for genomes >=10Mb. See [here](https://github.com/satoshikawato/gbdraw/blob/main/docs/FAQ.md#q-how-can-i-make-the-gc-content-graph-smootherfiner) for details.")
                c_adv_step = st.number_input("Step size (-s):", key="c_step", step=1, help="Step size for GC content and skew calculations. Default: 100 bp for genomes < 1Mb, 1kb for genomes <10Mb, 10kb for genomes >=10Mb. See [here](https://github.com/satoshikawato/gbdraw/blob/main/docs/FAQ.md#q-how-can-i-make-the-gc-content-graph-smootherfiner) for details.")
                st.markdown("##### Stroke Customization")
                c_adv_blk_color = st.color_picker("Block stroke color:", value="#808080", key="c_b_color", help="Color of the outline for feature blocks.")
                c_adv_blk_width = st.number_input("Block stroke width:", key="c_b_width", help="Block stroke width. Default: 2 pt for genomes <= 50 kb, 0 pt for genomes >= 50 kb.")
                c_adv_line_color = st.color_picker("Line stroke color:", value="#D3D3D3", key="c_l_color", help="Color of the lines representing introns.")
                c_adv_line_width = st.number_input("Line stroke width:", key="c_l_width", help="Width of the lines representing introns. Default: 5 pt for genomes <= 50 kb, 1 pt for genomes >= 50 kb.")
                st.markdown("##### Axis Customization")
                c_adv_axis_color = st.color_picker("Axis stroke color:", value="#808080", key="c_axis_color", help="Color of the main axis line.")
                c_adv_axis_width = st.number_input("Axis stroke width:", key="c_axis_width", help="Width of the main axis line. Default: 3 pt for genomes <= 50 kb, 1 pt for genomes >= 50 kb.")
            with adv_cols2:
                st.markdown("##### Scale Customization")
                c_adv_scale_interval = st.number_input("Scale interval (bp):", key="c_scale_interval", step=1, help="Manual scale interval (in bp). Overrides automatic calculation.")
                st.markdown("##### Font Sizes")
                c_adv_def_font_size = st.number_input("Definition font size (default: 18 pt):", value=18.0, min_value=1.0, step=0.5, key="c_def_font_size", help="Font size for the species and strain definition text. See [here](https://github.com/satoshikawato/gbdraw/blob/main/docs/TUTORIALS/3_Advanced_Customization.md#definition-font-size---definition_font_size) for details.")
                c_adv_label_font_size = st.number_input("Label font size (default: 8 pt (>=50 kb) or 16 pt (<50 kb):", key="c_label_font_size", help="Font size for feature labels. Default is 8 pt for genomes >= 50 kb, 16 pt for smaller genomes. See [here](https://github.com/satoshikawato/gbdraw/blob/main/docs/TUTORIALS/3_Advanced_Customization.md#label-font-size---label_font_size) for details.")
                st.markdown("##### Label Radius Offsets")
                st.markdown("Adjust the radius offsets for label placement. See [here](https://github.com/satoshikawato/gbdraw/blob/main/docs/TUTORIALS/3_Advanced_Customization.md#part-3-fine-tuning-plot-aesthetics) for details.")
                col_outer, col_inner = st.columns(2)
                with col_outer:
                    st.write("Outer Labels")
                    c_adv_outer_x_offset = st.number_input("X Radius Offset", value=1.0, key="c_outer_x_offset", min_value=0.5, max_value=2.0, step=0.1, help="Adjust the X radius offset for outer labels. Increasing this value moves labels further from the genome circle.")
                    c_adv_outer_y_offset = st.number_input("Y Radius Offset", value=1.0, key="c_outer_y_offset", min_value=0.5, max_value=2.0, step=0.1, help="Adjust the Y radius offset for outer labels. Increasing this value moves labels further from the genome circle.")
                with col_inner:
                    st.write("Inner Labels")
                    c_adv_inner_x_offset = st.number_input("X Radius Offset", value=1.0, key="c_inner_x_offset", min_value=0.5, max_value=2.0, step=0.1, help="Adjust the X radius offset for inner labels. Increasing this value moves labels further into the center.")
                    c_adv_inner_y_offset = st.number_input("Y Radius Offset", value=1.0, key="c_inner_y_offset", min_value=0.5, max_value=2.0, step=0.1, help="Adjust the Y radius offset for inner labels. Increasing this value moves labels further into the center.")

        c_submitted = st.form_submit_button("🚀 Run gbdraw Circular", type="primary")

    if "c_results" not in st.session_state:
        st.session_state.c_results = None
        
    if c_submitted:
        # Start execution in a secure temporary directory
        with tempfile.TemporaryDirectory() as tmpdirname:
            tmp_dir = Path(tmpdirname)
            
            circular_args = []
            selected_file_for_prefix = None

            if st.session_state.c_input_type == "GenBank":
                selected_gb_file = st.session_state.get("c_gb_manual", "")
                if not selected_gb_file:
                    st.error("Please select a GenBank file.")
                    st.stop()
                gb_path = st.session_state.uploaded_files[selected_gb_file]
                circular_args.extend(["--gbk", gb_path])
                selected_file_for_prefix = selected_gb_file
            else:  # GFF3 + FASTA
                selected_gff_file = st.session_state.get("c_gff_manual", "")
                selected_fasta_file = st.session_state.get("c_fasta_manual", "")
                if not selected_gff_file or not selected_fasta_file:
                    st.error("Please select both a GFF3 file and a FASTA file.")
                    st.stop()
                gff_path = st.session_state.uploaded_files[selected_gff_file]
                fasta_path = st.session_state.uploaded_files[selected_fasta_file]
                circular_args.extend(["--gff", gff_path, "--fasta", fasta_path])
                selected_file_for_prefix = selected_gff_file

            if not selected_file_for_prefix:
                st.error("Please select an input file.")
                st.stop()

            sanitized_prefix = sanitize_filename(c_prefix)
            prefix = sanitized_prefix or Path(selected_file_for_prefix).stem
            
            # Specify the output destination within the temporary directory
            full_output_prefix = tmp_dir / prefix
            circular_args.extend(["-o", str(full_output_prefix), "-f", c_fmt, "--track_type", c_track_type])
            
            if c_species:
                circular_args.extend(["--species", c_species])
            if c_strain:
                circular_args.extend(["--strain", c_strain])
            if c_show_labels: circular_args.append("--show_labels")
            if c_separate_strands: circular_args.append("--separate_strands")
            if c_allow_inner_labels:
                circular_args.extend(["--allow_inner_labels", "--suppress_gc", "--suppress_skew"])
            else:
                if c_suppress_gc: circular_args.append("--suppress_gc")
                if c_suppress_skew: circular_args.append("--suppress_skew")

            if c_adv_def_font_size:
                circular_args += ["--definition_font_size", str(c_adv_def_font_size)]
            if c_adv_label_font_size:
                circular_args += ["--label_font_size", str(c_adv_label_font_size)]

            if c_legend != "right": circular_args += ["-l", c_legend]
            
            selected_palette = st.session_state.get("c_palette_selector")
            if selected_palette: circular_args += ["--palette", selected_palette]

            if 'custom_circular_colors' in st.session_state:
                custom_color_filename = f"custom_colors_c_{uuid.uuid4().hex[:8]}.tsv"
                save_path = UPLOAD_DIR / custom_color_filename
                with open(save_path, "w") as f:
                    for feature, color in st.session_state.custom_circular_colors.items():
                        f.write(f"{feature}\t{color}\n")
                circular_args += ["-d", str(save_path)]
            circular_args += ["-k", ",".join(c_adv_feat), "-n", c_adv_nt]
            if c_adv_win: circular_args += ["--window", str(c_adv_win)]
            if c_adv_step: circular_args += ["--step", str(c_adv_step)]
            if c_adv_scale_interval: circular_args += ["--scale_interval", str(c_adv_scale_interval)]
            circular_args += ["--block_stroke_color", c_adv_blk_color]
            if c_adv_blk_width: circular_args += ["--block_stroke_width", str(c_adv_blk_width)]
            circular_args += ["--line_stroke_color", c_adv_line_color]
            if c_adv_line_width: circular_args += ["--line_stroke_width", str(c_adv_line_width)]
            circular_args += ["--axis_stroke_color", c_adv_axis_color]
            if c_adv_axis_width: circular_args += ["--axis_stroke_width", str(c_adv_axis_width)]
            circular_args += ["--outer_label_x_radius_offset", str(c_adv_outer_x_offset)]
            circular_args += ["--outer_label_y_radius_offset", str(c_adv_outer_y_offset)]
            circular_args += ["--inner_label_x_radius_offset", str(c_adv_inner_x_offset)]
            circular_args += ["--inner_label_y_radius_offset", str(c_adv_inner_y_offset)]

            # --- QUALIFIER PRIORITY & BLACKLIST LOGIC ---
            selected_prio_file = st.session_state.get("c_qualifier_priority_file_manual", "")
            if selected_prio_file:
                prio_path = st.session_state.uploaded_files[selected_prio_file]
                circular_args += ["--qualifier_priority", prio_path]
            else:
                prio_lines = []
                for row_state in st.session_state.manual_priorities:
                    feature_key = f"feature_{row_state['id']}"
                    qualifiers_key = f"qualifiers_{row_state['id']}"
                    feature_value = st.session_state.get(feature_key, "")
                    qualifiers_value = st.session_state.get(qualifiers_key, "")
                    if feature_value and qualifiers_value:
                        prio_lines.append(f"{feature_value}\t{qualifiers_value}")
                
                if prio_lines:
                    prio_content = "\n".join(prio_lines)
                    save_path = UPLOAD_DIR / f"qual_prio_c_{uuid.uuid4().hex}"
                    with open(save_path, "w", encoding="utf-8") as f:
                        f.write(prio_content)
                    circular_args += ["--qualifier_priority", str(save_path)]

            if c_filter_mode == "Blacklist (exclude keywords)":
                selected_blacklist_file = st.session_state.get("c_blacklist_file_manual", "")
                if selected_blacklist_file:
                    blacklist_path = st.session_state.uploaded_files[selected_blacklist_file]
                    circular_args += ["--label_blacklist", blacklist_path]
                else:
                    blacklist_keywords = st.session_state.get("c_blacklist_manual", "")
                    if blacklist_keywords:
                        circular_args += ["--label_blacklist", blacklist_keywords]
            elif c_filter_mode == "Whitelist (include keywords)":
                selected_whitelist_file = st.session_state.get("c_whitelist_file_manual", "")
                if selected_whitelist_file:
                    whitelist_path = st.session_state.uploaded_files[selected_whitelist_file]
                    circular_args += ["--label_whitelist", whitelist_path]
                whitelist_lines = []
                for row in st.session_state.manual_whitelist:
                    if row['feature'] and row['qualifier'] and row['keyword']:
                        whitelist_lines.append(f"{row['feature']}\t{row['qualifier']}\t{row['keyword']}")
                
                if whitelist_lines:
                    whitelist_content = "\n".join(whitelist_lines)
                    save_path = UPLOAD_DIR / f"label_whitelist_c_{uuid.uuid4().hex}.tsv"
                    with open(save_path, "w", encoding="utf-8") as f:
                        f.write(whitelist_content)
                    circular_args += ["--label_whitelist", str(save_path)]
            selected_d_color_file = st.session_state.get("c_d_color_manual", "")
            if selected_d_color_file: circular_args += ["-d", st.session_state.uploaded_files[selected_d_color_file]]
            selected_t_color_file = st.session_state.get("c_t_color_manual", "")
            if selected_t_color_file: circular_args += ["-t", st.session_state.uploaded_files[selected_t_color_file]]

            # --- Run gbdraw and handle output ---
            logger = logging.getLogger() 
            log_capture = io.StringIO()  
            command_str = f"gbdraw circular {' '.join(circular_args)}"
            log_capture.write(f"--- Executed Command ---\n{command_str}\n------------------------\n\n")
            stream_handler = logging.StreamHandler(log_capture)
            stream_handler.setLevel(logging.INFO) 
            logger.addHandler(stream_handler)
            
            start_time = time.time()
            had_exception = False
            exit_code = 0

            with st.spinner(f"Running: `{command_str}`"):
                try:
                    with redirect_stdout(log_capture), redirect_stderr(log_capture):
                        circular_main(circular_args)
                except SystemExit as e:
                    exit_code = e.code if e.code is not None else 1
                except Exception as e:
                    had_exception = True
                    log_capture.write(f"\n--- Streamlit App Exception ---\n{e}\n----------------------\n")
            
            end_time = time.time()
            duration = end_time - start_time
            log_capture.write(f"\n--- Execution Time ---\nTotal time: {duration:.2f} seconds\n----------------------")
            log_content = log_capture.getvalue()
            logger.removeHandler(stream_handler)
            
            # Find output files in the temporary directory
            output_files = sorted(list(tmp_dir.glob(f"{prefix}*.{c_fmt}")))
            
            is_successful = (
                not had_exception and
                exit_code == 0 and
                "ERROR:" not in log_content.upper() and
                bool(output_files)
            )

            if is_successful:
                results_data = []
                for out_path in output_files:
                    file_data = out_path.read_bytes()
                    results_data.append({
                        "name": out_path.name,
                        "data": file_data,
                        "suffix": out_path.suffix.lower()
                    })
                
                st.session_state.c_results = {
                    "success": True,
                    "files": results_data,
                    "log": log_content,
                    "fmt": c_fmt
                }
            else:
                st.session_state.c_results = {
                    "success": False,
                    "log": log_content
                }

    if st.session_state.c_results:
        res = st.session_state.c_results
        
        if res["success"]:
            st.success("✅ gbdraw finished successfully.")
            st.subheader("🌀 Circular Drawing Output")
            
            for f_info in res["files"]:
                file_name = f_info["name"]
                file_data = f_info["data"]
                file_extension = f_info["suffix"]
                fmt = res["fmt"]

                if file_extension == ".svg":
                    st.image(file_data.decode("utf-8"), caption=file_name)
                elif file_extension == ".png":
                    st.image(file_data, caption=file_name)
                else:
                    st.info(f"📄 Preview is not available for {file_extension.upper()} format.")
                
                st.download_button(
                    f"⬇️ Download {file_name}",
                    data=file_data,
                    file_name=file_name,
                    key=f"download_{file_name}",
                    mime=f"image/{fmt}" if fmt in ["png", "svg"] else "application/octet-stream"
                )
                st.markdown("---")
            
            with st.expander("Show Log"):
                st.text(res["log"])
        else:
            st.error("gbdraw execution failed. Please check the log for details.")
            with st.expander("Show Log", expanded=True):
                st.text(res["log"])

# --- LINEAR MODE ---
if selected_mode == "📏 Linear":
    st.header("Linear Mode")
    st.markdown("---")
    st.subheader("Input Genome Files")

    l_input_type = st.radio(
        "Input file type for all sequences",
        ("GenBank", "GFF3 + FASTA"),
        key="l_input_type",
        horizontal=True,
    )

    input_container = st.container()
    with input_container:
        for i in range(st.session_state.linear_seq_count):
            st.markdown(f"#### Sequence {i+1}")
            if l_input_type == "GenBank":
                cols = st.columns([3, 3])
                with cols[0]:
                    create_manual_selectbox("GenBank File", file_options, f"l_gb_{i}")
                if i < st.session_state.linear_seq_count - 1:
                    with cols[1]:
                        create_manual_selectbox(f"Comparison File {i+1} (BLAST)", file_options, f"l_blast_{i}")
            else: # GFF3 + FASTA
                cols = st.columns([2, 2, 2])
                with cols[0]:
                    create_manual_selectbox("GFF3 File", file_options, f"l_gff_{i}")
                with cols[1]:
                    create_manual_selectbox("FASTA File", file_options, f"l_fasta_{i}")
                if i < st.session_state.linear_seq_count - 1:
                    with cols[2]:
                        create_manual_selectbox(f"Comparison File {i+1} (BLAST)", file_options, f"l_blast_{i}")


    b_col1, b_col2, _ = st.columns([1.5, 2, 5])
    if b_col1.button("➕ Add Sequence"):
        st.session_state.linear_seq_count += 1
        st.rerun()
    if b_col2.button("➖ Remove Last Sequence") and st.session_state.linear_seq_count > 1:
        st.session_state.linear_seq_count -= 1
        # Clean up keys for the removed sequence
        keys_to_remove = [
            f"l_gb_{st.session_state.linear_seq_count}",
            f"l_gff_{st.session_state.linear_seq_count}",
            f"l_fasta_{st.session_state.linear_seq_count}",
            f"l_blast_{st.session_state.linear_seq_count -1}"
        ]
        for key in keys_to_remove:
            manual_key = f"{key}_manual"
            if manual_key in st.session_state: del st.session_state[manual_key]
            if key in st.session_state: del st.session_state[key]
        st.rerun()
    st.markdown("---")
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("🎨 Color Customization")
        create_manual_selectbox("Default-override color file (optional):", file_options, "l_d_color", help="Tab-separated value (TSV) file that overrides the color palette. Feature types not specified will use the selected palette colors. See [here](https://github.com/satoshikawato/gbdraw/blob/main/docs/TUTORIALS/3_Advanced_Customization.md#method-1-override-default-colors--d) for details.")
        create_manual_selectbox("Feature-specific color file (optional):", file_options, "l_t_color", help="TSV file that overrides the color palette for specific features. Features not specified will use the selected palette colors or the default-override colors if provided. See [here](https://github.com/satoshikawato/gbdraw/blob/main/docs/TUTORIALS/3_Advanced_Customization.md#method-2-feature-specific-colors--t) for details.")

        def linear_palette_changed():
            new_palette = st.session_state.l_palette_selector
            st.session_state.custom_linear_colors = get_palette_colors(new_palette).copy()
            st.session_state.color_picker_seed += 1

        if 'custom_linear_colors' not in st.session_state:
            default_palette = PALETTES[0] if PALETTES and PALETTES[0] != "" else "default"
            st.session_state.custom_linear_colors = get_palette_colors(default_palette).copy()

        st.selectbox(
            "Base color palette:",
            PALETTES, index=11,
            key="l_palette_selector",
            on_change=linear_palette_changed,
            help="Select a base color palette. You can further customize individual feature colors below. See [here](https://github.com/satoshikawato/gbdraw/blob/main/examples/color_palette_examples.md) for palette examples.",
        )

        st.write("Click on the color boxes below to customize the default colors for each feature.")
        l_cols = st.columns(5)
        
        l_color_keys = sorted(st.session_state.custom_linear_colors.keys())
        for i, feature in enumerate(l_color_keys):
            col = l_cols[i % 5]
            with col:
                st.session_state.custom_linear_colors[feature] = st.color_picker(
                    label=feature,
                    value=st.session_state.custom_linear_colors[feature],
                    key=f"l_color_picker_{feature}_{st.session_state.color_picker_seed}"
                )
    
    # --- Qualifier Priority & Label Filtering Section (OUTSIDE the form) ---
    with col2:
        st.subheader("✏️ Label Content Customization")
        
        # Qualifier Priority
        st.markdown("##### Qualifier Priority")
        st.markdown("Select a TSV file from the sidebar OR define priorities manually below. If a file is selected, manual entries are ignored.")
        create_manual_selectbox(
            "Qualifier Priority File (optional)",
            file_options,
            "l_qualifier_priority_file",
            help="Tab-separated value (TSV) file defining qualifier priorities for feature labels. See [here](https://github.com/satoshikawato/gbdraw/blob/main/docs/TUTORIALS/3_Advanced_Customization.md#changing-label-content---qualifier_priority) for details."
        )
        with st.expander("Or, define priorities manually:", expanded=False):
            col1, col2, _ = st.columns([3, 5, 1])
            col1.markdown("**Feature Type**")
            col2.markdown("**Qualifier Priority (comma-separated)**")

            for row in st.session_state.manual_priorities:
                col1, col2, col3 = st.columns([3, 5, 1])
                col1.text_input("Feature", value=row['feature'], key=f"feature_{row['id']}", label_visibility="collapsed")
                col2.text_input("Qualifiers", value=row['qualifiers'], key=f"qualifiers_{row['id']}", label_visibility="collapsed")
                col3.button("➖", key=f"remove_{row['id']}", on_click=remove_priority_row, args=(row['id'],))

            st.button("➕ Add Row", on_click=add_priority_row, use_container_width=True)

        # Label Content Filtering
        st.markdown("##### Label Content Filtering")
        l_filter_mode = st.radio(
            "Select label filtering mode:",
            ("None", "Blacklist (exclude keywords)", "Whitelist (include keywords)"),
            key="l_filter_mode",
            horizontal=True, help="Choose how to filter feature labels based on keywords. See [here](https://github.com/satoshikawato/gbdraw/blob/main/docs/TUTORIALS/3_Advanced_Customization.md#part-2-advanced-label-control) for details.")


        if l_filter_mode == "Blacklist (exclude keywords)":
            st.info("Select a file with blacklist keywords (one per line) OR enter them manually below. If a file is selected, manual entry is ignored. See [here](https://github.com/satoshikawato/gbdraw/blob/main/docs/TUTORIALS/3_Advanced_Customization.md#part-2-advanced-label-control) for details.")
            create_manual_selectbox(
                "Blacklist File (optional)",
                file_options,
                "l_blacklist_file"
            )
            st.text_area(
                "Or, enter keywords manually (comma-separated):",
                value="hypothetical, uncharacterized, putative, unknown",
                help="Features with these keywords in their labels will be hidden.",
                key="l_blacklist_manual"
            )
        elif l_filter_mode == "Whitelist (include keywords)":
            st.info("Select a file with label whitelist (one per line) OR enter them manually below. If a file is selected, manual entry is ignored. See [here](https://github.com/satoshikawato/gbdraw/blob/main/docs/TUTORIALS/3_Advanced_Customization.md#whitelisting-labels---label_whitelist) for details.")  
            create_manual_selectbox(
                "Whitelist File (optional)",
                file_options,
                "l_whitelist_file"
            )
            st.info("Define rules to ONLY show labels containing specific keywords. For example, show 'CDS' features where the 'product' contains 'DNA polymerase'.")
            with st.container():
                wl_col1, wl_col2, wl_col3, _ = st.columns([3, 3, 4, 1])
                wl_col1.markdown("**Feature Type**")
                wl_col2.markdown("**Qualifier**")
                wl_col3.markdown("**Keyword to Include**")
                for row in st.session_state.manual_whitelist:
                    wl_col1, wl_col2, wl_col3, wl_col4 = st.columns([3, 3, 4, 1])
                    row['feature'] = wl_col1.selectbox("Feature", options=FEATURE_KEYS, index=FEATURE_KEYS.index(row['feature']) if row['feature'] in FEATURE_KEYS else 0, key=f"l_wl_feature_{row['id']}", label_visibility="collapsed")
                    row['qualifier'] = wl_col2.selectbox("Qualifier", options=QUALIFIER_KEYS, index=QUALIFIER_KEYS.index(row['qualifier']) if row['qualifier'] in QUALIFIER_KEYS else 0, key=f"l_wl_qualifier_{row['id']}", label_visibility="collapsed")
                    row['keyword'] = wl_col3.text_input("Keyword", value=row['keyword'], key=f"l_wl_keyword_{row['id']}", label_visibility="collapsed")
                    wl_col4.button("➖", key=f"l_wl_remove_{row['id']}", on_click=remove_whitelist_row, args=(row['id'],))

            st.button("➕ Add Whitelist Row", on_click=add_whitelist_row, use_container_width=True, key="l_add_wl")
    st.markdown("---")

    with st.form("linear_form"):
        st.header("Drawing Options")
        col1, col2 = st.columns(2)
        with col1:
            st.subheader("Basic Settings")
            st.markdown("##### Output Settings")
            l_prefix = st.text_input("Output prefix:", value="linear", key="l_prefix", help="Prefix for output files. Default is 'linear'.")
            l_fmt = st.selectbox("Output format:", ["svg", "png", "pdf", "eps", "ps"], index=0, key="l_fmt", help="Output file format. Default is SVG, which is the fastest and most flexible for web display.")
            st.markdown("##### Legend position")
            l_legend = st.selectbox("Legend:", ["right", "left", "top", "bottom", "none"], index=0, key="l_legend", help="Position of the legend. 'none' hides the legend.")
            st.markdown("##### Scale layout")
            l_scale_style = st.selectbox("Scale style (--scale_style):", ["bar", "ruler"], index=0, key="l_scale_style", help="Style of the scale bar on the linear map. 'bar' draws a simple line, 'ruler' draws a ruler with ticks and labels.")
        with col2:
            st.subheader("Display Options")
            st.markdown("##### Track layout")
            l_separate_strands = st.checkbox("Separate strands", value=True, key="l_strands", help="Display features on separate strands for better distinction of forward and reverse strands.")
            l_normalize_length = st.checkbox("Normalize sequence lengths", value=False, key="l_normalize", help="Normalize the lengths of all sequences to be equal.The length bar will be suppressed when this option is enabled.")
            l_align_center = st.checkbox("Align center", value=False, key="l_align", help="Align the linear map to the center of the page. This can help with aesthetics, especially for long sequences.")
            l_resolve_overlaps = st.checkbox("Resolve overlaps (experimental)", value=False, key="l_overlaps", help="Attempt to resolve label overlaps. This is experimental and may not work well for all genomes.")
            st.markdown("##### Label layout")
            l_show_labels_mode = st.radio(
                "Show labels", 
                options=["None", "All", "First"], 
                index=0, 
                key="l_labels_mode", 
                horizontal=True,
                help="Display feature labels on the linear map. 'First' shows labels only for the top track."
            )
            st.markdown("##### Dinucleotide Tracks")     
            l_show_gc = st.checkbox("Show GC content", value=False, key="l_gc", help="Display the GC content track on the linear map.")
            l_show_skew = st.checkbox("Show GC skew", value=False, key="l_skew", help="Display the GC skew track on the linear map.")
            

        with st.expander("🔧 Advanced Options"):

            st.subheader("Advanced Drawing")
            adv_cols1, adv_cols2 = st.columns(2)
            with adv_cols1:
                st.markdown("##### Feature Selection")
                l_adv_feat = st.multiselect("Features (-k):", options=FEATURE_KEYS, default=["CDS","rRNA","tRNA","tmRNA","ncRNA","misc_RNA","repeat_region"], key="l_feat", help="Select which features to include in the linear map. Default includes CDS, tRNA, rRNA, and repeat regions.")
                l_adv_feat_height = st.number_input("Feature height:", key="l_feat_height", help="Height of the feature blocks in pixels. Default: 20 pixels for genomes >= 50 kb, 80 pixels for smaller genomes.")
                st.markdown("##### Dinucleotide and Window/Step Size")
                l_adv_gc_height = st.number_input("GC content and skew height:", key="l_gc_height", help="Height of the GC content track in pixels. Default: 20 pixels.")
                l_adv_nt = st.text_input("nt (--nt):", value="GC", key="l_nt", help="Dinucleotide to use for GC content and skew calculations. Default is 'GC'.")
                l_adv_win = st.number_input("Window size (-w):", key="l_win", step=1, help="Window size for GC content and skew calculations. Default: 1kb for genomes < 1Mb, 10kb for genomes <10Mb, 100kb for genomes >=10Mb")
                l_adv_step = st.number_input("Step size (-s):", key="l_step", step=1, help="Step size for GC content and skew calculations. Default: 100 bp for genomes < 1Mb, 1kb for genomes <10Mb, 10kb for genomes >=10Mb")
                st.markdown("##### Legend")
                l_adv_legend_font_size = st.number_input("Legend font size:", key="l_legend_font_size", help="Font size for the legend. Default: 20 (pixels, 96 dpi) for genomes <= 50 kb, 24 for genomes >= 50 kb")
                l_adv_legend_box_size = st.number_input("Legend box size:", key="l_legend_box_size", help="Box size for the legend. Default: 16 (pt) for genomes <= 50 kb, 20 for genomes >= 50 kb)")
                st.markdown("##### Font Sizes")
                l_adv_def_font_size = st.number_input("Definition font size:", key="l_def_font_size", help="Font size for the definition text beside each sequence. Default: 24 pt for genomes <= 50 kb, 10 pt for genomes >= 50 kb")
                l_adv_label_font_size = st.number_input("Label font size:", key="l_label_font_size", help="Font size for feature labels. Default: 24 pt for genomes <= 50 kb, 5 pt for genomes >= 50 kb")

            with adv_cols2:
                st.markdown("##### Scale Customization")
                l_adv_scale_interval = st.number_input("Scale interval (bp):", key="l_scale_interval", step=1, help="Manual scale interval (in bp). Overrides automatic calculation.")
                l_adv_scale_font_size = st.number_input("Scale font size (--scale_font_size):", key="l_scale_font_size", help="Font size for the scale labels on the linear map. Default: 24 (pt) for genomes <= 50 kb, 16 for genomes >= 50 kb")
                l_adv_scale_stroke_color = st.color_picker("Scale stroke color:", value="#000000", key="l_scale_color", help="Color of the scale bar.")
                l_adv_scale_stroke_width = st.number_input("Scale stroke width:", key="l_scale_width", help="Width of the scale bar. Default: 3 pt.")
                st.markdown("##### Stroke Customization")
                l_adv_blk_color = st.color_picker("Block stroke color:", value="#808080", key="l_b_color", help="Color of the outline for feature blocks.")
                l_adv_blk_width = st.number_input("Block stroke width:", key="l_b_width", help="Width of the outline for feature blocks. Default: 2 pt for genomes <= 50 kb, 0 pt for genomes >= 50 kb")
                l_adv_line_color = st.color_picker("Line stroke color:", value="#D3D3D3", key="l_l_color", help="Color of the lines representing introns.")
                l_adv_line_width = st.number_input("Line stroke width:", key="l_l_width", help="Width of the lines representing introns. Default: 5 pt for genomes <= 50 kb, 1 pt for genomes >= 50 kb")
                st.markdown("##### Axis Customization")
                l_adv_axis_color = st.color_picker("Axis stroke color:", value="#808080", key="l_axis_color", help="Color of the main axis line.")
                l_adv_axis_width = st.number_input("Axis stroke width:", key="l_axis_width", help="Width of the main axis line. Default: 5 pt for genomes <= 50 kb, 2 pt for genomes >= 50 kb")
            st.subheader("Comparison Filters")
            l_adv_comp_height = st.number_input("Comparison height:", key="l_comp_height", help="Height of the comparison blocks in pixels. Default: 60 pixels.")
            l_adv_bitscore = st.number_input("Min bitscore:", value=50.0, key="l_bitscore", help="Minimum bitscore for BLAST comparisons. Default is 50.0.")
            l_adv_evalue = st.text_input("Max E-value:", value="1e-2", key="l_evalue", help="Maximum E-value for BLAST comparisons. Default is '1e-2'.")
            l_adv_identity = st.number_input("Min identity (%):", value=0.0, key="l_identity", help="Minimum identity percentage for BLAST comparisons. Default is 0.0%.")
            
        l_submitted = st.form_submit_button("🚀 Run gbdraw Linear", type="primary")

    if "l_results" not in st.session_state:
        st.session_state.l_results = None
        
    if l_submitted:
        # Use a secure temporary directory for execution
        with tempfile.TemporaryDirectory() as tmpdirname:
            tmp_dir = Path(tmpdirname)
            
            linear_args = []
            num_sequences = 0

            if st.session_state.l_input_type == "GenBank":
                selected_gb = [st.session_state.get(f"l_gb_{i}_manual", "") for i in range(st.session_state.linear_seq_count)]
                selected_gb = [f for f in selected_gb if f]
                if not selected_gb:
                    st.error("Please select at least one GenBank file.")
                    st.stop()
                gb_paths = [st.session_state.uploaded_files[f] for f in selected_gb]
                linear_args.extend(["--gbk", *gb_paths])
                num_sequences = len(gb_paths)
            else: # GFF3 + FASTA
                gff_paths = []
                fasta_paths = []
                for i in range(st.session_state.linear_seq_count):
                    gff_file = st.session_state.get(f"l_gff_{i}_manual", "")
                    fasta_file = st.session_state.get(f"l_fasta_{i}_manual", "")
                    if gff_file and fasta_file:
                        gff_paths.append(st.session_state.uploaded_files[gff_file])
                        fasta_paths.append(st.session_state.uploaded_files[fasta_file])
                    elif gff_file or fasta_file:
                        st.error(f"For Sequence {i+1}, both GFF3 and FASTA files must be provided as a pair.")
                        st.stop()
                
                if not gff_paths:
                    st.error("Please select at least one pair of GFF3 and FASTA files.")
                    st.stop()
                
                linear_args.extend(["--gff", *gff_paths, "--fasta", *fasta_paths])
                num_sequences = len(gff_paths)
            
            selected_blast = [st.session_state.get(f"l_blast_{i}_manual", "") for i in range(st.session_state.linear_seq_count - 1)]
            selected_blast = [f for f in selected_blast if f]
            
            if selected_blast and len(selected_blast) != num_sequences - 1:
                st.error(f"Please provide {num_sequences - 1} comparison file(s) for {num_sequences} sequence files.")
                st.stop()

            sanitized_prefix = sanitize_filename(l_prefix) 
            prefix = sanitized_prefix or "linear"
            
            # Specify output path inside the temporary directory
            full_output_prefix = tmp_dir / prefix
            output_path = full_output_prefix.with_suffix(f".{l_fmt}")
            
            linear_args.extend(["-o", str(full_output_prefix), "-f", l_fmt])

            if selected_blast:
                blast_paths = [st.session_state.uploaded_files[f] for f in selected_blast]
                linear_args += ["-b", *blast_paths]
            if l_show_labels_mode == "All":
                linear_args.extend(["--show_labels", "all"])
            elif l_show_labels_mode == "First":
                linear_args.extend(["--show_labels", "first"])
            if l_separate_strands: linear_args.append("--separate_strands")
            if l_align_center: linear_args.append("--align_center")
            if l_show_gc: linear_args.append("--show_gc")
            if l_show_skew: linear_args.append("--show_skew")
            if l_resolve_overlaps: linear_args.append("--resolve_overlaps")
            if l_legend != "right": linear_args += ["-l", l_legend]
            if l_normalize_length: linear_args.append("--normalize_length")
            if l_adv_gc_height:
                linear_args += ["--gc_height", str(l_adv_gc_height)]
            if l_adv_def_font_size:
                linear_args += ["--definition_font_size", str(l_adv_def_font_size)]
            if l_adv_label_font_size:
                linear_args += ["--label_font_size", str(l_adv_label_font_size)]
            if l_adv_feat_height:
                linear_args += ["--feature_height", str(l_adv_feat_height)]
            if l_adv_comp_height:
                linear_args += ["--comparison_height", str(l_adv_comp_height)]
            if l_scale_style != "bar":
                linear_args += ["--scale_style", l_scale_style]
            if l_adv_scale_font_size:
                linear_args += ["--scale_font_size", str(l_adv_scale_font_size)]
            if l_adv_scale_interval:
                linear_args += ["--scale_interval", str(l_adv_scale_interval)]
            if l_adv_scale_stroke_width:
                linear_args += ["--scale_stroke_width", str(l_adv_scale_stroke_width)]
            if l_adv_scale_stroke_color:
                linear_args += ["--scale_stroke_color", str(l_adv_scale_stroke_color)]
            if l_adv_legend_box_size:
                linear_args += ["--legend_box_size", str(l_adv_legend_box_size)]
            if l_adv_legend_font_size:
                linear_args += ["--legend_font_size", str(l_adv_legend_font_size)]
            if l_adv_axis_color:
                linear_args += ["--axis_stroke_color", l_adv_axis_color]
            if l_adv_axis_width:
                linear_args += ["--axis_stroke_width", str(l_adv_axis_width)]
            if l_adv_blk_width:
                linear_args += ["--block_stroke_width", str(l_adv_blk_width)]
            selected_palette = st.session_state.get("l_palette_selector")
            if selected_palette: linear_args += ["--palette", selected_palette]

            if 'custom_linear_colors' in st.session_state:
                custom_color_filename = f"custom_colors_l_{uuid.uuid4().hex[:8]}.tsv"
                save_path = UPLOAD_DIR / custom_color_filename
                with open(save_path, "w") as f:
                    for feature, color in st.session_state.custom_linear_colors.items():
                        f.write(f"{feature}\t{color}\n")
                linear_args += ["-d", str(save_path)]
            
            linear_args += ["-k", ",".join(l_adv_feat), "-n", l_adv_nt]
            if l_adv_win: linear_args += ["--window", str(l_adv_win)]
            if l_adv_step: linear_args += ["--step", str(l_adv_step)]
            
            linear_args += ["--bitscore", str(l_adv_bitscore), "--evalue", l_adv_evalue, "--identity", str(l_adv_identity)]
            linear_args += ["--block_stroke_color", l_adv_blk_color]
            
            linear_args += ["--line_stroke_color", l_adv_line_color]
            if l_adv_line_width:
                linear_args += ["--line_stroke_width", str(l_adv_line_width)]

            
            # --- QUALIFIER PRIORITY & BLACKLIST LOGIC ---
            selected_prio_file = st.session_state.get("l_qualifier_priority_file_manual", "")
            if selected_prio_file:
                prio_path = st.session_state.uploaded_files[selected_prio_file]
                linear_args += ["--qualifier_priority", prio_path]
            else:
                prio_lines = []
                for row_state in st.session_state.manual_priorities:
                    feature_key = f"feature_{row_state['id']}"
                    qualifiers_key = f"qualifiers_{row_state['id']}"
                    feature_value = st.session_state.get(feature_key, "")
                    qualifiers_value = st.session_state.get(qualifiers_key, "")
                    if feature_value and qualifiers_value:
                        prio_lines.append(f"{feature_value}\t{qualifiers_value}")
                
                if prio_lines:
                    prio_content = "\n".join(prio_lines)
                    save_path = UPLOAD_DIR / f"qual_prio_l_{uuid.uuid4().hex}"
                    with open(save_path, "w", encoding="utf-8") as f:
                        f.write(prio_content)
                    linear_args += ["--qualifier_priority", str(save_path)]
            
            if l_filter_mode == "Blacklist (exclude keywords)":
                selected_blacklist_file = st.session_state.get("l_blacklist_file_manual", "")
                if selected_blacklist_file:
                    blacklist_path = st.session_state.uploaded_files[selected_blacklist_file]
                    linear_args += ["--label_blacklist", blacklist_path]
                else:
                    blacklist_keywords = st.session_state.get("l_blacklist_manual", "")
                    if blacklist_keywords:
                        linear_args += ["--label_blacklist", blacklist_keywords]
            elif l_filter_mode == "Whitelist (include keywords)":
                selected_whitelist_file = st.session_state.get("l_whitelist_file_manual", "")
                if selected_whitelist_file:
                    whitelist_path = st.session_state.uploaded_files[selected_whitelist_file]
                    linear_args += ["--label_whitelist", whitelist_path]
                whitelist_lines = []
                for row in st.session_state.manual_whitelist:
                    if row['feature'] and row['qualifier'] and row['keyword']:
                        whitelist_lines.append(f"{row['feature']}\t{row['qualifier']}\t{row['keyword']}")
                if whitelist_lines:
                    whitelist_content = "\n".join(whitelist_lines)
                    save_path = UPLOAD_DIR / f"label_whitelist_l_{uuid.uuid4().hex}.tsv"
                    with open(save_path, "w", encoding="utf-8") as f:
                        f.write(whitelist_content)
                    linear_args += ["--label_whitelist", str(save_path)]
            selected_d_color_file = st.session_state.get("l_d_color_manual", "")
            if selected_d_color_file: linear_args += ["-d", st.session_state.uploaded_files[selected_d_color_file]]
            selected_t_color_file = st.session_state.get("l_t_color_manual", "")
            if selected_t_color_file: linear_args += ["-t", st.session_state.uploaded_files[selected_t_color_file]]
            
            logger = logging.getLogger()
            log_capture = io.StringIO()
            command_str = f"gbdraw linear {' '.join(map(str, linear_args))}"
            stream_handler = logging.StreamHandler(log_capture)
            stream_handler.setLevel(logging.INFO)
            logger.addHandler(stream_handler)
            log_capture.write(f"--- Executed Command ---\n{command_str}\n------------------------\n\n")
            
            start_time = time.time()
            had_exception = False
            exit_code = 0
            
            with st.spinner(f"Running: `{command_str}`"):
                try:
                    with redirect_stdout(log_capture), redirect_stderr(log_capture):
                        linear_main(linear_args)
                except SystemExit as e:
                    exit_code = e.code if e.code is not None else 1
                except Exception as e:
                    had_exception = True
                    log_capture.write(f"\n--- Streamlit App Exception ---\n{e}\n----------------------\n")

            end_time = time.time()
            duration = end_time - start_time
            log_capture.write(f"\n--- Execution Time ---\nTotal time: {duration:.2f} seconds\n----------------------")
            log_content = log_capture.getvalue()
            logger.removeHandler(stream_handler)

            is_successful = (
                not had_exception and
                exit_code == 0 and
                "ERROR:" not in log_content.upper() and
                output_path.exists()
            )

            # Display results (load file into memory for download)
            if is_successful:
                file_data = output_path.read_bytes()
                
                st.session_state.l_results = {
                    "success": True,
                    "file_name": output_path.name,
                    "data": file_data,
                    "suffix": output_path.suffix.lower(),
                    "fmt": l_fmt,
                    "log": log_content
                }
            else:
                st.session_state.l_results = {
                    "success": False,
                    "log": log_content
                }

    if st.session_state.l_results:
        res = st.session_state.l_results

        if res["success"]:
            st.success("✅ gbdraw finished successfully.")
            st.subheader("📏 Linear Drawing Output")

            file_name = res["file_name"]
            file_data = res["data"]
            file_extension = res["suffix"]
            fmt = res["fmt"]
            
            if file_extension == ".svg":
                st.image(file_data.decode("utf-8"), caption=file_name)
            elif file_extension == ".png":
                st.image(file_data, caption=file_name)
            else:
                st.info(f"📄 Preview is not available for {file_extension.upper()} format. Please use the download button below.")
            
            st.download_button(
                f"⬇️ Download {file_name}",
                data=file_data,
                file_name=file_name,
                mime=f"image/{fmt}" if fmt in ["png", "svg"] else "application/octet-stream"
            )
            
            with st.expander("Show Log"):
                st.text(res["log"])
        else:
            st.error("gbdraw execution failed. Please check the log for details.")
            with st.expander("Show Log", expanded=True):
                st.text(res["log"])

# --- Footer ---
st.markdown("---")
st.markdown(
    f"""
    Author: [Satoshi Kawato](https://github.com/satoshikawato)  |
    Source: [gbdraw](https://github.com/satoshikawato/gbdraw)  |
    Version: {VERSION}-{COMMIT_ID}
    """,
    unsafe_allow_html=True
)