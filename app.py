#!/usr/bin/env python
# coding: utf-8

import subprocess
import os
import io
import logging
import shutil
import uuid
import tomllib
import time
import streamlit as st
from pathlib import Path
from importlib import resources
from contextlib import redirect_stdout, redirect_stderr

from gbdraw.circular import circular_main
from gbdraw.linear import linear_main


# --- Basic Application Settings ---
st.set_page_config(layout="wide")

st.title("🧬 gbdraw Web App")
st.caption("A genome diagram generator for microbes and organelles")


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


st.markdown(
    """
    <p>
    <a href="https://anaconda.org/bioconda/gbdraw"><img src="https://anaconda.org/bioconda/gbdraw/badges/version.svg" alt="version"></a>
    <a href="https://anaconda.org/bioconda/gbdraw"><img src="https://anaconda.org/bioconda/gbdraw/badges/platforms.svg" alt="platforms"></a>
    <a href="http://bioconda.github.io/recipes/gbdraw/README.html"><img src="https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat" alt="install with bioconda"></a>
    <a href="https://anaconda.org/bioconda/gbdraw"><img src="https://anaconda.org/bioconda/gbdraw/badges/license.svg" alt="license"></a>
    <a href="https://deepwiki.com/satoshikawato/gbdraw"><img src="https://deepwiki.com/badge.svg" alt="Ask DeepWiki"></a>
    </p>
    """,
    unsafe_allow_html=True
)

# --- Temporary Directory and Session State Initialization ---
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

def create_manual_selectbox(label, options, key):
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
        args=(key,)
    )

PALETTES = get_palettes()

# --- Sidebar (File Management) ---
with st.sidebar:
    st.header("📂 File Management")
    uploaded_files_list = st.file_uploader(
        "Upload GenBank, Color, or BLAST files",
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

# --- Main Content (Mode Selection) ---
file_options = [""] + sorted(st.session_state.uploaded_files.keys())

# Replace st.tabs with st.radio to maintain state across reruns
selected_mode = st.radio(
    "Select Mode",
    ["🔵 Circular", "📏 Linear"],
    horizontal=True,
    label_visibility="collapsed"
)
st.markdown("---")


# --- CIRCULAR MODE ---
if selected_mode == "🔵 Circular":
    st.header("Circular Genome Map")
    st.subheader("Input Files")

    create_manual_selectbox("GenBank file:", file_options, "c_gb")
    create_manual_selectbox("Feature-specific color file (optional):", file_options, "c_t_color")
    st.markdown("---")

    st.subheader("🎨 Color Customization")
    
    # Callback for circular palette change
    def circular_palette_changed():
        new_palette = st.session_state.c_palette_selector
        st.session_state.custom_circular_colors = get_palette_colors(new_palette).copy()

    # Initialize state for circular colors if it doesn't exist
    if 'custom_circular_colors' not in st.session_state:
        default_palette = PALETTES[0] if PALETTES and PALETTES[0] != "" else "default"
        st.session_state.custom_circular_colors = get_palette_colors(default_palette).copy()

    st.selectbox(
        "Base color palette:", 
        PALETTES, 
        key="c_palette_selector",
        on_change=circular_palette_changed
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
                key=f"c_color_picker_{feature}"
            )
    st.markdown("---")

    with st.form("circular_form"):
        st.header("Drawing Options")
        col1, col2 = st.columns(2)
        with col1:
            st.subheader("Basic Settings")
            c_prefix = st.text_input("Output prefix (optional):", help="Default is input file name")
            c_species = st.text_input("Species name (optional):", help='e.g., "<i>Escherichia coli</i>"')
            c_strain = st.text_input("Strain/isolate name (optional):", help='e.g., "K-12"')
            c_fmt = st.selectbox("Output format:", ["svg", "png", "pdf", "eps", "ps"], index=0, key="c_fmt")
            c_track_type = st.selectbox("Track type:", ["tuckin", "middle", "spreadout"], index=0, key="c_track")
            c_legend = st.selectbox("Legend:", ["right", "left", "upper_left", "upper_right", "lower_left", "lower_right", "none"], index=0, key="c_legend")
        with col2:
            st.subheader("Display Options")
            c_show_labels = st.checkbox("Show labels", value=False, key="c_labels")
            c_separate_strands = st.checkbox("Separate strands", value=False, key="c_strands")
            c_allow_inner_labels = st.checkbox("Allow inner labels", value=False, key="c_inner_labels")
            if c_allow_inner_labels:
                st.info("💡 Inner labels are enabled. GC content and GC skew tracks will be automatically suppressed to avoid overlap.")
                c_suppress_gc = True
                c_suppress_skew = True
            else:
                c_suppress_gc = st.checkbox("Suppress GC content track", value=False, key="c_gc_suppress")
                c_suppress_skew = st.checkbox("Suppress GC skew track", value=False, key="c_skew_suppress")
            with st.expander("🔧 Advanced Options"):
                c_adv_feat = st.multiselect("Features (-k):", options=FEATURE_KEYS, default=["CDS", "tRNA", "rRNA", "repeat_region"], key="c_feat")
                c_adv_nt = st.text_input("Dinucleotide (--nt):", value="GC", key="c_nt")
                c_adv_win = st.number_input("Window size:", value=1000, key="c_win")
                c_adv_step = st.number_input("Step size:", value=100, key="c_step")
                c_adv_blk_color = st.color_picker("Block stroke color:", value="#808080", key="c_b_color")
                c_adv_blk_width = st.number_input("Block stroke width:", 0.0, key="c_b_width")
                c_adv_line_color = st.color_picker("Line stroke color:", value="#808080", key="c_l_color")
                c_adv_line_width = st.number_input("Line stroke width:", 1.0, key="c_l_width")
                c_adv_label_font_size = st.number_input("Label font size (default: 8 pt (>=50 kb) or 16 pt (<50 kb):", key="c_label_font_size")
        c_submitted = st.form_submit_button("🚀 Run gbdraw Circular", type="primary")

    if c_submitted:
        selected_gb_file = st.session_state.get("c_gb_manual", "")
        if not selected_gb_file:
            st.error("Please select a GenBank file.")
        else:
            gb_path = st.session_state.uploaded_files[selected_gb_file]
            sanitized_prefix = os.path.basename(c_prefix.strip())
            prefix = sanitized_prefix or Path(selected_gb_file).stem
            circular_args = ["-i", gb_path, "-o", prefix, "-f", c_fmt, "--track_type", c_track_type]
            if c_species:
                circular_args.extend(["--species", c_species])
            if c_strain:
                circular_args.extend(["--strain", c_strain])
            if c_show_labels: circular_args.append("--show_labels")
            if c_separate_strands: circular_args.append("--separate_strands")
            if c_allow_inner_labels:
                circular_args.extend(["--allow_inner_labels", "--suppress_gc", "--suppress_skew"])
            if c_adv_label_font_size:
                circular_args += ["--label_font_size", str(c_adv_label_font_size)]
            else:
                if c_suppress_gc: circular_args.append("--suppress_gc")
                if c_suppress_skew: circular_args.append("--suppress_skew")
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
            circular_args += ["-k", ",".join(c_adv_feat), "-n", c_adv_nt, "-w", str(c_adv_win), "-s", str(c_adv_step)]
            circular_args += ["--block_stroke_color", c_adv_blk_color, "--block_stroke_width", str(c_adv_blk_width)]
            circular_args += ["--line_stroke_color", c_adv_line_color, "--line_stroke_width", str(c_adv_line_width)]
            selected_t_color_file = st.session_state.get("c_t_color_manual", "")
            if selected_t_color_file: circular_args += ["-t", st.session_state.uploaded_files[selected_t_color_file]]
            
            logger = logging.getLogger() 
            log_capture = io.StringIO()  
            command_str = f"gbdraw circular {' '.join(circular_args)}"
            log_capture.write(f"--- Executed Command ---\n{command_str}\n------------------------\n\n")
            stream_handler = logging.StreamHandler(log_capture)
            stream_handler.setLevel(logging.INFO) 
            logger.addHandler(stream_handler)
            start_time = time.time()
            with st.spinner(f"Running: `gbdraw circular {' '.join(circular_args)}`"):
                try:
                    with redirect_stdout(log_capture), redirect_stderr(log_capture):
                        circular_main(circular_args)
                    st.success("✅ gbdraw finished successfully.")
                    end_time = time.time()
                    duration = end_time - start_time
                    log_capture.write(f"\n--- Execution Time ---\nTotal time: {duration:.2f} seconds\n----------------------")
                    st.session_state.circular_result = {"prefix": prefix, "fmt": c_fmt, "log": log_capture.getvalue()}

                except SystemExit as e:
                    if e.code != 0:
                        st.error(f"Error running gbdraw (exit code {e.code}):\n{log_capture.getvalue()}")
                        st.session_state.circular_result = None
                    else: 
                        st.success("✅ gbdraw finished successfully.")
                        end_time = time.time()
                        duration = end_time - start_time
                        log_capture.write(f"\n--- Execution Time ---\nTotal time: {duration:.2f} seconds\n----------------------")
                        st.session_state.circular_result = {"prefix": prefix, "fmt": c_fmt, "log": log_capture.getvalue()}

                except Exception as e:
                    st.error(f"An unexpected error occurred:\n{e}\n\nLog:\n{log_capture.getvalue()}")
                    st.session_state.circular_result = None
                finally:
                    logger.removeHandler(stream_handler)

    if st.session_state.circular_result:
        st.subheader("🌀 Circular Drawing Output")
        res = st.session_state.circular_result
        prefix = res["prefix"] 
        fmt = res["fmt"]
        output_files = sorted(list(Path(".").glob(f"{prefix}*.{fmt}")))
        if output_files:
            for out_path in output_files:
                file_extension = out_path.suffix.lower()
                if file_extension == ".svg":
                    st.image(out_path.read_text(), caption=str(out_path.name))
                elif file_extension == ".png":
                    st.image(str(out_path), caption=str(out_path.name))
                else:
                    st.info(f"📄 Preview is not available for {out_path.suffix.upper()} format. Please use the download button below.")
                with open(out_path, "rb") as f:
                    st.download_button(
                        f"⬇️ Download {out_path.name}",
                        data=f,
                        file_name=out_path.name,
                        key=f"download_{out_path.name}"
                    )
                st.markdown("---")
            with st.expander("Show Log"):
                st.text(res["log"])
        else:
            st.warning("Output file(s) seem to be missing. Please run again.")
            with st.expander("Show Log"):
                st.text(res["log"])

# --- LINEAR MODE ---
if selected_mode == "📏 Linear":
    st.header("Linear Genome Map")
    st.subheader("Input Files")
    input_container = st.container()
    with input_container:
        for i in range(st.session_state.linear_seq_count):
            cols = st.columns([3, 3])
            with cols[0]:
                create_manual_selectbox(f"Sequence File {i+1}", file_options, f"l_gb_{i}")
            if i < st.session_state.linear_seq_count - 1:
                with cols[1]:
                    create_manual_selectbox(f"Comparison File {i+1}", file_options, f"l_blast_{i}")

    b_col1, b_col2, _ = st.columns([1, 2, 5])
    if b_col1.button("➕ Add Pair"):
        st.session_state.linear_seq_count += 1
        st.rerun()
    if b_col2.button("➖ Remove Last Pair") and st.session_state.linear_seq_count > 1:
        st.session_state.linear_seq_count -= 1
        for key in [f"l_gb_{st.session_state.linear_seq_count}", f"l_blast_{st.session_state.linear_seq_count - 1}"]:
            manual_key = f"{key}_manual"
            if manual_key in st.session_state:
                del st.session_state[manual_key]
            if key in st.session_state:
                del st.session_state[key]
        st.rerun()

    st.subheader("Color Options")
    create_manual_selectbox("Feature-specific color file (optional):", file_options, "l_t_color")

    def linear_palette_changed():
        new_palette = st.session_state.l_palette_selector
        st.session_state.custom_linear_colors = get_palette_colors(new_palette).copy()

    if 'custom_linear_colors' not in st.session_state:
        default_palette = PALETTES[0] if PALETTES and PALETTES[0] != "" else "default"
        st.session_state.custom_linear_colors = get_palette_colors(default_palette).copy()

    st.selectbox(
        "Base color palette:",
        PALETTES,
        key="l_palette_selector",
        on_change=linear_palette_changed,
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
                key=f"l_color_picker_{feature}"
            )
    st.markdown("---")

    with st.form("linear_form"):
        st.header("Drawing Options")
        col1, col2 = st.columns(2)
        with col1:
            l_prefix = st.text_input("Output prefix:", value="linear", key="l_prefix")
            l_fmt = st.selectbox("Output format:", ["svg", "png", "pdf", "eps", "ps"], index=0, key="l_fmt")
            l_legend = st.selectbox("Legend:", ["right", "left", "none"], index=0, key="l_legend")
        with col2:
            l_show_labels = st.checkbox("Show labels", value=False, key="l_labels")
            l_separate_strands = st.checkbox("Separate strands", value=False, key="l_strands")
            l_align_center = st.checkbox("Align center", value=False, key="l_align")
            l_show_gc = st.checkbox("Show GC content", value=False, key="l_gc")
            l_resolve_overlaps = st.checkbox("Resolve overlaps (experimental)", value=False, key="l_overlaps")
            with st.expander("🔧 Advanced Options"):
                l_adv_feat = st.multiselect("Features (-k):", options=FEATURE_KEYS, default=["CDS", "tRNA", "rRNA", "repeat_region"], key="l_feat")
                l_adv_nt = st.text_input("nt (--nt):", value="GC", key="l_nt")
                l_adv_win = st.number_input("Window size:", value=1000, key="l_win")
                l_adv_step = st.number_input("Step size:", value=100, key="l_step")
                st.markdown("---")
                st.write("Comparison Filters:")
                l_adv_bitscore = st.number_input("Min bitscore:", value=50.0, key="l_bitscore")
                l_adv_evalue = st.text_input("Max E-value:", value="1e-2", key="l_evalue")
                l_adv_identity = st.number_input("Min identity (%):", value=0.0, key="l_identity")
                st.markdown("---")
                l_adv_blk_color = st.color_picker("Block stroke color:", value="#808080", key="l_b_color")
                l_adv_blk_width = st.number_input("Block stroke width:", 0.0, key="l_b_width")
                l_adv_line_color = st.color_picker("Line stroke color:", value="#808080", key="l_l_color")
                l_adv_line_width = st.number_input("Line stroke width:", 1.0, key="l_l_width")
                l_adv_label_font_size = st.number_input("Label font size (default: 5 pt (>=50 kb) or 16 pt (<50 kb):", key="l_label_font_size")
        l_submitted = st.form_submit_button("🚀 Run gbdraw Linear", type="primary")

    if l_submitted:
        selected_gb = [st.session_state.get(f"l_gb_{i}_manual", "") for i in range(st.session_state.linear_seq_count)]
        selected_gb = [f for f in selected_gb if f]
        selected_blast = [st.session_state.get(f"l_blast_{i}_manual", "") for i in range(st.session_state.linear_seq_count - 1)]
        selected_blast = [f for f in selected_blast if f]

        if not selected_gb:
            st.error("Please select at least one Sequence file.")
        elif selected_blast and len(selected_blast) != len(selected_gb) - 1:
            st.error(f"Please provide {len(selected_gb) - 1} comparison file(s) for {len(selected_gb)} sequence files.")
        else:
            gb_paths = [st.session_state.uploaded_files[f] for f in selected_gb]
            sanitized_prefix = os.path.basename(l_prefix.strip())
            prefix = sanitized_prefix or "linear"
            output_path = Path(f"{prefix}.{l_fmt}")
            linear_args = ["-i", *gb_paths, "-o", prefix, "-f", l_fmt]
            if selected_blast:
                blast_paths = [st.session_state.uploaded_files[f] for f in selected_blast]
                linear_args += ["-b", *blast_paths]
            if l_show_labels: linear_args.append("--show_labels")
            if l_separate_strands: linear_args.append("--separate_strands")
            if l_align_center: linear_args.append("--align_center")
            if l_show_gc: linear_args.append("--show_gc")
            if l_resolve_overlaps: linear_args.append("--resolve_overlaps")
            if l_legend != "right": linear_args += ["-l", l_legend]
            if l_adv_label_font_size:
                linear_args += ["--label_font_size", str(l_adv_label_font_size)]
            selected_palette = st.session_state.get("l_palette_selector")
            if selected_palette: linear_args += ["--palette", selected_palette]

            if 'custom_linear_colors' in st.session_state:
                custom_color_filename = f"custom_colors_l_{uuid.uuid4().hex[:8]}.tsv"
                save_path = UPLOAD_DIR / custom_color_filename
                with open(save_path, "w") as f:
                    for feature, color in st.session_state.custom_linear_colors.items():
                        f.write(f"{feature}\t{color}\n")
                linear_args += ["-d", str(save_path)]
            
            # FIX: Use linear-specific variables (l_adv_win, l_adv_step) instead of circular ones.
            linear_args += ["-k", ",".join(l_adv_feat), "-n", l_adv_nt, "-w", str(l_adv_win), "-s", str(l_adv_step)]
            
            linear_args += ["--bitscore", str(l_adv_bitscore), "--evalue", l_adv_evalue, "--identity", str(l_adv_identity)]
            linear_args += ["--block_stroke_color", l_adv_blk_color, "--block_stroke_width", str(l_adv_blk_width)]
            linear_args += ["--line_stroke_color", l_adv_line_color, "--line_stroke_width", str(l_adv_line_width)]
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
            with st.spinner(f"Running: `{command_str}`"):
                try:
                    with redirect_stdout(log_capture), redirect_stderr(log_capture):
                        linear_main(linear_args)
                    st.success("✅ gbdraw finished successfully.")
                    end_time = time.time()
                    duration = end_time - start_time
                    log_capture.write(f"\n--- Execution Time ---\nTotal time: {duration:.2f} seconds\n----------------------")
                    st.session_state.linear_result = {"path": output_path, "log": log_capture.getvalue()}

                except SystemExit as e:
                    if e.code != 0:
                        st.error(f"Error running gbdraw (exit code {e.code}):\n{log_capture.getvalue()}")
                        st.session_state.linear_result = None
                    else:
                        st.success("✅ gbdraw finished successfully.")
                        end_time = time.time()
                        duration = end_time - start_time
                        log_capture.write(f"\n--- Execution Time ---\nTotal time: {duration:.2f} seconds\n----------------------")
                        st.session_state.linear_result = {"path": output_path, "log": log_capture.getvalue()}

                except Exception as e:
                    st.error(f"An unexpected error occurred:\n{e}\n\nLog:\n{log_capture.getvalue()}")
                    st.session_state.linear_result = None
                finally:
                    logger.removeHandler(stream_handler)

    if st.session_state.linear_result:
        st.subheader("📏 Linear Drawing Output")
        res = st.session_state.linear_result
        out_path = res["path"]
        if out_path.exists():
            file_extension = out_path.suffix.lower()
            if file_extension == ".svg":
                st.image(out_path.read_text(), caption=str(out_path.name))
            elif file_extension == ".png":
                st.image(str(out_path), caption=str(out_path.name))
            else:
                st.info(f"📄 Preview is not available for {out_path.suffix.upper()} format. Please use the download button below.")
            with open(out_path, "rb") as f:
                st.download_button(
                    f"⬇️ Download {out_path.name}",
                    data=f,
                    file_name=out_path.name
                )
            with st.expander("Show Log"):
                st.text(res["log"])
        else:
            st.warning("Output file seems to be missing. Please run again.")
            with st.expander("Show Log"):
                st.text(res["log"])

# --- Footer ---
st.markdown("---")
st.markdown(
    "Author: [Satoshi Kawato](https://github.com/satoshikawato)  |  "
    "Source: [gbdraw](https://github.com/satoshikawato/gbdraw)",
    unsafe_allow_html=True
)
