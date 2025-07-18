import io
import logging
import shutil
import uuid
import tomllib
import streamlit as st
from pathlib import Path
from importlib import resources
from contextlib import redirect_stderr

from gbdraw.circular import circular_main
from gbdraw.linear import linear_main

# --- Basic Application Settings ---
st.set_page_config(layout="wide")

st.title("🧬 gbdraw Web App")
st.caption("A genome diagram generator for microbes and organelles")

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
        return ["", "default", "paired", "pastel1", "pastel2", "set1", "set2", "set3", "tab10", "tab20", "tab20b", "tab20c"]

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
        for fname in list(st.session_state.uploaded_files.keys()):
            st.write(f"- `{fname}`")

    if st.button("Clear All Uploaded Files"):
        if UPLOAD_DIR.exists():
            shutil.rmtree(UPLOAD_DIR)
        UPLOAD_DIR.mkdir(exist_ok=True, parents=True)
        st.session_state.uploaded_files = {}
        st.success("All uploaded files have been cleared.")
        st.rerun()

# --- Main Content (Tabs) ---
file_options = [""] + sorted(st.session_state.uploaded_files.keys())
tab_circular, tab_linear = st.tabs(["🔵 Circular", "📏 Linear"])

# --- CIRCULAR TAB ---
with tab_circular:
    st.header("Circular Genome Map")
    st.subheader("Input Files")
    c_gb_file = st.selectbox("GenBank file:", file_options, key="c_gb_selector")
    c_mod_default_colors = st.selectbox("Custom default color file (optional):", file_options, key="c_d_color_selector")
    c_feature_specific_color_table = st.selectbox("Feature-specific color file (optional):", file_options, key="c_t_color_selector")
    st.markdown("---")

    with st.form("circular_form"):
        st.header("Drawing Options")
        col1, col2 = st.columns(2)
        with col1:
            st.subheader("Basic Settings")
            c_prefix = st.text_input("Output prefix (optional):", help="Default is input file name. Directory traversal is disabled.")
            c_fmt = st.selectbox("Output format:", ["svg", "png", "pdf", "eps", "ps"], index=0, key="c_fmt")
            c_track_type = st.selectbox("Track type:", ["tuckin", "middle", "spreadout"], index=0, key="c_track")
            c_legend = st.selectbox("Legend:", ["right", "left", "upper_left", "upper_right", "lower_left", "lower_right", "none"], index=0, key="c_legend")
            c_palette = st.selectbox("Color palette:", PALETTES, key="c_palette")

        with col2:
            st.subheader("Display Options")
            c_show_labels = st.checkbox("Show labels", value=False, key="c_labels")
            c_separate_strands = st.checkbox("Separate strands", value=False, key="c_strands")
            with st.expander("🔧 Advanced Options"):
                c_adv_feat = st.text_input("Features (-k):", value="CDS,tRNA,rRNA,repeat_region", key="c_feat")
                c_adv_nt = st.text_input("Dinucleotide (--nt):", value="GC", key="c_nt")
                c_adv_win = st.number_input("Window size:", value=1000, key="c_win")
                c_adv_step = st.number_input("Step size:", value=100, key="c_step")
                c_adv_blk_color = st.text_input("Block stroke color:", "gray", key="c_b_color")
                c_adv_blk_width = st.number_input("Block stroke width:", 0.0, key="c_b_width")
                c_adv_line_color = st.text_input("Line stroke color:", "gray", key="c_l_color")
                c_adv_line_width = st.number_input("Line stroke width:", 1.0, key="c_l_width")
        
        c_submitted = st.form_submit_button("🚀 Run gbdraw Circular", type="primary")

    if c_submitted:
        if not c_gb_file:
            st.error("Please select a GenBank file.")
        else:
            safe_prefix = Path(c_prefix.strip() or Path(c_gb_file).stem).name
            output_path = TEMP_DIR / f"{safe_prefix}.{c_fmt}"
            
            # Explicit security check to satisfy static analysis tools
            try:
                trusted_dir = TEMP_DIR.resolve()
                absolute_output_path = output_path.resolve()
                if not str(absolute_output_path).startswith(str(trusted_dir)):
                    raise Exception("Attempted to write to a disallowed path.")
            except Exception as e:
                st.error(f"Security check failed: {e}")
                st.stop()
            
            gb_path = st.session_state.uploaded_files[c_gb_file]
            circular_args = ["-i", gb_path, "-o", str(output_path.with_suffix('')), "-f", c_fmt, "--track_type", c_track_type]
            if c_show_labels: circular_args.append("--show_labels")
            if c_separate_strands: circular_args.append("--separate_strands")
            if c_legend != "right": circular_args += ["-l", c_legend]
            if c_palette: circular_args += ["--palette", c_palette]
            circular_args += ["-k", c_adv_feat, "-n", c_adv_nt, "-w", str(c_adv_win), "-s", str(c_adv_step)]
            circular_args += ["--block_stroke_color", c_adv_blk_color, "--block_stroke_width", str(c_adv_blk_width)]
            circular_args += ["--line_stroke_color", c_adv_line_color, "--line_stroke_width", str(c_adv_line_width)]
            if c_mod_default_colors: circular_args += ["-d", st.session_state.uploaded_files[c_mod_default_colors]]
            if c_feature_specific_color_table: circular_args += ["-t", st.session_state.uploaded_files[c_feature_specific_color_table]]

            logger = logging.getLogger()
            log_capture = io.StringIO()
            command_str = f"gbdraw circular {' '.join(circular_args)}"
            log_capture.write(f"--- Executed Command ---\n{command_str}\n------------------------\n\n")
            stream_handler = logging.StreamHandler(log_capture)
            stream_handler.setLevel(logging.INFO)
            logger.addHandler(stream_handler)

            with st.spinner(f"Running: {command_str}"):
                try:
                    with redirect_stderr(log_capture):
                        circular_main(circular_args)
                    st.success("✅ gbdraw finished successfully.")
                    st.session_state.circular_result = {"path": output_path, "log": log_capture.getvalue()}
                except SystemExit as e:
                    if e.code != 0:
                        st.error(f"Error running gbdraw (exit code {e.code}):\n{log_capture.getvalue()}")
                        st.session_state.circular_result = None
                    else:
                        st.success("✅ gbdraw finished successfully.")
                        st.session_state.circular_result = {"path": output_path, "log": log_capture.getvalue()}
                except Exception as e:
                    st.error(f"An unexpected error occurred:\n{e}\n\nLog:\n{log_capture.getvalue()}")
                    st.session_state.circular_result = None
                finally:
                    logger.removeHandler(stream_handler)

    if st.session_state.circular_result:
        st.subheader("🌀 Circular Drawing Output")
        res = st.session_state.circular_result
        out_path = res["path"]
        if out_path.exists():
            if out_path.suffix.lower() == ".svg":
                st.image(out_path.read_text(), caption=str(out_path.name))
            elif out_path.suffix.lower() == ".png":
                st.image(str(out_path), caption=str(out_path.name))
            else:
                st.info(f"📄 Preview not available for {out_path.suffix.upper()}. Please download.")
            with open(out_path, "rb") as f:
                st.download_button(f"⬇️ Download {out_path.name}", data=f, file_name=out_path.name)
            with st.expander("Show Log"):
                st.text(res["log"])
        else:
            st.warning("Output file seems to be missing. Please run again.")

# --- LINEAR TAB ---
with tab_linear:
    st.header("Linear Genome Map")
    st.subheader("Input Files")
    for i in range(st.session_state.linear_seq_count):
        cols = st.columns([3, 3])
        with cols[0]:
            st.selectbox(f"Sequence File {i+1}", file_options, key=f"l_gb_{i}")
        if i < st.session_state.linear_seq_count - 1:
            with cols[1]:
                st.selectbox(f"Comparison File {i+1}", file_options, key=f"l_blast_{i}")
    b_col1, b_col2, _ = st.columns([1, 2, 5])
    if b_col1.button("➕ Add Pair"):
        st.session_state.linear_seq_count += 1
        st.rerun()
    if b_col2.button("➖ Remove Last Pair") and st.session_state.linear_seq_count > 1:
        st.session_state.linear_seq_count -= 1
        st.rerun()

    st.subheader("Custom Color Files (Optional)")
    l_mod_default_colors = st.selectbox("Custom default color file:", file_options, key="l_d_color_selector")
    l_feature_specific_color_table = st.selectbox("Feature-specific color file:", file_options, key="l_t_color_selector")
    st.markdown("---")

    with st.form("linear_form"):
        st.header("Drawing Options")
        col1, col2 = st.columns(2)
        with col1:
            st.subheader("Basic Settings")
            l_prefix = st.text_input("Output prefix:", value="linear", key="l_prefix", help="Directory traversal is disabled.")
            l_fmt = st.selectbox("Output format:", ["svg", "png", "pdf", "eps", "ps"], index=0, key="l_fmt")
            l_legend = st.selectbox("Legend:", ["right", "left", "none"], index=0, key="l_legend")
            l_palette = st.selectbox("Color palette:", PALETTES, key="l_palette")
        with col2:
            st.subheader("Display Options")
            l_show_labels = st.checkbox("Show labels", value=False, key="l_labels")
            l_separate_strands = st.checkbox("Separate strands", value=False, key="l_strands")
            l_align_center = st.checkbox("Align center", value=False, key="l_align")
            l_show_gc = st.checkbox("Show GC content", value=False, key="l_gc")
            l_resolve_overlaps = st.checkbox("Resolve overlaps (experimental)", value=False, key="l_overlaps")
            with st.expander("🔧 Advanced Options"):
                l_adv_feat = st.text_input("Features (-k):", value="CDS,tRNA,rRNA,repeat_region", key="l_feat")
                l_adv_nt = st.text_input("nt (--nt):", value="GC", key="l_nt")
                l_adv_win = st.number_input("Window size:", value=1000, key="l_win")
                l_adv_step = st.number_input("Step size:", value=100, key="l_step")
                st.markdown("---")
                st.write("Comparison Filters:")
                l_adv_bitscore = st.number_input("Min bitscore:", value=50.0, key="l_bitscore")
                l_adv_evalue = st.text_input("Max E-value:", value="1e-2", key="l_evalue")
                l_adv_identity = st.number_input("Min identity (%):", value=0.0, key="l_identity")
                st.markdown("---")
                l_adv_blk_color = st.text_input("Block stroke color:", "gray", key="l_b_color")
                l_adv_blk_width = st.number_input("Block stroke width:", 0.0, key="l_b_width")
                l_adv_line_color = st.text_input("Line stroke color:", "gray", key="l_l_color")
                l_adv_line_width = st.number_input("Line stroke width:", 1.0, key="l_l_width")
        
        l_submitted = st.form_submit_button("🚀 Run gbdraw Linear", type="primary")

    if l_submitted:
        selected_gb = [st.session_state[f"l_gb_{i}"] for i in range(st.session_state.linear_seq_count) if st.session_state.get(f"l_gb_{i}")]
        selected_blast = [st.session_state[f"l_blast_{i}"] for i in range(st.session_state.linear_seq_count - 1) if st.session_state.get(f"l_blast_{i}")]

        if not selected_gb:
            st.error("Please select at least one Sequence file.")
        elif selected_blast and len(selected_blast) != len(selected_gb) - 1:
            st.error(f"Please provide {len(selected_gb) - 1} comparison file(s) for {len(selected_gb)} sequence files.")
        else:
            safe_prefix = Path(l_prefix.strip() or "linear").name
            output_path = TEMP_DIR / f"{safe_prefix}.{l_fmt}"

            # Explicit security check to satisfy static analysis tools
            try:
                trusted_dir = TEMP_DIR.resolve()
                absolute_output_path = output_path.resolve()
                if not str(absolute_output_path).startswith(str(trusted_dir)):
                    raise Exception("Attempted to write to a disallowed path.")
            except Exception as e:
                st.error(f"Security check failed: {e}")
                st.stop()

            gb_paths = [st.session_state.uploaded_files[f] for f in selected_gb]
            linear_args = ["-i", *gb_paths, "-o", str(output_path.with_suffix('')), "-f", l_fmt]
            if selected_blast:
                blast_paths = [st.session_state.uploaded_files[f] for f in selected_blast]
                linear_args += ["-b", *blast_paths]
            if l_show_labels: linear_args.append("--show_labels")
            if l_separate_strands: linear_args.append("--separate_strands")
            if l_align_center: linear_args.append("--align_center")
            if l_show_gc: linear_args.append("--show_gc")
            if l_resolve_overlaps: linear_args.append("--resolve_overlaps")
            if l_legend != "right": linear_args += ["-l", l_legend]
            if l_palette: linear_args += ["--palette", l_palette]
            linear_args += ["-k", l_adv_feat, "-n", l_adv_nt, "-w", str(l_adv_win), "-s", str(l_adv_step)]
            linear_args += ["--bitscore", str(l_adv_bitscore), "--evalue", l_adv_evalue, "--identity", str(l_adv_identity)]
            linear_args += ["--block_stroke_color", l_adv_blk_color, "--block_stroke_width", str(l_adv_blk_width)]
            linear_args += ["--line_stroke_color", l_adv_line_color, "--line_stroke_width", str(l_adv_line_width)]
            if l_mod_default_colors: linear_args += ["-d", st.session_state.uploaded_files[l_mod_default_colors]]
            if l_feature_specific_color_table: linear_args += ["-t", st.session_state.uploaded_files[l_feature_specific_color_table]]
            
            logger = logging.getLogger()
            log_capture = io.StringIO()
            command_str = f"gbdraw linear {' '.join(linear_args)}"
            log_capture.write(f"--- Executed Command ---\n{command_str}\n------------------------\n\n")
            stream_handler = logging.StreamHandler(log_capture)
            stream_handler.setLevel(logging.INFO)
            logger.addHandler(stream_handler)
            
            with st.spinner(f"Running: {command_str}"):
                try:
                    with redirect_stderr(log_capture):
                        linear_main(linear_args)
                    st.success("✅ gbdraw finished successfully.")
                    st.session_state.linear_result = {"path": output_path, "log": log_capture.getvalue()}
                except SystemExit as e:
                    if e.code != 0:
                        st.error(f"Error running gbdraw (exit code {e.code}):\n{log_capture.getvalue()}")
                        st.session_state.linear_result = None
                    else:
                        st.success("✅ gbdraw finished successfully.")
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
            if out_path.suffix.lower() == ".svg":
                st.image(out_path.read_text(), caption=str(out_path.name))
            elif out_path.suffix.lower() == ".png":
                st.image(str(out_path), caption=str(out_path.name))
            else:
                st.info(f"📄 Preview not available for {out_path.suffix.upper()}. Please download.")
            with open(out_path, "rb") as f:
                st.download_button(f"⬇️ Download {out_path.name}", data=f, file_name=out_path.name)
            with st.expander("Show Log"):
                st.text(res["log"])
        else:
            st.warning("Output file seems to be missing. Please run again.")

# --- Footer ---
st.markdown("---")
st.markdown(
    "Author: [Satoshi Kawato](https://github.com/satoshikawato)  |  "
    "Source: [gbdraw](https://github.com/satoshikawato/gbdraw)",
    unsafe_allow_html=True
)