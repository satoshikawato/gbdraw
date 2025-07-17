import streamlit as st
import subprocess
import os
import shutil
import uuid
from pathlib import Path
import tomllib  
from importlib import resources

# --- アプリケーションの基本設定 ---
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

# --- 一時ディレクトリとセッション状態の初期化 ---
TEMP_DIR = Path("gbdraw_temp")
UPLOAD_DIR = TEMP_DIR / "uploads"
UPLOAD_DIR.mkdir(exist_ok=True, parents=True)

# セッション状態の初期化
if 'uploaded_files' not in st.session_state:
    st.session_state.uploaded_files = {}
if 'circular_result' not in st.session_state:
    st.session_state.circular_result = None
if 'linear_result' not in st.session_state:
    st.session_state.linear_result = None
if 'linear_seq_count' not in st.session_state:
    st.session_state.linear_seq_count = 1

# --- ヘルパー関数 ---
@st.cache_data
def get_palettes():
    """gbdrawの内部ファイルから動的にカラーパレットのリストを取得する"""
    try:
        # tomllib を使用
        with resources.files("gbdraw").joinpath("data").joinpath("color_palettes.toml").open("rb") as fh:
            doc = tomllib.load(fh)
        return [""] + sorted(k for k in doc if k != "title")
    except (FileNotFoundError, ModuleNotFoundError, AttributeError):
        st.warning("Could not dynamically load palettes from gbdraw. Using a default list.")
        return ["", "default", "paired", "pastel1", "pastel2", "set1", "set2", "set3", "tab10", "tab20", "tab20b", "tab20c"]

PALETTES = get_palettes()

# --- サイドバー (ファイル管理) ---
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
        if UPLOAD_DIR.exists():
            shutil.rmtree(UPLOAD_DIR)
        UPLOAD_DIR.mkdir(exist_ok=True, parents=True)
        st.session_state.uploaded_files = {}
        st.success("All uploaded files have been cleared.")
        st.rerun() # experimental_rerun を rerun に変更

# --- メインコンテンツ（タブ） ---
file_options = [""] + sorted(st.session_state.uploaded_files.keys())
tab_circular, tab_linear = st.tabs(["🔵 Circular", "📏 Linear"])

# --- CIRCULARタブ ---
with tab_circular:
    st.header("Circular Genome Map")
    with st.form("circular_form"):
        # (Circularモードのフォーム部分は変更なし)
        col1, col2 = st.columns(2)
        with col1:
            st.subheader("Basic Settings")
            c_gb_file = st.selectbox("GenBank file:", file_options, key="c_gb")
            c_prefix = st.text_input("Output prefix (optional):", help="Default is input file name")
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
                c_adv_blk_color = st.text_input("Block stroke color:", "black", key="c_b_color")
                c_adv_blk_width = st.number_input("Block stroke width:", 0.0, key="c_b_width")
                c_adv_line_color = st.text_input("Line stroke color:", "gray", key="c_l_color")
                c_adv_line_width = st.number_input("Line stroke width:", 1.0, key="c_l_width")
                c_mod_default_colors = st.selectbox("Custom default color file:", file_options, key="c_d_color")
                c_feature_specific_color_table = st.selectbox("Feature-specific color file:", file_options, key="c_t_color")
        c_submitted = st.form_submit_button("🚀 Run gbdraw Circular", type="primary")

    if c_submitted:
        # (Circularモードの実行ロジックは変更なし)
        if not c_gb_file:
            st.error("Please select a GenBank file.")
        else:
            gb_path = st.session_state.uploaded_files[c_gb_file]
            prefix = c_prefix.strip() or Path(c_gb_file).stem
            output_path = Path(f"{prefix}.{c_fmt}")
            cmd = ["gbdraw", "circular", "-i", gb_path, "-o", prefix, "-f", c_fmt, "--track_type", c_track_type]
            if c_show_labels: cmd.append("--show_labels")
            if c_separate_strands: cmd.append("--separate_strands")
            if c_legend != "right": cmd += ["-l", c_legend]
            if c_palette: cmd += ["--palette", c_palette]
            cmd += ["-k", c_adv_feat, "-n", c_adv_nt, "-w", str(c_adv_win), "-s", str(c_adv_step)]
            cmd += ["--block_stroke_color", c_adv_blk_color, "--block_stroke_width", str(c_adv_blk_width)]
            cmd += ["--line_stroke_color", c_adv_line_color, "--line_stroke_width", str(c_adv_line_width)]
            if c_mod_default_colors: cmd += ["-d", st.session_state.uploaded_files[c_mod_default_colors]]
            if c_feature_specific_color_table: cmd += ["-t", st.session_state.uploaded_files[c_feature_specific_color_table]]
            with st.spinner(f"Running: `{' '.join(cmd)}`"):
                result = subprocess.run(cmd, capture_output=True, text=True)
                if result.returncode != 0:
                    st.error(f"Error running gbdraw:\n{result.stderr}")
                    st.session_state.circular_result = None
                else:
                    st.success("✅ gbdraw finished successfully.")
                    st.session_state.circular_result = {"path": output_path, "log": result.stdout}

    if st.session_state.circular_result:
        # (Circularモードの結果表示ロジックは変更なし)
        st.subheader("🌀 Circular Drawing Output")
        res = st.session_state.circular_result
        out_path = res["path"]
        if out_path.exists():
            if out_path.suffix.lower() == ".svg": st.image(out_path.read_text(), caption=str(out_path.name))
            else: st.image(str(out_path), caption=str(out_path.name))
            with open(out_path, "rb") as f: st.download_button(f"⬇️ Download {out_path.name}", data=f, file_name=out_path.name)
            with st.expander("Show Log"): st.text(res["log"])
        else: st.warning("Output file seems to be missing. Please run again.")

# --- LINEAR TAB ---
with tab_linear:
    st.header("Linear Genome Map")
    st.subheader("Input Files")
    input_container = st.container()

    with input_container:
        for i in range(st.session_state.linear_seq_count):
            cols = st.columns([3, 3])
            with cols[0]:
                # --- START: 選択状態を確実に維持するための修正 ---
                gb_key = f"l_gb_{i}"
                # session_stateから現在の選択値を取得
                current_gb_selection = st.session_state.get(gb_key, "")
                try:
                    # 現在利用可能な選択肢リスト(file_options)から、現在の選択値が何番目にあるかを探す
                    gb_index = file_options.index(current_gb_selection)
                except ValueError:
                    # もし選択値がリストになければ（例：ファイルを全削除した場合）、先頭（空白）を選択する
                    gb_index = 0
                
                st.selectbox(
                    f"Sequence File {i+1}",
                    file_options,
                    index=gb_index,  # インデックスを明示的に指定
                    key=gb_key,
                )
                # --- END ---

            if i < st.session_state.linear_seq_count - 1:
                with cols[1]:
                    # --- START: comparisonファイルも同様に修正 ---
                    blast_key = f"l_blast_{i}"
                    current_blast_selection = st.session_state.get(blast_key, "")
                    try:
                        blast_index = file_options.index(current_blast_selection)
                    except ValueError:
                        blast_index = 0

                    st.selectbox(
                        f"Comparison File {i+1}",
                        file_options,
                        index=blast_index, # インデックスを明示的に指定
                        key=blast_key,
                    )
                    # --- END ---

    b_col1, b_col2, _ = st.columns([1, 2, 5])
    if b_col1.button("➕ Add Pair"):
        st.session_state.linear_seq_count += 1
        st.rerun()
    if b_col2.button("➖ Remove Last Pair") and st.session_state.linear_seq_count > 1:
        last_seq_key = f"l_gb_{st.session_state.linear_seq_count - 1}"
        last_blast_key = f"l_blast_{st.session_state.linear_seq_count - 2}"
        if last_seq_key in st.session_state:
            del st.session_state[last_seq_key]
        if last_blast_key in st.session_state:
            del st.session_state[last_blast_key]
        
        st.session_state.linear_seq_count -= 1
        st.rerun()

    with st.form("linear_form"):
        col1, col2 = st.columns(2)
        with col1:
            st.subheader("Basic Settings")
            l_prefix = st.text_input("Output prefix:", value="linear", key="l_prefix")
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
                l_mod_default_colors = st.selectbox("Custom default color file:", file_options, key="l_d_color")
                l_feature_specific_color_table = st.selectbox("Feature-specific color file:", file_options, key="l_t_color")
        l_submitted = st.form_submit_button("🚀 Run gbdraw Linear", type="primary")

    if l_submitted:
        selected_gb = [
            st.session_state[f"l_gb_{i}"]
            for i in range(st.session_state.linear_seq_count)
            if st.session_state.get(f"l_gb_{i}")
        ]
        selected_blast = [
            st.session_state[f"l_blast_{i}"]
            for i in range(st.session_state.linear_seq_count - 1)
            if st.session_state.get(f"l_blast_{i}")
        ]

        if not selected_gb:
            st.error("Please select at least one Sequence file.")
        elif selected_blast and len(selected_blast) != len(selected_gb) - 1:
            st.error(f"Please provide {len(selected_gb) - 1} comparison file(s) for {len(selected_gb)} sequence files.")
        else:
            gb_paths = [st.session_state.uploaded_files[f] for f in selected_gb]
            prefix = l_prefix.strip() or "linear"
            output_path = Path(f"{prefix}.{l_fmt}")
            cmd = ["gbdraw", "linear", "-i", *gb_paths, "-o", prefix, "-f", l_fmt]
            if selected_blast:
                blast_paths = [st.session_state.uploaded_files[f] for f in selected_blast]
                cmd += ["-b", *blast_paths]
            if l_show_labels: cmd.append("--show_labels")
            if l_separate_strands: cmd.append("--separate_strands")
            if l_align_center: cmd.append("--align_center")
            if l_show_gc: cmd.append("--show_gc")
            if l_resolve_overlaps: cmd.append("--resolve_overlaps")
            if l_legend != "right": cmd += ["-l", l_legend]
            if l_palette: cmd += ["--palette", l_palette]
            cmd += ["-k", l_adv_feat, "-n", l_adv_nt, "-w", str(l_adv_win), "-s", str(l_adv_step)]
            cmd += ["--bitscore", str(l_adv_bitscore), "--evalue", l_adv_evalue, "--identity", str(l_adv_identity)]
            cmd += ["--block_stroke_color", l_adv_blk_color, "--block_stroke_width", str(l_adv_blk_width)]
            cmd += ["--line_stroke_color", l_adv_line_color, "--line_stroke_width", str(l_adv_line_width)]
            if l_mod_default_colors: cmd += ["-d", st.session_state.uploaded_files[l_mod_default_colors]]
            if l_feature_specific_color_table: cmd += ["-t", st.session_state.uploaded_files[l_feature_specific_color_table]]
            with st.spinner(f"Running: `{' '.join(cmd)}`"):
                result = subprocess.run(cmd, capture_output=True, text=True)
                if result.returncode != 0:
                    st.error(f"Error running gbdraw:\n{result.stderr}")
                    st.session_state.linear_result = None
                else:
                    st.success("✅ gbdraw finished successfully.")
                    st.session_state.linear_result = {"path": output_path, "log": result.stdout}

    if st.session_state.linear_result:
        st.subheader("📏 Linear Drawing Output")
        res = st.session_state.linear_result
        out_path = res["path"]
        if out_path.exists():
            if out_path.suffix.lower() == ".svg": st.image(out_path.read_text(), caption=str(out_path.name))
            else: st.image(str(out_path), caption=str(out_path.name))
            with open(out_path, "rb") as f: st.download_button(f"⬇️ Download {out_path.name}", data=f, file_name=out_path.name)
            with st.expander("Show Log"): st.text(res["log"])
        else: st.warning("Output file seems to be missing. Please run again.")

# --- フッター ---
st.markdown("---")
st.markdown(
    "Author: [Satoshi Kawato](https://github.com/satoshikawato)  |  "
    "Source: [gbdraw](https://github.com/satoshikawato/gbdraw)",
    unsafe_allow_html=True
)