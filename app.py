import streamlit as st
import subprocess
import os
import tempfile
import shutil
import uuid
from pathlib import Path

# Title
st.title("🧬 gbdraw: GenBank Genome Visualizer")
st.markdown(
    """
    Upload a GenBank file to visualize your genome in **linear** or **circular** mode.
    Customize colors, tracks, and advanced options.
    """
)

# Temporary workspace
WORKDIR = tempfile.mkdtemp(prefix="gbdraw_")

# Upload GenBank
uploaded_file = st.file_uploader("📄 Upload GenBank file", type=["gb", "gbk", "genbank"])
if uploaded_file:
    gbk_filename = os.path.join(WORKDIR, f"{uuid.uuid4().hex}.gb")
    with open(gbk_filename, "wb") as f:
        f.write(uploaded_file.getbuffer())
    st.success(f"✅ Uploaded {uploaded_file.name}")

# Visualization options
mode = st.radio("Genome Visualization Mode", ["circular", "linear"])
track_type = st.selectbox("Track Type", ["tuckin", "middle", "spreadout"], index=0)
output_format = st.selectbox("Output Format", ["svg", "png", "pdf", "eps"], index=1)
show_labels = st.checkbox("Show Feature Labels", value=True)
separate_strands = st.checkbox("Separate Forward/Reverse Strands", value=False)
legend_position = st.selectbox(
    "Legend Position",
    ["right", "left", "upper_left", "upper_right", "lower_left", "lower_right", "none"],
    index=0
)

# Advanced options
with st.expander("🔧 Advanced Options"):
    palette = st.text_input("Color Palette (leave blank for default):", "")
    features = st.text_input("Features to include (comma-separated)", "CDS,tRNA,rRNA")
    nt = st.text_input("GC Skew dinucleotide (--nt)", "GC")
    window = st.number_input("GC Skew window size", value=1000, min_value=100, step=100)
    step = st.number_input("GC Skew step size", value=100, min_value=10, step=10)
    block_stroke_color = st.text_input("Block Stroke Color", "black")
    block_stroke_width = st.number_input("Block Stroke Width", value=0.0, step=0.1)
    line_stroke_color = st.text_input("Line Stroke Color", "gray")
    line_stroke_width = st.number_input("Line Stroke Width", value=1.0, step=0.1)

# Run gbdraw
if st.button("🚀 Run gbdraw") and uploaded_file:
    # Build command
    output_prefix = os.path.join(WORKDIR, Path(uploaded_file.name).stem)
    cmd = [
        "gbdraw", mode,
        "-i", gbk_filename,
        "-f", output_format,
        "--track_type", track_type,
        "-o", output_prefix
    ]

    if show_labels:
        cmd.append("--show_labels")
    if separate_strands:
        cmd.append("--separate_strands")
    if legend_position != "right":
        cmd += ["-l", legend_position]
    if palette.strip():
        cmd += ["--palette", palette]
    cmd += [
        "-n", nt,
        "-w", str(window),
        "-s", str(step),
        "--block_stroke_color", block_stroke_color,
        "--block_stroke_width", str(block_stroke_width),
        "--line_stroke_color", line_stroke_color,
        "--line_stroke_width", str(line_stroke_width),
        "-k", features
    ]

    st.code(" ".join(cmd), language="bash")

    # Execute
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            st.error(f"❌ gbdraw failed:\n{result.stderr}")
        else:
            st.success("✅ gbdraw finished successfully")
            output_file = f"{output_prefix}.{output_format}"
            if output_format in ("png", "svg"):
                st.image(output_file, caption="Genome Visualization")
            st.download_button(
                "⬇ Download Output File",
                open(output_file, "rb").read(),
                file_name=f"{Path(uploaded_file.name).stem}.{output_format}"
            )
    except FileNotFoundError:
        st.error("❌ gbdraw not found. Make sure it is installed in this environment.")

