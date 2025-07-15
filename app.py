import streamlit as st
import subprocess
import os
import shutil
import uuid
from pathlib import Path

# 一時ディレクトリの作成
TEMP_DIR = Path("gbdraw_temp")
TEMP_DIR.mkdir(exist_ok=True)

# ユーザーファイル保存用ディレクトリ
UPLOAD_DIR = TEMP_DIR / "uploads"
UPLOAD_DIR.mkdir(exist_ok=True)

# 一時ファイル掃除
def cleanup():
    shutil.rmtree(TEMP_DIR, ignore_errors=True)

# タイトル
st.title("🧬 gbdraw Web App")
st.caption("Draw genome maps in circular or linear mode")

# --- Tabs for circular / linear ---
tab_circular, tab_linear = st.tabs(["🔵 Circular", "📏 Linear"])

# --- CIRCULAR TAB ---
with tab_circular:
    st.header("Circular Genome Map")
    gb_file_circular = st.file_uploader("Upload GenBank file", type=["gb", "gbk", "genbank"], key="circular_upload")
    prefix_circular = st.text_input("Output prefix (optional)", key="circular_prefix")
    fmt_circular = st.selectbox("Output format", ["svg", "png", "pdf", "eps", "ps"], index=0, key="circular_fmt")

    track_type = st.selectbox("Track type", ["tuckin", "middle", "spreadout"], index=0)
    show_labels = st.checkbox("Show labels", value=True, key="circular_labels")
    separate_strands = st.checkbox("Separate strands", value=False, key="circular_strands")
    palette = st.text_input("Palette name (optional)", value="default", key="circular_palette")

    if st.button("🚀 Run gbdraw Circular"):
        if gb_file_circular is None:
            st.error("Please upload a GenBank file first.")
        else:
            # Save uploaded file
            gb_path = UPLOAD_DIR / f"{uuid.uuid4().hex}_{gb_file_circular.name}"
            with open(gb_path, "wb") as f:
                f.write(gb_file_circular.read())

            # Build gbdraw command
            cmd = ["gbdraw", "circular", "-i", str(gb_path), "-f", fmt_circular]
            if prefix_circular.strip():
                cmd += ["-o", prefix_circular.strip()]
            cmd += ["--track_type", track_type]
            if show_labels:
                cmd.append("--show_labels")
            if separate_strands:
                cmd.append("--separate_strands")
            if palette and palette != "default":
                cmd += ["--palette", palette]

            st.info(f"🛠️ Running: {' '.join(cmd)}")

            # Run gbdraw
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                st.error(f"Error running gbdraw:\n{result.stderr}")
            else:
                st.success("✅ gbdraw finished successfully.")
                # Circular drawing output
                st.subheader("🌀 Circular Drawing Output")
                output_files = [p for p in Path(".").iterdir() if p.is_file()]
                if not output_files:
                    st.warning("No output files found.")
                else:
                    for out in output_files:
                        if out.suffix.lower() in [".png", ".svg", ".pdf", ".eps", ".ps"]:
                            st.image(str(out), caption=str(out.name))
                            st.download_button(
                                "⬇️ Download " + out.name,
                                data=open(out, "rb"),
                                file_name=out.name
                            )

    # --- LINEAR TAB ---
with tab_linear:
    st.header("Linear Genome Map")
    gb_file_linear = st.file_uploader("Upload GenBank file", type=["gb", "gbk", "genbank"], key="linear_upload")
    prefix_linear = st.text_input("Output prefix (optional)", key="linear_prefix")
    fmt_linear = st.selectbox("Output format", ["svg", "png", "pdf", "eps", "ps"], index=0, key="linear_fmt")

    show_labels_linear = st.checkbox("Show labels", value=True, key="linear_labels")
    palette_linear = st.text_input("Palette name (optional)", value="default", key="linear_palette")

    if st.button("🚀 Run gbdraw Linear"):
        if gb_file_linear is None:
            st.error("Please upload a GenBank file first.")
        else:
            # Save uploaded file
            gb_path = UPLOAD_DIR / f"{uuid.uuid4().hex}_{gb_file_linear.name}"
            with open(gb_path, "wb") as f:
                f.write(gb_file_linear.read())

            # Build gbdraw command
            cmd = ["gbdraw", "linear", "-i", str(gb_path), "-f", fmt_linear]
            if prefix_linear.strip():
                cmd += ["-o", prefix_linear.strip()]
            if show_labels_linear:
                cmd.append("--show_labels")
            if palette_linear and palette_linear != "default":
                cmd += ["--palette", palette_linear]

            st.info(f"🛠️ Running: {' '.join(cmd)}")

            # Run gbdraw
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                st.error(f"Error running gbdraw:\n{result.stderr}")
            else:
                st.success("✅ gbdraw finished successfully.")
                # Linear drawing output

                st.subheader("Linear Drawing Output")
                output_files = [p for p in Path(".").iterdir() if p.is_file()]
                if not output_files:
                    st.warning("No output files found.")
                else:
                    for out in output_files:
                        if out.suffix.lower() in [".png", ".svg", ".pdf", ".eps", ".ps"]:
                            st.image(str(out), caption=str(out.name))
                            st.download_button(
                                "⬇️ Download " + out.name,
                                data=open(out, "rb"),
                                file_name=out.name
                            )