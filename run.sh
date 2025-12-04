#!/bin/bash
set -e

# Start Nginx in the background
nginx

# Start Streamlit on the internal port Nginx proxies to
streamlit run app.py \
    --server.port=8501 \
    --server.address=0.0.0.0 \
    --server.enableWebsocketCompression=false

