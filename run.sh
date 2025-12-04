#!/bin/bash

# --server.enableCORS=false と --server.enableXsrfProtection=false を追加
streamlit run app.py \
    --server.port=8501 \
    --server.address=127.0.0.1 \
    --server.headless=true \
    --server.enableCORS=false \
    --server.enableXsrfProtection=false \
    --server.fileWatcherType=none \
    --browser.gatherUsageStats=false &

sleep 2

nginx -g 'daemon off;'

