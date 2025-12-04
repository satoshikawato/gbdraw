#!/bin/bash

# Nginxをバックグラウンドで起動
nginx -c /app/nginx.conf &

# Streamlitを起動（ポートは8501に固定）
# maxUploadSizeはここで指定
streamlit run app.py \
    --server.port=8501 \
    --server.address=127.0.0.1 \
    --server.maxUploadSize=500 \
    --server.enableCORS=false \
    --server.enableXsrfProtection=false
