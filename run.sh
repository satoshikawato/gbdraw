#!/bin/bash

# 1. Streamlit をバックグラウンドで起動 (ポート8501, ローカル接続のみ)
#    --server.headless=true をつけると非対話モードで安定します
streamlit run app.py \
    --server.port=8501 \
    --server.address=127.0.0.1 \
    --server.headless=true \
    --server.fileWatcherType=none \
    --browser.gatherUsageStats=false &

# 2. Streamlitが起きるまで少し待つ (簡易的なヘルスチェック待ち)
sleep 2

# 3. Nginx をフォアグラウンドで起動 (ポート8080)
#    Cloud Run はこのプロセスを見張り続けます
nginx -g 'daemon off;'

