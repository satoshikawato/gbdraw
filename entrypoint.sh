#!/bin/bash
# Streamlit をバックグラウンド(ポート8501)で起動
streamlit run app.py --server.port=8501 --server.address=0.0.0.0 --server.headless=true &

# Nginx をフォアグラウンド(ポート8080)で起動
nginx -g 'daemon off;'

