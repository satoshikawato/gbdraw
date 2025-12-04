FROM python:3.13-slim

ENV PYTHONUNBUFFERED=1
WORKDIR /app

# 1. Nginx を追加インストール
RUN apt-get update && apt-get install -y \
    git \
    libcairo2 \
    libpango-1.0-0 \
    libpangocairo-1.0-0 \
    gdk-pixbuf2.0-0 \
    libffi-dev \
    shared-mime-info \
    curl \
    fonts-liberation \
    fontconfig \
    nginx \
    && rm -rf /var/lib/apt/lists/*

COPY . .

# 2. フォント登録 (元のまま)
RUN mkdir -p /usr/share/fonts/truetype/gbdraw \
    && cp gbdraw/data/*.ttf /usr/share/fonts/truetype/gbdraw/ 2>/dev/null || true \
    && fc-cache -f -v

RUN pip3 install --no-cache-dir -r requirements.txt

# Google Analytics (元のまま)
RUN sed -i 's~<head>~<head><script async src="https://www.googletagmanager.com/gtag/js?id=G-GG6JMKM02Y"></script><script>window.dataLayer = window.dataLayer || [];function gtag(){dataLayer.push(arguments);}gtag("js", new Date());gtag("config", "G-GG6JMKM02Y");</script>~' /usr/local/lib/python3.13/site-packages/streamlit/static/index.html

# 3. Nginxの設定ファイルを配置
COPY nginx.conf /etc/nginx/nginx.conf


RUN mkdir -p /root/.streamlit
COPY .streamlit/config.toml /root/.streamlit/config.toml

# 4. 起動スクリプトに実行権限を付与
COPY run.sh /app/run.sh
RUN chmod +x /app/run.sh

# Cloud Run は 8080 を見に来る
EXPOSE 8080

# ヘルスチェック: Nginx経由で確認するように変更
HEALTHCHECK CMD curl --fail http://localhost:8080/_stcore/health || exit 1

# 5. エントリーポイントをシェルスクリプトに変更
ENTRYPOINT ["/app/run.sh"]


