# Use Python 3.13 as the base image
FROM python:3.13-slim

ENV PYTHONUNBUFFERED=1
WORKDIR /app

# -----------------------------------------------------------
# 1. Install latest Nginx (Official Repo) & dependencies
# デフォルトの古いNginxではなく、公式の最新版を入れるための手順
# -----------------------------------------------------------
RUN apt-get update && apt-get install -y \
    curl gnupg2 ca-certificates lsb-release debian-archive-keyring \
    && curl https://nginx.org/keys/nginx_signing.key | gpg --dearmor \
    | tee /usr/share/keyrings/nginx-archive-keyring.gpg >/dev/null \
    && echo "deb [signed-by=/usr/share/keyrings/nginx-archive-keyring.gpg] \
    http://nginx.org/packages/debian bookworm nginx" \
    | tee /etc/apt/sources.list.d/nginx.list \
    && apt-get update && apt-get install -y nginx \
    && rm -rf /var/lib/apt/lists/*

# 2. Install library dependencies (gbdraw用)
RUN apt-get update && apt-get install -y \
    git \
    libcairo2 \
    libpango-1.0-0 \
    libpangocairo-1.0-0 \
    gdk-pixbuf2.0-0 \
    libffi-dev \
    shared-mime-info \
    fonts-liberation \
    fontconfig \
    && rm -rf /var/lib/apt/lists/*

# 3. Copy application files
COPY . .

# 4. Copy Nginx config (前回作成した「偽造キー付き」の設定ファイルをそのまま使います)
COPY nginx.conf /etc/nginx/nginx.conf

# 5. Register gbdraw bundled fonts
RUN mkdir -p /usr/share/fonts/truetype/gbdraw \
    && cp gbdraw/data/*.ttf /usr/share/fonts/truetype/gbdraw/ 2>/dev/null || true \
    && fc-cache -f -v

# 6. Install Python libraries
RUN pip3 install --no-cache-dir -r requirements.txt

# 7. Streamlit Config (GA tag insertion)
RUN sed -i 's~<head>~<head><script async src="https://www.googletagmanager.com/gtag/js?id=G-GG6JMKM02Y"></script><script>window.dataLayer = window.dataLayer || [];function gtag(){dataLayer.push(arguments);}gtag("js", new Date());gtag("config", "G-GG6JMKM02Y");</script>~' /usr/local/lib/python3.13/site-packages/streamlit/static/index.html

# 8. Expose Cloud Run port
EXPOSE 8080

# 9. Start Script
# Nginxを起動しつつ、Streamlitも起動
CMD service nginx start && \
    streamlit run app.py \
    --server.port=8501 \
    --server.address=0.0.0.0 \
    --server.enableCORS=false \
    --server.enableXsrfProtection=false \
    --server.maxUploadSize=1024
    