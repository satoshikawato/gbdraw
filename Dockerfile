FROM python:3.13-slim

ENV PYTHONUNBUFFERED=1
WORKDIR /app

# Nginx 等のインストール
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

# フォント設定
RUN mkdir -p /usr/share/fonts/truetype/gbdraw \
    && cp gbdraw/data/*.ttf /usr/share/fonts/truetype/gbdraw/ 2>/dev/null || true \
    && fc-cache -f -v

RUN pip3 install --no-cache-dir -r requirements.txt

# GAタグ埋め込み
RUN sed -i 's~<head>~<head><script async src="https://www.googletagmanager.com/gtag/js?id=G-GG6JMKM02Y"></script><script>window.dataLayer = window.dataLayer || [];function gtag(){dataLayer.push(arguments);}gtag("js", new Date());gtag("config", "G-GG6JMKM02Y");</script>~' /usr/local/lib/python3.13/site-packages/streamlit/static/index.html

# 設定ファイルの配置
COPY nginx.conf /etc/nginx/nginx.conf
COPY run.sh /app/run.sh
RUN chmod +x /app/run.sh

EXPOSE 8080

# 【重要】ヘルスチェックは Nginx(8080) ではなく Streamlit(8501) を直接叩く
# NginxはHTTP/2専用にしているため、curlでのチェックが失敗するのを防ぐためです
HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health || exit 1

ENTRYPOINT ["/app/run.sh"]

