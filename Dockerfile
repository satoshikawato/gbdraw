# Use Python 3.13 as the base image
FROM python:3.13-slim

ENV PYTHONUNBUFFERED=1
WORKDIR /app

# 1. Install required libraries + Nginx
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

# 2. Copy application files
COPY . .

# 3. Copy Nginx config
# 作成した nginx.conf を所定の位置にコピー
COPY nginx.conf /etc/nginx/nginx.conf

# 4. Register gbdraw bundled fonts
RUN mkdir -p /usr/share/fonts/truetype/gbdraw \
    && cp gbdraw/data/*.ttf /usr/share/fonts/truetype/gbdraw/ 2>/dev/null || true \
    && fc-cache -f -v

# 5. Install Python libraries
RUN pip3 install --no-cache-dir -r requirements.txt

# 6. Streamlit Config (GA tag insertion)
RUN sed -i 's~<head>~<head><script async src="https://www.googletagmanager.com/gtag/js?id=G-GG6JMKM02Y"></script><script>window.dataLayer = window.dataLayer || [];function gtag(){dataLayer.push(arguments);}gtag("js", new Date());gtag("config", "G-GG6JMKM02Y");</script>~' /usr/local/lib/python3.13/site-packages/streamlit/static/index.html

# 7. Expose Cloud Run port
EXPOSE 8080

# 8. Start Script
# Nginxをバックグラウンドで起動し、Streamlitをフォアグラウンドで起動
# Streamlitは 8501 で待ち受ける（Nginxがそこへ転送する）
CMD service nginx start && \
    streamlit run app.py \
    --server.port=8501 \
    --server.address=0.0.0.0 \
    --server.enableCORS=false \
    --server.enableXsrfProtection=false \
    --server.maxUploadSize=500
    