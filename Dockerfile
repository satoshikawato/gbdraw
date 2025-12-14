# 1. Use Python 3.13 as the base image
FROM python:3.13-slim

ENV PYTHONUNBUFFERED=1
WORKDIR /app

# 2. Install required libraries 
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
    && rm -rf /var/lib/apt/lists/*

# 3. Copy application files
COPY . .

# 4. Important: Register gbdraw bundled fonts to the system
RUN mkdir -p /usr/share/fonts/truetype/gbdraw \
    && cp gbdraw/data/*.ttf /usr/share/fonts/truetype/gbdraw/ 2>/dev/null || true \
    && fc-cache -f -v

# 5. Install Python libraries
RUN pip3 install --no-cache-dir -r requirements.txt

# 6. Add Google Analytics to Streamlit
RUN sed -i 's~<head>~<head><script async src="https://www.googletagmanager.com/gtag/js?id=G-GG6JMKM02Y"></script><script>window.dataLayer = window.dataLayer || [];function gtag(){dataLayer.push(arguments);}gtag("js", new Date());gtag("config", "G-GG6JMKM02Y");</script>~' /usr/local/lib/python3.13/site-packages/streamlit/static/index.html

# 7. Expose port (削除または8501に変更)
EXPOSE 8501

# 8. Healthcheck (ポートを8501に変更)
HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health || exit 1

# 9. Entry point command
ENTRYPOINT ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0", "--server.enableCORS=false", "--server.enableXsrfProtection=false", "--server.enableWebsocketCompression=false"]