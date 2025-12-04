FROM python:3.13-slim

ENV PYTHONUNBUFFERED=1
WORKDIR /app

# Install Nginx (official) & dependencies
RUN apt-get update && apt-get install -y \
    curl gnupg2 ca-certificates lsb-release debian-archive-keyring \
    supervisor \
    git \
    libcairo2 libpango-1.0-0 libpangocairo-1.0-0 \
    gdk-pixbuf2.0-0 libffi-dev shared-mime-info \
    fonts-liberation fontconfig \
    && curl https://nginx.org/keys/nginx_signing.key | gpg --dearmor \
       | tee /usr/share/keyrings/nginx-archive-keyring.gpg >/dev/null \
    && echo "deb [signed-by=/usr/share/keyrings/nginx-archive-keyring.gpg] \
       http://nginx.org/packages/debian bookworm nginx" \
       | tee /etc/apt/sources.list.d/nginx.list \
    && apt-get update && apt-get install -y nginx \
    && rm -rf /var/lib/apt/lists/*

COPY . .

# Register fonts
RUN mkdir -p /usr/share/fonts/truetype/gbdraw \
    && cp gbdraw/data/*.ttf /usr/share/fonts/truetype/gbdraw/ 2>/dev/null || true \
    && fc-cache -f -v

# Install Python libs
RUN pip install --no-cache-dir -r requirements.txt

# Google Analytics inject
RUN sed -i 's~<head>~<head><script async src="https://www.googletagmanager.com/gtag/js?id=G-GG6JMKM02Y"></script><script>window.dataLayer = window.dataLayer || [];function gtag(){dataLayer.push(arguments);}gtag("js", new Date());gtag("config", "G-GG6JMKM02Y");</script>~' \
    /usr/local/lib/python3.13/site-packages/streamlit/static/index.html

# Copy nginx and supervisor configs
COPY nginx.conf /etc/nginx/nginx.conf
COPY supervisord.conf /etc/supervisor/conf.d/supervisord.conf

# --- ADD THIS LINE ---
# Copy Streamlit config to disable CORS/XSRF
COPY config.toml /root/.streamlit/config.toml
# ---------------------

ENV PORT=8080
EXPOSE 8080

CMD ["/usr/bin/supervisord", "-n"]

