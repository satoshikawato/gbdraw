# Dockerfile for deploying the application with Nginx on Cloud Run
FROM python:3.13-slim

ENV PYTHONUNBUFFERED=1
WORKDIR /app

# Install necessary system packages including Nginx
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

# Install custom fonts for gbdraw
RUN mkdir -p /usr/share/fonts/truetype/gbdraw \
    && cp gbdraw/data/*.ttf /usr/share/fonts/truetype/gbdraw/ 2>/dev/null || true \
    && fc-cache -f -v

# Copy application files including nginx.conf and run.sh
COPY . .

# Give execute permission to run.sh
RUN chmod +x run.sh

# Install Python libraries
RUN pip3 install --no-cache-dir -r requirements.txt

# Add Google Analytics to Streamlit
RUN python3 -c "import streamlit; import os; \
    p = os.path.join(os.path.dirname(streamlit.__file__), 'static', 'index.html'); \
    content = open(p).read(); \
    ga_code = '<head><script async src=\"https://www.googletagmanager.com/gtag/js?id=G-GG6JMKM02Y\"></script><script>window.dataLayer = window.dataLayer || [];function gtag(){dataLayer.push(arguments);}gtag(\"js\", new Date());gtag(\"config\", \"G-GG6JMKM02Y\");</script>'; \
    new_content = content.replace('<head>', ga_code); \
    open(p, 'w').write(new_content);"

# Set environment variables for Streamlit configuration
ENV STREAMLIT_SERVER_CORS_ALLOWED_ORIGINS='["https://gbdraw.app", "http://localhost:8501"]'

# Expose port 8080 for Cloud Run
EXPOSE 8080

# Add health check for Cloud Run
HEALTHCHECK CMD curl --fail http://localhost:8080/_stcore/health || exit 1

# Set entrypoint to run.sh
ENTRYPOINT ["./run.sh"]