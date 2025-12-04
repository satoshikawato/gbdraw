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
    nginx \
    && rm -rf /var/lib/apt/lists/*

# 3. Copy application files
COPY . .
# Fix line endings and set execute permission for run.sh
RUN sed -i 's/\r$//' run.sh && chmod +x run.sh

# 4. Important: Register gbdraw bundled fonts to the system
RUN mkdir -p /usr/share/fonts/truetype/gbdraw \
    && cp gbdraw/data/*.ttf /usr/share/fonts/truetype/gbdraw/ 2>/dev/null || true \
    && fc-cache -f -v

# 5. Install Python libraries
RUN pip3 install --no-cache-dir -r requirements.txt

# 6. Add Google Analytics to Streamlit
RUN python3 -c "import streamlit; import os; \
    p = os.path.join(os.path.dirname(streamlit.__file__), 'static', 'index.html'); \
    content = open(p).read(); \
    ga_code = '<head><script async src=\"https://www.googletagmanager.com/gtag/js?id=G-GG6JMKM02Y\"></script><script>window.dataLayer = window.dataLayer || [];function gtag(){dataLayer.push(arguments);}gtag(\"js\", new Date());gtag(\"config\", \"G-GG6JMKM02Y\");</script>'; \
    new_content = content.replace('<head>', ga_code); \
    open(p, 'w').write(new_content);"

# Set Streamlit configuration via environment variables
ENV STREAMLIT_SERVER_CORS_ALLOWED_ORIGINS='["https://gbdraw.app", "http://localhost:8501"]'

# 7. Expose port
EXPOSE 8080

# 8. Healthcheck
HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health || exit 1

# 9. Entry point command
ENTRYPOINT ["./run.sh"]