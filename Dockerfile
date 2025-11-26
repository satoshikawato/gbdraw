# 1. Use Python 3.13
FROM python:3.13-slim

# 2. Set unbuffered output for logs
ENV PYTHONUNBUFFERED=1

# 3. Set working directory
WORKDIR /app

# 4. Install necessary OS libraries
# These are necessary OS libraries for cairosvg, etc.
RUN apt-get update && apt-get install -y \
    git \
    libcairo2 \
    libpango-1.0-0 \
    libpangocairo-1.0-0 \
    gdk-pixbuf2.0-0 \
    libffi-dev \
    shared-mime-info \
    curl \
    && rm -rf /var/lib/apt/lists/*

# 5. Copy project files
COPY . .

# 6. Install Python libraries
RUN pip3 install --no-cache-dir -r requirements.txt

# 7. Expose port
EXPOSE 8080

# 8. Health check
HEALTHCHECK CMD curl --fail http://localhost:8080/_stcore/health || exit 1

# 9. Set entrypoint command
ENTRYPOINT ["streamlit", "run", "app.py", "--server.port=8080", "--server.address=0.0.0.0"]
