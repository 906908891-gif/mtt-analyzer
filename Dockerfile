FROM python:3.11-slim

WORKDIR /app

# Install dependencies first for better caching
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
# Scipy is optional but recommended for better 4PL fit
RUN pip install --no-cache-dir scipy

# Copy application code
COPY . .

# Environment
ENV HOST=0.0.0.0
ENV PORT=5000
ENV FLASK_APP=app.py
ENV PYTHONUNBUFFERED=1

EXPOSE 5000

# Healthcheck
HEALTHCHECK --interval=30s --timeout=5s --start-period=10s \
    CMD python -c "import urllib.request; urllib.request.urlopen("http://localhost:5000/api/health").read()" || exit 1

# Run with gunicorn in production mode
CMD ["gunicorn", "--bind", "0.0.0.0:5000", "--workers", "2", "--timeout", "120", "app:create_app()"]
