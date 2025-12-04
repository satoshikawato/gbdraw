#!/bin/bash

streamlit run app.py --server.port=8501 --server.address=0.0.0.0 --server.maxUploadSize=500 --server.maxMessageSize=1000 --server.enableXsrfProtection=false --server.enableCORS=true &

nginx -c /app/nginx.conf -g 'daemon off;'