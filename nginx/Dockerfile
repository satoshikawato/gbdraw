# Dockerfile.nginx
FROM nginx:alpine

# 前回の設定ファイルを nginx.conf としてコピー
# (前回のテキストを nginx.conf という名前で保存しておく必要があります)
COPY nginx.conf /etc/nginx/nginx.conf

# 8080ポートを開放
EXPOSE 8080

# Nginxをフォアグラウンドで実行
CMD ["nginx", "-g", "daemon off;"]
