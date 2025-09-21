# Dockerfile for web interface using nginx
FROM nginx:alpine

# Copy web files to nginx default directory
COPY web/ /usr/share/nginx/html/

# Explicitly copy the build info file
COPY web/build-info.json /usr/share/nginx/html/build-info.json

# Copy nginx configuration
COPY web/nginx.conf /etc/nginx/conf.d/default.conf

# Expose port 80
EXPOSE 80

# Start nginx
CMD ["nginx", "-g", "daemon off;"]
