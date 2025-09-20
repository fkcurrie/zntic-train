# Dockerfile for web interface using nginx
FROM nginx:alpine

# Copy web files to nginx default directory
COPY web/ /usr/share/nginx/html/

# Copy nginx configuration
COPY nginx.conf /etc/nginx/conf.d/default.conf

# Expose port 80
EXPOSE 80

# Start nginx
CMD ["nginx", "-g", "daemon off;"]