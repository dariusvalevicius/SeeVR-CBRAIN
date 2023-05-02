FROM dvalev/seevr-base:latest

USER root

# Copy application
COPY ./standalone/ /opt/seevr/

# Set execute permissions
RUN chmod -R a+x /opt/seevr/

# Set working directory
WORKDIR /opt/seevr/

USER nonroot
