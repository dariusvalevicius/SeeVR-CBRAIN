FROM dvalev/seevr-base:latest

USER root

# Install dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
  && apt-get install -y dc

# Copy application
COPY ./standalone/ /opt/seevr/
COPY preprocess_BIDS.sh /opt/seevr/preprocess_BIDS.sh

# Set execute permissions
#RUN chmod -R a+x /opt/seevr/
RUN chown -R nonroot: /opt/seevr

# Set working directory
WORKDIR /opt/seevr

USER nonroot
