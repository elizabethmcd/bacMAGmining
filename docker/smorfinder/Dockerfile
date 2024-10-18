# python 3.7 full image
FROM python:3.7

# Set the working directory in the container
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    gcc \
    g++ \
    make \
    libpython3-dev \
    && rm -rf /var/lib/apt/lists/*

# Clone the SmORFinder repository
RUN git clone https://github.com/elizabethmcd/SmORFinder.git

# Install SmORFinder
RUN pip install ./SmORFinder

# Set the PATH to include SmORFinder scripts
ENV PATH="/app/SmORFinder/scripts:${PATH}"

# Command to run when starting the container
CMD ["smorf", "--help"]