# Builder stage
FROM rust:1.82-slim-bookworm AS builder

# Install build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /build

# Clone the repository
RUN git clone https://github.com/bguo068/hmmibd-rs.git .

# Build hmmibd-rs (Rust version)
RUN cargo build --release

# Build hmmIBD (C version)
WORKDIR /build/c
RUN cc -O3 -Wall hmmIBD.c -o hmmIBD -lm

# Final stage - minimal runtime image
FROM debian:bookworm-slim

# Install only runtime dependencies (if any needed)
# For these binaries, we might only need basic C libraries which are already in debian-slim
RUN apt-get update && apt-get install -y \
    libgcc-s1 \
    libc6 \
    && rm -rf /var/lib/apt/lists/*

# Copy only the compiled binaries from builder stage
COPY --from=builder /build/target/release/hmmibd-rs /usr/local/bin/
COPY --from=builder /build/c/hmmIBD /usr/local/bin/

# Set working directory for runtime
WORKDIR /data

# Verify binaries are executable
RUN chmod +x /usr/local/bin/hmmibd-rs /usr/local/bin/hmmIBD

# Set default command to show help
CMD ["hmmibd-rs"]
