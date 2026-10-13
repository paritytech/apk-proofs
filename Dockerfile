# syntax=docker/dockerfile:1
#
# Builds the crate and `docker run` runs the full test suite and both
# light-client simulations — APK-377 (BLS12-377/BW6-761) and APK-381 (BLS12-381/BW6-767).
#
#   docker build -t apk-proofs .
#   docker run --rm apk-proofs                 # tests + both light clients, defaults
#   docker run --rm -e VALIDATORS=1500 -e ERAS=3 apk-proofs
#   docker run --rm apk-proofs tests           # tests only
#   docker run --rm apk-proofs light-clients   # light clients only

FROM rust:1.96-slim-bookworm

WORKDIR /apk-proofs
COPY . .

# Compile everything up front so `docker run` measures proving, not compiling. `--locked` holds
# the build to Cargo.lock.
RUN cargo build --locked --release -p apk-proofs --features print-trace --examples \
 && cargo test  --locked --release -p apk-proofs --no-run

COPY --chmod=755 <<'EOF' /usr/local/bin/apk-proofs-run
#!/bin/sh
set -e
VALIDATORS="${VALIDATORS:-1000}"
ERAS="${ERAS:-2}"

run_tests() {
  echo "==> Unit, integration and doc tests (both configurations)"
  cargo test --locked --release -p apk-proofs
  echo "==> Rader's algorithm against a naive DFT at 10177 and 20354 points"
  cargo test --locked --release -p apk-proofs --lib -- --ignored rader
}

run_light_clients() {
  for example in apk_377 apk_381; do
    echo "==> Light client: ${example}, ${VALIDATORS} validators, ${ERAS} eras"
    cargo run --locked --release -q -p apk-proofs --features print-trace \
      --example "${example}" -- "${VALIDATORS}" "${ERAS}"
  done
}

case "${1:-all}" in
  tests)         run_tests ;;
  light-clients) run_light_clients ;;
  all)           run_tests; run_light_clients ;;
  *)             exec "$@" ;;
esac
EOF

ENTRYPOINT ["apk-proofs-run"]
CMD ["all"]
