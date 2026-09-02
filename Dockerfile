FROM gcc:14-bookworm AS usalign-builder

ARG USALIGN_SHA256="5e05ddcd66d7ca7553954b7df94e938fcfbf1a5d94a81edf5b04e80f8d8913f5"

WORKDIR /build
COPY vendor/USalign/USalign.cpp ./USalign.cpp
RUN set -eu; \
    echo "$USALIGN_SHA256  USalign.cpp" | sha256sum --check -; \
    g++ -static -O3 -ffast-math -o USalign USalign.cpp; \
    status=0; \
    ./USalign -v > /tmp/usalign-version.txt 2>&1 || status=$?; \
    cat /tmp/usalign-version.txt; \
    test "$status" -eq 1; \
    grep -q "Version 20260527" /tmp/usalign-version.txt

FROM python:3.12-slim-bookworm AS runtime

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    USALIGN_BIN=/usr/local/bin/USalign

WORKDIR /app
COPY pyproject.toml README.md LICENSE THIRD_PARTY_NOTICES.md ./
COPY src ./src
RUN pip install --no-cache-dir .

COPY --from=usalign-builder /build/USalign /usr/local/bin/USalign
COPY app.py ./app.py
COPY .streamlit ./.streamlit
COPY examples/input ./examples/input

RUN useradd --create-home --uid 10001 structdiff \
    && chown -R structdiff:structdiff /app
USER structdiff

EXPOSE 8501
HEALTHCHECK --interval=30s --timeout=5s --start-period=20s --retries=3 \
  CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8501/_stcore/health', timeout=3)" || exit 1

CMD ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]