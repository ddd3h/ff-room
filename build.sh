#!/usr/bin/env bash
# Build the C++ extension module and place it in the Python package.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="$ROOT/cpp/build"
PYTHON_PKG="$ROOT/python/ff_room"

# Resolve Python: PYTHON_EXECUTABLE env var > conda > active venv > system.
if [ -n "${PYTHON_EXECUTABLE:-}" ]; then
    PYTHON_EXE="$PYTHON_EXECUTABLE"
elif [ -n "${CONDA_PREFIX:-}" ] && [ -f "$CONDA_PREFIX/bin/python3" ]; then
    PYTHON_EXE="$CONDA_PREFIX/bin/python3"
elif [ -f "$HOME/miniconda3/bin/python3" ]; then
    PYTHON_EXE="$HOME/miniconda3/bin/python3"
elif [ -f "$HOME/anaconda3/bin/python3" ]; then
    PYTHON_EXE="$HOME/anaconda3/bin/python3"
elif [ -f "$HOME/miniforge3/bin/python3" ]; then
    PYTHON_EXE="$HOME/miniforge3/bin/python3"
else
    PYTHON_EXE="$(which python3 || which python)"
fi

echo "=== ff-room build ==="
echo "Root:     $ROOT"
echo "Build:    $BUILD_DIR"
echo "Target:   $PYTHON_PKG"
echo "Python:   $PYTHON_EXE"
echo

mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

cmake "$ROOT/cpp" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PYTHON_PKG" \
    -DPython3_EXECUTABLE="$PYTHON_EXE"

make -j"$(nproc 2>/dev/null || sysctl -n hw.logicalcpu)"
make install

echo
echo "=== Build done ==="
echo "Extension installed to: $PYTHON_PKG"
ls "$PYTHON_PKG"/_ffroom_core*.so 2>/dev/null && echo "OK" || echo "WARNING: .so not found"
