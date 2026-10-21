#!/usr/bin/env bash
# Build the C++ extension module and place it in the Python package.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="$ROOT/cpp/build"
PYTHON_PKG="$ROOT/python/ff_room"

# Resolve Python: explicit override > active venv > conda > PATH.
if [ -n "${PYTHON_EXECUTABLE:-}" ]; then
    PYTHON_EXE="$PYTHON_EXECUTABLE"
elif [ -n "${VIRTUAL_ENV:-}" ] && [ -x "$VIRTUAL_ENV/bin/python" ]; then
    PYTHON_EXE="$VIRTUAL_ENV/bin/python"
elif [ -n "${CONDA_PREFIX:-}" ] && [ -f "$CONDA_PREFIX/bin/python3" ]; then
    PYTHON_EXE="$CONDA_PREFIX/bin/python3"
else
    PYTHON_EXE="$(command -v python3 || command -v python)"
fi

echo "=== ff-room build ==="
echo "Root:     $ROOT"
echo "Build:    $BUILD_DIR"
echo "Target:   $PYTHON_PKG"
echo "Python:   $PYTHON_EXE"
echo

cmake -S "$ROOT/cpp" -B "$BUILD_DIR" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PYTHON_PKG" \
    -DPython3_EXECUTABLE="$PYTHON_EXE"

cmake --build "$BUILD_DIR" --config Release --parallel "${CMAKE_BUILD_PARALLEL_LEVEL:-2}"
cmake --install "$BUILD_DIR" --config Release --component python

echo
echo "=== Build done ==="
echo "Extension installed to: $PYTHON_PKG"
PYTHONPATH="$ROOT/python${PYTHONPATH:+:$PYTHONPATH}" "$PYTHON_EXE" -c \
    'from ff_room import _ffroom_core; print("OK:", _ffroom_core.__file__)'
