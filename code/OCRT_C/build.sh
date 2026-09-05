#!/bin/sh
set -e
cd "$(dirname "$0")"
make -j2 all
echo BUILD_OK build/ocrt
