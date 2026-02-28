#!/usr/bin/env bash

set -e

source .venv/bin/activate
export PYTHONPATH="$(pwd)"
python src/main.py
