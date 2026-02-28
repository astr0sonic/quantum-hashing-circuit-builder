#!/usr/bin/env bash

set -e

python3.12 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install -r requirements-dev.txt
pre-commit install
deactivate
