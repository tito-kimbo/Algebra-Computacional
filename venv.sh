#!/bin/bash
virtualenv venv -p `which python3`
. venv/bin/activate
pip install .
