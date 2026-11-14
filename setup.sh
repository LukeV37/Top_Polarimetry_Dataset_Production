#!/bin/bash
WORK_DIR=$(pwd)
if [ ! -f ./submodules/venv/bin/activate ]; then
    cd submodules
    python3 -m venv venv
    cd venv
    source ./bin/activate
    pip install --upgrade pip
    pip install -r pip_requirements.txt
else
    source ./submodules/venv/bin/activate
fi
cd $WORK_DIR
