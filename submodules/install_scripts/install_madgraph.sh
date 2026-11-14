#!/bin/sh

set -eu

url="https://launchpad.net/mg5amcnlo/3.0/3.8.x/+download/MG5_aMC_v3.5.17.tar.gz"
script_dir=$(CDPATH= cd "$(dirname "$0")" && pwd)
install_dir="$script_dir/../mg5amcnlo-v3.5.17"
archive="$install_dir/MG5_aMC_v3.5.17.tar.gz"

mkdir -p "$install_dir"
curl -fL "$url" -o "$archive"
tar -xzf "$archive" --strip-components=1 -C "$install_dir"
rm "$archive"
