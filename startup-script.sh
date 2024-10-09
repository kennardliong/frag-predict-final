#!/bin/bash

wget https://github.com/git-lfs/git-lfs/releases/download/v3.5.1/git-lfs-linux-amd64-v3.5.1.tar.gz
tar xvf git-lfs-linux-amd64-v3.5.1.tar.gz

cd git-lfs-3.5.1/
./install.sh

cd ..

git lfs fetch --all
git lfs checkout
git pull