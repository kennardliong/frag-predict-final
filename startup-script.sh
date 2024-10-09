#!/bin/bash

apt update
apt install git-lfs
git lfs install

git lfs fetch --all
git lfs checkout
git pull