#!/bin/bash
# This script is a shorthand method of running configure with
# a specific set of options, including a custom version string
# which indicates the current git branch name and commit count.
# Copy this file to ./config and edit the options as required.
# This allows you to specify your own site-specific prefix and
# other preferences in a file not tracked by git.
./configure \
--prefix=/usr/local/spinhawk \
--enable-setuid-hercifc=hercules \
--enable-multi-cpu=32 \
--enable-custom=`awk -F"[(,)]" '/^AM_INIT_AUTOMAKE/ {print $3}' configure.ac`\
-spinhawk-`git rev-parse --abbrev-ref HEAD`-`git log --oneline | wc -l`
