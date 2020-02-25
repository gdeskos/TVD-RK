#!/bin/bash

# Example configuration script for TIOGA Utilities

COMPILER=clang

EXTRA_ARGS="$@"

cmake \
    $EXTRA_ARGS \
    ..
