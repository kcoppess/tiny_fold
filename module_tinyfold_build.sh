#!/bin/bash

cmake -G Ninja
ninja
mv libtinyfold.so tinyfold.so
