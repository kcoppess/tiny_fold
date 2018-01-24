#!/bin/bash

python /Users/kcoppess/yggdrasil/gprof2dot/gprof2dot.py -f pstats debug.statout > debug.dot
dot debug.dot -Tpng > debug.png
open debug.png
