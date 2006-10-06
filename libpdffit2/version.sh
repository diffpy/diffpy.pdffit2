#!/bin/sh

# print defines for DATE and VERSION macros

MYDIR=`dirname "$0"`

svn info "$MYDIR" 2>/dev/null | awk 'BEGIN { q = "\\\""; }
    /^Last Changed Rev: / && !lcrev++ { print "-DVERSION=" q "2.0." $4 q }
    /^Last Changed Date: / && !lcdate++ { print "-DDATE=" q $4 q }'
