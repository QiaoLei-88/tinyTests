#!/bin/bash
FILENAME="output"
if [ -n "$1" ]; then
  FILENAME=$1
fi

darwin=false
case "`uname`" in
  Darwin*) darwin=true ;;
esac

if $darwin ; then
  sedi="sed -i ''"
else
  sedi="sed -i"
fi

eval "$sedi -e 's/JobId.*$//g' $FILENAME"
