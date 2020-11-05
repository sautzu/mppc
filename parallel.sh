#!/bin/bash
if [ $# -ne 2 ]; then
  echo "引数を正しく入力してください"
  echo "引数として並列にする数とマクロが必要です"
  exit 1
fi

seq -w $1 | xargs -P$1 -I{} ./LXe ../result/test{}.txt $2
