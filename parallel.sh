#!/bin/bash
if [ $# -ne 2 ]; then
  echo "引数を正しく入力してください"
  echo "引数として並列にする数とマクロが必要です"
  exit 1
fi

seq -w $1 | xargs -P$1 -I{} ./LXe ../result/test{}.txt $2
seq -w $1 | xargs -I{}  cat ../result/test{}.txt >> ../result/test.txt
seq -w $1 | xargs -I{} rm ../result/test{}.text
