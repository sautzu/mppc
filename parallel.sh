#!/bin/bash
#並列に実行する
#例:sh parallel.sh 4 vis.mac
OUTPUT="../result/test.txt"

if [ $# -ne 2 ]; then
  echo "引数を正しく入力してください"
  echo "引数として並列にする数とマクロが必要です"
  exit 1
fi

seq -w $1 | xargs -I{} rm ../result/test{}.txt > /dev/null
seq -w $1 | xargs -P$1 -I{} sh child.sh {} $2
seq -w $1 | xargs -I{}  cat ../result/test{}.txt >> ${OUTPUT}
seq -w $1 | xargs -I{} rm ../result/test{}.txt
