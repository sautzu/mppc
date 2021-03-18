#!/bin/bash
#乱数のシード値を時間から取っているので、並列に実行する際は遅延を挟む必要がある
sleep ${1}0s
./LXe ../result/test${1}.txt $2
