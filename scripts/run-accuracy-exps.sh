# Script for experiments in the section 6.1 and 6.2

seed=0
now="$(date +'%Y%m%d%H%M')"
result_dir=../data/result/accuracy/${now}

methods=( wj cset cs jsub impr bsk sumrdf )
datasets=( lubm80 human aids yago )

p=0.03

for method in ${methods[@]}; do
  if [ ${method} == cset ] || [ ${method} == bsk ] || [ ${method} == sumrdf ] # deterministic estimator
  then
    repeat=1
  else
    repeat=30
  fi
  if [ ${method} == bsk ] # Bound sketch requires one more parameter: the budget of BSK
  then
    export GCARE_BSK_BUDGET=4096
  fi
  for data in ${datasets[@]}; do
    ./run-exp.sh ${method} ${data} ${p} ${seed} ${repeat} ${result_dir}
  done
done
