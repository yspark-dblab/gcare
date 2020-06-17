# Script for experiments in the section 6.3

seed=0
repeat=30
now="$(date +'%Y%m%d%H%M')"
result_dir=../data/result/varying-sampling-ratios/${now}

for method in wj cs impr jsub; do
  for data in aids yago; do
    for p in 0.0001 0.0003 0.001 0.003 0.01 0.03; do
      ./run-exp.sh ${method} ${data} ${p} ${seed} ${repeat} ${result_dir}
    done
  done
done
