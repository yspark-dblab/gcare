# Script for experiments in the section 6.4

mode=$1 # it has 2 modes: 1) b - build mode and 2) q - query mode
commit_id=$(git rev-parse HEAD)
now="$(date +'%Y%m%d%H%M')"

datasets=( lubm80 aids )
p=0.03

if [ ${mode} == b ] # this is a mode to measure summary build times
then

  methods=( cset bsk sumrdf )
  result_dir=../data/result/build/${now}/
  mkdir -p ${result_dir}

  for method in ${methods[@]}; do
    if [ ${method} == cs ] || [ ${method} == bsk ]
    then
      bin=../build/gcare_relation
    else
      bin=../build/gcare_graph
    fi

    if [ ${method} == sumrdf ] # SumRDF requires one more parameter: the threshold to determine whether merging two summary vertices
    then
      export GCARE_SUMRDF_THRESHOLD=0.15
    fi
    if [ ${method} == bsk ] # Bound sketch requires one more parameter: the budget of BSK
    then
      export GCARE_BSK_BUDGET=4096
    fi

    for data in ${datasets[@]}; do
      input_dir=../data/dataset/${data}/
      data_dir=../data/dataset/${data}/
      for seed in 0 1 2 3 4 5; do
        if [ ${method} == bsk ]
        then
          output_file=${result_dir}/${data}_${method}_b${GCARE_BSK_BUDGET}_s${seed}_build_result.txt
        else
          output_file=${result_dir}/${data}_${method}_p${p}_s${seed}_build_result.txt
        fi
        cmd="${bin} -b -m ${method} -i ${input_dir}/${data}.txt -d ${data_dir}/${data} -p ${p} -s ${seed} -o ${output_file}"
        echo "${cmd}"
        ${cmd}
        echo "#${commit_id}" >> ${output_file}
      done
    done
  done
  exit
fi
\
if [ ${mode} == q ] # this is a mode to measure query times
then
  methods=( wj cset cs impr jsub bsk sumrdf )
  seed=0
  repeat=30
  result_dir=../data/result/efficiency/${now}

  for method in ${methods[@]}; do
    if [ ${method} == bsk ]
    then
      export GCARE_BSK_BUDGET=4096
    fi
    for data in ${datasets[@]}; do
      ./run-exp.sh ${method} ${data} ${p} ${seed} ${repeat} ${result_dir}
    done
  done
fi
