seed=0
commit_id=$(git rev-parse HEAD)
now="$(date +'%Y%m%d%H%M')"
result_dir=../data/result/build/${now}/

mkdir -p ${result_dir}

p=0.03

for data in lubm80 human aids yago; do
  input_dir=../data/dataset/${data}/
  data_dir=../data/dataset/${data}/

  for method in wj cs cset impr jsub sumrdf bsk; do

    if [ ${method} == cs ] || [ ${method} == bsk ]
    then
      bin=../build/gcare_relation
    else
      bin=../build/gcare_graph
    fi

    if [ ${method} == bsk ] # Bound sketch requires one more parameter: the budget of BSK
    then
      export GCARE_BSK_BUDGET=4096
    fi

    if [ ${method} == sumrdf ]  # SumRDF requires one more parameter: the threshold to determine whether merging two summary vertices
    then
      export GCARE_SUMRDF_THRESHOLD=0.15
    fi

    output_file=${result_dir}/${data}_${method}_p${p}_s${seed}_build_result.txt

    cmd="${bin} -b -m ${method} -i ${input_dir}/${data}.txt -d ${data_dir}/${data} -p ${p} -s ${seed} -o ${output_file}"
    echo ${cmd}
    ${cmd}

    echo "#${commit_id}" >> ${output_file}
  done
done
