method=$1
data=$2
p=$3
seed=$4
repeat=$5
result_dir=$6

echo ${method}

if [ ${method} == cs ] || [ ${method} == bsk ]
then
bin=../build/gcare_relation
else
bin=../build/gcare_graph
fi

input_dir=../data/queryset/${data}/
data_dir=../data/dataset/${data}/

if [ ${method} == bsk ]
then
  output_file=${result_dir}/${data}_${method}_b${GCARE_BSK_BUDGET}_s${seed}_query_result.txt
else
  output_file=${result_dir}/${data}_${method}_p${p}_s${seed}_query_result.txt
fi
commit_id=$(git rev-parse HEAD)

mkdir -p ${result_dir}

cmd="${bin} -q -m ${method} -i ${input_dir} -d ${data_dir}/${data} -p ${p} -n ${repeat} -s ${seed} -o ${output_file}"
echo ${cmd}
${cmd}

echo "#${commit_id}" >> ${output_file}
