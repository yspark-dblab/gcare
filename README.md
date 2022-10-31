

# G-CARE


## How to install and run G-CARE

1. Download the repository
    ```
    git clone https://github.com/yspark-dblab/gcare.git
    ```
2. Compile the framework
    * Compiler: g++ version 7.3.1 20180303 (Red Hat 7.3.1-5) (GCC)
    * Libraries: boost 1.58.0
    * CMake: minimum 2.8
    
    ```
    mkdir build
    cd build 
    cmake ..
    make -j  
    ```
    
3. Run G-CARE
    1. Convert a raw input graph to a binary format graph and build a summary of estimator
        * Command to build a binary format data graph and a summary
        
            ```
            ./gcare_graph -b -m ${method} -i ${input_data_graph}.txt -d ${output_data_graph} -o ${output_file}
            ./gcare_relation -b -m ${method} -i ${input_data_graph}.txt -d ${output_data_graph} -o ${output_file}
            ```
        * Example for the toy data and toy query
        
            ```
            ./gcare_graph -b -m wj -i ../data/dataset/toy/toy.txt -d ../data/dataset/toy/toy -o ./toy_wj_build_result.txt
            ./gcare_relation -b -m bsk -i ../data/dataset/toy/toy.txt -d ../data/dataset/toy/toy -o ./toy_bsk_build_result.txt
            ```
            
    2. Run an estimator
        * Command to run estimator for given queries, a data_graph, and a summary
        
            ```
            ./gcare_graph -q -m ${method} -i ${input_dir} -d ${data_dir}/${data} -o ${output_file}            
            ./gcare_relation -q -m ${method} -i ${input_dir} -d ${data_dir}/${data} -o ${output_file}
            ```
        * Example for the toy data and toy query
        
            ```
            ./gcare_graph -q -m wj -i ../data/queryset/toy -d ../data/dataset/toy/toy -o ./toy_wj_query_result.txt
            ./gcare_relation -q -m bsk -i ../data/queryset/toy -d ../data/dataset/toy/toy -o ./toy_bsk_query_result.txt
            ```

    3. You can run the experiments of the paper [1] using following scripts. You also need to download the dataset and queryset. See "Dataset and Queryset" section of this document.
        * ./scripts/run-build.sh : a script to build binary-format data graphs and summaries of estimators
        * ./scripts/run-accuracy-exps.sh : a script for the section 6.1 and 6.2
        * ./scripts/run-varying-sampling-ratios-exps.sh : a script for the section 6.3
        * ./scripts/run-efficiency-exps.sh : a script for the section 6.4

## Output file

1. Stdout and a query result file 
    * File name: ${output_file}
    * For each query, the framework prints 2 lines to the file and stdout:
        * The 1st line: {query_file} {the average of estimated cardinality} {the average of elapsed time} 
        {the variance of estimated cardinality}
        * The 2nd line: {the number of repetitions} {the 1st estimated cardinality}, 
        {the 2nd estimated cardinality} ... {the last estimated cardinality}
    * The following is an example
    
        ```
        ../queryset/wj/aids/Tree_9/uf_Q_4_3.txt 63373.639640640926 28.860840633333325 13339858036.646763
        30 4709.8505140316756 4709.8505140316756 246645.60600166713 55959.877743817728 5052.4167824395663 393795.06307307584 4528.7024173381496 5062.6551819949982 4343.4465129202554 43089.49819394276 4902.1083634342876 4750.2439566546263 453 8.2250625173656 14012.968046679633 4300.3145318143925 53451.118644067799 4285.0969713809391 47396.281189219226 4759.8599611003056 5097.8827452070018 4487.9666574048342 4457.8427340928038 458384.97360377881 42835.934426229505 198388.1     833842734 15053.510419560989 199227.79438732981 53734.464017782717 4506.1405946096138 4741.3125868296747
        ```
2. A query error file
    * File name: ${output_file}.err
    * For each failed query, the framework prints a line
        * {query_file} error with code {error_code}
        * {query_file} error with signal {signal}
    * Error code (you can add your error code)
        * Time Out Error - 1
    * Signal: https://www.man7.org/linux/man-pages/man7/signal.7.html 
        
    * The following is an example
    
        ```
        ../queryset/sumrdf/yago_vm/Tree_3/f_Q_0_6.txt error with code 1
        ../queryset/sumrdf/yago_vm/Tree_3/f_Q_0_6.txt error with signal 6
        ```

3. A build result file
    * File name: ${output_file} 
    * The framework prints the elapsed time of a summary build.
    
     
## Dataset and Queryset

1. Query graph
    * The query graph file is in text format describing a query graph with M vertices, N edges.

    * The first line is the header and should be t # s \<query_id\> where \<query_id\> 
    is an arbitrary integer which is no smaller than 0.
    
    * Next N lines represent N query vertices sorted by query vertex id. 
    Each line should be v \<id\> \<label\> \<dvid\>, where \<id\> 
    is the query vertex id starting from 0, 
    \<label\> is the query vertex label, and \<dvid\> is the 
    data vertex id of a data vertex bound to this query vertex.
    The default value of \<label\> and \<dvid\> is -1.
    The meaning of the default value is as follows.
    A query vertex with \<label\>=-1 can be matched with any data vertices regardless of their data vertex labels.
    A query vertex with \<dvid\>=-1 can be matched with any data vertices regardless of their data vertex ids.
    

    
    * Next M lines represent M query edges. 
    Each line should be e \<id1\> \<id2\> \<label\>, representing a directed, 
    labeled edge from query vertex with id <id1> to the query vertex with 
    id \<id2\> whose label is \<label\>.
    
    * The query graph files used in the paper's experiments are in querysets directory.
        - LUBM80 [2] - Files in querysets/lubm80/ directory: Figure 6 and 10
        - HUMAN [3] - Files in querysets/human/ directory: Figure 7 and 8
        - AIDS [4] - Files in querysets/aids/ directory: Figure 7, 8, 9, and 10
        - YAGO [5] - Files in querysets/yago/ directory: Figure 6
    
    * The following is a typical example.
        ```
        t # s 54
        v 0 2 -1
        v 1 -1 -1
        v 2 -1 -1
        v 3 -1 5
        e 1 0 33
        e 2 1 19
        e 1 3 0
        ```

2. Data graph
    * The data graph file is in text format describing a labeled/unlabeled graph with M vertices, N edges.
    
    * The first line is the header and should be t # \<graph_id\> where \<graph_id\> 
    is an arbitrary integer which is no smaller than 0.
    
    * Next N lines represent N data vertices in increasing order of data vertex id. 
    Each line should be v \<id\> \<label1\> \<label2\> …  where \<id\> is the data vertex id 
    and \<label1\>, \<label2\>, … are data vertex labels 
    (i.e., a data vertex can have multiple labels). The default value of this set of labels is 0.
    
    * Next M lines represent M data edges and have the same format as the query edges.
    
    * The data graph files used in the paper's experiments are in datasets directory.
        - LUBM80 [2] - datasets/lubm80/lubm80.txt : Figure 6 and 10
        - HUMAN [3] - datasets/human/human.txt : Figure 7 and 8
        - AIDS [4] - datasets/aids/aids.txt : Figure 7, 8, 9, and 10
        - YAGO [5] - datasets/yago/yago.txt : Figure 6

        For RDF datasets, the data graph files are the results of type-aware transformation [6], thus each label represents a certain type.

    * The following is a typical example.
       ```
       t # 16
       v 0 0
       v 1 0
       v 2 1 2 5
       v 3 2 7
       v 4 9
       e 2 3 3
       e 2 1 19
       e 1 0 0
       e 2 4 5
       ```
    
3. How to download them?
    * Download files from the following links using Chrome
        * [queryset.tar.gz](https://drive.google.com/file/d/1Dlj43rBAOVPAsfzKlYxIbZ9RsqeGM_MN/view?usp=sharing)
        consists of query graphs for LUBM80 [2], HUMAN [3], AIDS [4], and YAGO [5]
            * The path of each query has the following meaning.  
                * LUBM80 - ${data}/lubm80_Q${query_id}.txt
                * HUMAN - ${data}/${topology}\_${topology_size}/uf\_Q\_${<img src="https://render.githubusercontent.com/render/math?math=\lfloor \log(max(1,true\_cardinality-1)) \rfloor">}\_${query_id}.txt
                * AIDS - ${data}/${topology}\_${topology_size}/uf\_Q\_${<img src="https://render.githubusercontent.com/render/math?math=\lfloor \log(max(1,true\_cardinality-1)) \rfloor">}\_${query_id}.txt
                * YAGO - ${data}/${topology}\_${topology_size}/f\_Q\_${<img src="https://render.githubusercontent.com/render/math?math=\lfloor \log(max(1,true\_cardinality-1)) \rfloor">}\_${query_id}.txt
                
            * You can also download the true cardinalities of the queries from the [true_cardinality.tar.gz](https://drive.google.com/file/d/1Bc6Q2RZQTcIB8IfOw5KafNYwPhq2BO94/view?usp=sharing) using Chrome
        * [dataset.tar.gz](https://drive.google.com/file/d/1HAgSVE-24NOap6_Q1_twH56Dkb2kPvGU/view?usp=sharing)
        consists of data graphs for LUBM80 [2], HUMAN [3], AIDS [4], and YAGO [5] 
    * Decompress the files under the ./data/dataset and ./data/queryset
    
        ```
        cd ./data/queryset
        # download the queryset.tar.gz from the link
        tar -zxvf queryset.tar.gz
        cd ./data/dataset
        # download the dataset.tar.gz from the link
        tar -zxvf dataset.tar.gz
        ```
## Figures of Experimental Results
We follow the next procedures to report the accuracy of techniques in the paper.
    1) We run each query 30 times and compute q-errors for each estimate in a set of 30 estimates.
    2) We obtain a distribution consisting of 30 * (# of queries in each group on the x-axis) q-error values. For example, in Figure 6(c), the q-error distribution for star queries consists of 30 * 318 q-error values since the number of star queries for YAGO dataset is 318.
    3) From this distribution, we report the average and standard deviation for LUBM and report the 5%, 25%, 50%, 75%, and 95% percentiles for the other datasets in Figures 6-8.
    
In Figures 6-8, we differentiate over- and under-estimation cases for better explanation of estimation behaviors. However, for the same query and estimator, some out of 30 trials may fall in under-estimation cases while others fall in over-estimation cases. Therefore, averaging the 30 q-errors for each query may hard to answer whether an estimator under-estimates for this query. In that sense, for each query group Q, we plot 5%, 25%, 50%, 75%, and 95% percentiles of 30 * |Q| trials in Figures 6-8.
    
If one does not want to differentiate between over- and under-estimation cases but just compare q-error between queries (i.e., how large is the average relative error on a specific query?), one could aggregate the q-errors over the trials.

## Errata
* For LUBM, CharacteristicSets, SumRDF, and BoundSketch spend 0.96, 12.26, and <a><img src="https://drive.google.com/uc?export=view&id=16bziI4RL0NWtB7JrylVo2t8cLtQrgqpN"></a><a><img src="https://drive.google.com/uc?export=view&id=1uQ3D5nZ1o_c5FbGpX9-NpBXtaVphK8xB"></a>seconds, respectively.
* For AIDS, CharacteristicSets, SumRDF, and BoundSketch spend 0.07, 0.64, and <a><img src="https://drive.google.com/uc?export=view&id=1X1T3kXnnDdSnCoW2zt5U49tvgBl8CLOx"></a><a><img src="https://drive.google.com/uc?export=view&id=1Rq-dSewWqaOKy49qjgTTIJ8QpZrFVDmy"></a>seconds, respectively.
<br>(In Section 6.4, Paragraph 2)
  
## Contact us
If you have any inquiry or bug report, please send emails to me at yspark@dblab.postech.ac.kr.

## References
[1] Yeonsu Park,  Seongyun Ko,  Sourav S. Bhowmick,  Kyoungmin Kim,  Kijae Hong,  Wook-Shin Han. 2020. 
G-CARE: A Framework for Performance Benchmarking of Cardinality Estimation Techniques for Subgraph Matching.
In Proceedings of the 2020 International Conference on Management of Data. 1099-1114.

[2] Yuanbo Guo, Zhengxiang Pan, and Jeff Heflin. 2005. LUBM: A benchmark for OWL knowledge base systems. 
    Journal of Web Semantics 3, 2-3 (2005), 158–182.
    
[3] Shijie Zhang, Shirong Li, and Jiong Yang. 2009. GADDI: distance index
     based subgraph matching in biological networks. In Proceedings of
     the 12th International Conference on Extending Database Technology:
     Advances in Database Technology. 192–203.
     
[4] Haichuan Shang, Ying Zhang, Xuemin Lin, et al. 2008. 
     Taming verification hardness: an efficient algorithm for testing subgraph isomorphism.
     Proceedings of the VLDB Endowment 1, 1 (2008), 364–375.

[5] Fabian M Suchanek, Gjergji Kasneci, and Gerhard Weikum. 2008. Yago:
     A large ontology from wikipedia and wordnet. Journal of Web Semantics 6, 3 (2008), 203–217.   

[6] Jinha Kim, Hyungyu Shin, Wook-Shin Han, Sungpack Hong, Hassan Chafi. 2015. Taming Subgraph Isomorphism for RDF Query Processing. Proceedings of the VLDB Endowment, 8(11).
