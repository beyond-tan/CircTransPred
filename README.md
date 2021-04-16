# CircTransPred
CircTransPred是一款识别环状RNA翻译潜能的工具。首先我们提取了一些常用的序列编码特征（开放阅读框、Fickett分数和六聚体频率）。然后使用改进的Tri-training方法来预测环状RNA的编码能力。

环境依赖：
1、Python3
2、Numpy
3、sklearn
4、joblib

使用方法：
python CircTransPred.py –input_file data/test.fa –output_file data/pred_result.txt
--input_file 指定输入文件（输入文件需为标准的FASTA格式文件）
--output_file 指定输出文件目录

输出文件格式如下：
circRNA-id	prediction-label	prediction-probability
circ1	0	0.1983548
circ2	1	0.8192133
