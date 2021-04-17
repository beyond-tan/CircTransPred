# CircTransPred

CircTransPred is a tool for identifying the protein-coding potential of circRNAs. First, we extracted some common sequence coding features(ORF, Fickett scores, and hexamer frequency). Then used an improved Tri-training method to predict the coding ability of circRNAs.

Requirements:
1. Python3
2. Numpy
3. sklearn
4. joblib

Usage:  
python CircTransPred.py --input_file data/test.fa --output_file data/pred_result.txt  
--input_file: input file must be standard FASTA format file  
--output_file: the output file used to store the prediction label and probability of input data  

Output file:
(circRNA-ID prediciton-label prediciton-probability)  
circ1 0 0.1983548  
circ2 1 0.8192133  
