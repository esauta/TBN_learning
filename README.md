# TBN_learning
TBN_learning is a MATLAB implemented hybrid algorithm for learning the structure of transcriptional Bayesian networks on a genome-wide scale. It exploits ChIP-seq derived transcriptional regulations as prior information and gene expression data during the learning phase. Thanks to its search and score schema applied at local and global network level, it can handle genomic transcriptional models. 

**** Installation *****

Recommended Software:

A functional Matlab installation (MATLAB 2015b or higher).

Main scripts:

<Master.m>: It represents the main script to perform the hybrid structure learning on your Transcriptional Bayesian Network.
<conf.m>: It is a configuration script in which some parameters (e.g. threshold as end criterion) and variables (e.g. number of algorithm global iterations) can be set.
<Initial_learning.m>: A starting script to prepare input data for the learning phase, estimating the initial scores of the network.

<Complementary functions (*.m)>: all required functions to allow the learning on the network model at local and global level.

***** Folder organization *****

Main Folder
The main folder contains the algorithm main script, its configuration script and all required functions. The user can define a sub-directory to save his/her results, setting the full path in the configuration script.

Input folder
All input data should be created as tab separated files and stored in this sub-directory of the main folder. An example of dataset is here available representing how valid input data should be created and to test how the algorithm works.
The structure of example data should be mantained in the user input data. 

***** Example Data *****

<All_TFTF_edges.txt>: All weighted interactions among core regulators TFs. If interactions are not weighted, set an equal weight (e.g. 0.5) to all edges.

<whitelist.txt>: All excluded TF-TF arcs from the initial DAG definition, constitued of only TF-TF interactions.

<edges_TFTF_initialDAG_binding.txt>: TF-TF weighted interactions which define the initial DAG. If interactions are not weighted, set an equal weight (e.g. 0.5) to all edges.

<TF_genes_binding.txt>: All weighted interactions from core TFs to their target genes. If interactions are not weighted, set an equal weight (e.g. 0.5) to all edges.

<Expression_matrix.txt>: Gene Expression Matrix containing expression values for all nodes of the initial transcriptional Bayesian Network; samples in columns, genes/TFs names on rows.
This experimental dataset is created from this normalized data [Spellman PT, Sherlock G, Zhang MQ, Iyer VR, Anders K, Eisen MB, et al. Comprehensive Identification of Cell Cycle regulated Genes of the Yeast Saccharomyces cerevisiae by Microarray Hybridization. Mol Biol Cell 1998;9:3273-97].

<Starting_wksp.mat>: The initial workspace with processed input data created by Initial_learning.m script; if learning variable (in the conf.m script) is set to zero, Master.m will load this workspace.

***** Run the algorithm on the example data *******

- In the conf.m script:
  -- the specified threshold works for example data;
  -- n_run can be set to 2 (1:2) to quickly see how the algorithm works and what it returns as output.
  
- In the Master.m script:
  -- gs_prov is initialized to work on example data; it must be lower that the global network score (gs_old) calculated on the starting Bayesian Network during the initial learning phase
  -- numb_to_extract variable at lines 59, 64 and 86 can be uncommented to quickly run the algorithm on example data.
  
***** Results *****

The code will save a global workspace at the end of each iteration, and the final one (globalwkspace_*_exiting.mat) is created when the algorithm ends its run, in which the final Bayesian model is saved and can be exported.
Local and global network parameters are saved into a report file. All model changes during the algorithm run are reported in a log file. 



Contributor

Elisabetta Sauta, BMS Lab, University of Pavia, via Ferrata 5 (Italy).
