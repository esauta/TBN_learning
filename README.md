# License
Copyright 2020 Code by **Elisabetta Sauta**

"A Bayesian data fusion based approach for learning genome-wide transcriptional regulatory networks" is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

"A Bayesian data fusion based approach for learning genome-wide transcriptional regulatory networks" is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with "A Bayesian data fusion based approach for learning genome-wide transcriptional regulatory networks"

If not, see http://www.gnu.org/licenses/.

# TBN_learning
_TBN_learning is a MATLAB implemented hybrid algorithm for learning the structure of transcriptional Bayesian networks on a genome-wide scale. It jointly integrates ChIP-seq derived transcriptional regulations as prior information and gene expression data both as evidence and as sampling probability tied to each transcriptional dependency. Thanks to its search and score schema applied at local and global network level, the algorithm efficiently handles genome-scale networks_

# Cite
If you use this code in your research work, please cite:

Sauta, E., Demartini, A., Vitali, F. et al. A Bayesian data fusion based approach for learning genome-wide transcriptional regulatory networks. BMC Bioinformatics 21, 219 (2020). https://doi.org/10.1186/s12859-020-3510-1

# Installation 

**Recommended Software**:
A functional Matlab installation (MATLAB 2015b or higher).

## **Main scripts**:

* `Master.m`: It represents the main script to perform the hybrid structure learning on your Transcriptional Bayesian Network.
* `conf.m`: It is a configuration script in which some parameters (e.g. threshold as end criterion) and variables (e.g. number of algorithm global iterations) can be set.
* `Initial_learning.m`: A starting script to prepare input data for the learning phase, estimating the initial scores of the network.

* `Complementary functions(*.m)`: all required functions to allow the learning on the network model at local and global level.

## Folder organization

* **Main Folder**:
The main folder contains the algorithm main script, its configuration script and all required functions. The user can define a sub-directory to save his/her results, setting the full path in the configuration script.

* **Input folder**:
All input data should be created as tab separated files and stored in this sub-directory of the main folder. An example of dataset is here available representing how valid input data should be created and to test how the algorithm works.
The structure of example data should be mantained in the user input data. 

## Example Data
The transcriptional Bayesian model used as algorithm's input is directed acyclic graph (DAG).

* _All_TFTF_edges.txt_: All weighted interactions among core regulators TFs. If interactions are not weighted, set an equal weight (e.g. 0.5) to all edges.

* _whitelist.txt_: All excluded TF-TF arcs from the initial DAG definition, constitued of only TF-TF interactions.

* _edges_TFTF_initialDAG_binding.txt_: TF-TF weighted interactions which define the initial DAG. If interactions are not weighted, set an equal weight (e.g. 0.5) to all edges.

* _TF_genes_binding.txt_: All weighted interactions from core TFs to their target genes. If interactions are not weighted, set an equal weight (e.g. 0.5) to all edges.

* _Expression_matrix.txt_: Gene Expression Matrix containing expression values for all nodes of the initial transcriptional Bayesian Network; samples in columns, genes/TFs names on rows.
This experimental dataset is created from this normalized data [Spellman PT, Sherlock G, Zhang MQ, Iyer VR, Anders K, Eisen MB, et al. Comprehensive Identification of Cell Cycle regulated Genes of the Yeast Saccharomyces cerevisiae by Microarray Hybridization. Mol Biol Cell 1998;9:3273-97].

* _Starting_wksp.mat_: The initial workspace with processed input data created by Initial_learning.m script; if learning variable (in the conf.m script) is set to zero, Master.m will load this workspace.

## Run the algorithm on the example data

* In the conf.m script:

  --The specified threshold works for example data

  -- n_run can be set to 2 (1:2) to quickly see how the algorithm works and what it returns as output
  
* In the Master.m script:

  -- gs_prov is initialized to work on example data; it must be lower that the global network score (gs_old) calculated on the starting Bayesian Network during the initial learning phase

  -- numb_to_extract variable at lines 59, 64 and 86 can be uncommented to quickly run the algorithm on example data.
  
## Results

The code will save a global workspace at the end of each iteration, and the final one (globalwkspace_*_exiting.mat) is created when the algorithm ends its run, in which the final Bayesian model is saved and can be exported.
Local and global network parameters are saved into a report file. All model changes during the algorithm run are reported in a log file. 



### Contributor

Elisabetta Sauta, BMS Lab, University of Pavia, via Ferrata 5 (Italy).
# # 
