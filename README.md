# AFP2-RNA-Seq
RNA-Seq data and analysis of an Arabidopsis mutant overexpressing AFP2, a transcription factor that represses abscisic acid response by interacting with ABI5. 
Mutants and WT seeds are observed in three conditions: control media, exposure to abscisic acid (ABA), and exposure to trichostatin A (TSA).
Some of the replicates are large outliers, so medians were used for clustering.

All files are reproducible with main_script.R

The original hypothesis was that AFP2 represses ABI5 in part by recruiting histone deactylases. If HDACs are inhibited by trichostatin A, AFP2 repression would also be inhibited.
The original hypothesis has mostly been rejected, but there are some clusters that do fit the pattern 
(differential expression between Y-AFP2 and WT in GM and ABA, but not in TSA)

In general TSA has very broad effects, and the experiment might benefit from more specific HDAC inhibitors.

As a whole we saw what we would expect to see: 
ABA specifically regulates lipid storage, response to cold, cell homeostasis, and response to water. 
TSA (an HDAC inhibitor) upregulates genes related to chromatin organization, negative regulation of gene expression, and cell response to DNA damage stimulus.

Interesting findings:

AFP2 overexpression upregulates glutathione related genes, nutrient reservoir activity, and water transport. 

TSA significantly upregulates ribosome biogenesis genes, even though the seeds are not in a state of growth. (see KEGG maps)
This might be a compensatory mechanism for the increased transcription caused by inhibition of histone deacetylases.

This experiment is still in progress. See "results_of_interest" in /k12_cluster_analysis for a short list of GO term analysis.

