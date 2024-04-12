using QCBase
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using Printf

@load "data_bst.jld2"
display(clusters)
display(init_fspace)

clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.compute_pt2_energy(v, cluster_ops, clustered_ham, thresh_foi=1e-5)

FermiCG.compute_pt2_energy(v, cluster_ops, clustered_ham, thresh_foi=1e-6)

FermiCG.compute_pt2_energy(v, cluster_ops, clustered_ham, thresh_foi=1e-7)

FermiCG.compute_pt2_energy(v, cluster_ops, clustered_ham, thresh_foi=1e-8)
