using QCBase
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using Printf

@load "data_bst_4e44.jld2"
display(clusters)
display(init_fspace)

clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);
FermiCG.compute_pt2_energy(v, cluster_ops, clustered_ham, thresh_foi=5e-5)
FermiCG.compute_pt2_energy(v, cluster_ops, clustered_ham, thresh_foi=1e-5)

FermiCG.compute_pt2_energy(v, cluster_ops, clustered_ham, thresh_foi=1e-6)

FermiCG.compute_pt2_energy(v, cluster_ops, clustered_ham, thresh_foi=1e-7)

FermiCG.compute_pt2_energy(v, cluster_ops, clustered_ham, thresh_foi=1e-8)
