using QCBase
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using Printf


M = 100 
@load "data_bst_3e34.jld2"

clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);

nroots=4
bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,max_iter=500, thresh_var=8e-3, thresh_foi=1e-5,thresh_pt=1e-3,thresh_spin=6e-4,do_pt=false);
@save "data_bst_8e35.jld2" clusters init_fspace ints cluster_bases cluster_ops v 
bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,max_iter=500, thresh_var=8e-3, thresh_foi=5e-6,thresh_pt=1e-3,thresh_spin=5e-4,do_pt=false);
@save "data_bst_8e356.jld2" clusters init_fspace ints cluster_bases cluster_ops v
bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,max_iter=500, thresh_var=8e-3, thresh_foi=5e-6,thresh_pt=1e-3,thresh_spin=4e-4,do_pt=false);
@save "data_bst.jld2" clusters init_fspace ints cluster_bases cluster_ops v
