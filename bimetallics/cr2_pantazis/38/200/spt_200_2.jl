using QCBase
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using Printf

@load "data_cmf.jld2"

M = 200
display(clusters)
display(init_fspace)
init_fspace = FermiCG.FockConfig(init_fspace)
#cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3,3,3,3], init_fspace, max_roots=M, verbose=1);
@load "data_bst_200.jld2"
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

nroots=4


bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,max_iter=5, thresh_var=7e-3, thresh_foi=5e-5,thresh_pt=1e-3,thresh_spin=9e-4);
@save "/home/arnab22/SPF-data/cr2_pantazis/rohf.def2svp/38__3d4d_2p3p/spin_corrected/data_bst_200.jld2" clusters init_fspace ints cluster_bases v
bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,max_iter=5, thresh_var=7e-3, thresh_foi=5e-5,thresh_pt=1e-3,thresh_spin=9e-4);
@save "/home/arnab22/SPF-data/cr2_pantazis/rohf.def2svp/38__3d4d_2p3p/spin_corrected/data_bst_200.jld2" clusters init_fspace ints cluster_bases v
bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,max_iter=5, thresh_var=7e-3, thresh_foi=4e-5,thresh_pt=1e-3,thresh_spin=9e-4);
@save "/home/arnab22/SPF-data/cr2_pantazis/rohf.def2svp/38__3d4d_2p3p/spin_corrected/data_bst_200.jld2" clusters init_fspace ints cluster_bases v
bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,max_iter=5, thresh_var=7e-3, thresh_foi=3e-5,thresh_pt=1e-3,thresh_spin=9e-4);
@save "/home/arnab22/SPF-data/cr2_pantazis/rohf.def2svp/38__3d4d_2p3p/spin_corrected/data_bst_200.jld2" clusters init_fspace ints cluster_bases v
bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,max_iter=5, thresh_var=7e-3, thresh_foi=2e-5,thresh_pt=1e-3,thresh_spin=9e-4);
@save "/home/arnab22/SPF-data/cr2_pantazis/rohf.def2svp/38__3d4d_2p3p/spin_corrected/data_bst_200.jld2" clusters init_fspace ints cluster_bases v
bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,max_iter=2, thresh_var=7e-3, thresh_foi=1e-5,thresh_pt=1e-3,thresh_spin=9e-4);
@save "/home/arnab22/SPF-data/cr2_pantazis/rohf.def2svp/38__3d4d_2p3p/spin_corrected/data_bst_200.jld2" clusters init_fspace ints cluster_bases v
bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,max_iter=50, thresh_var=7e-3, thresh_foi=1e-5,thresh_pt=1e-3,thresh_spin=9e-4);
@save "/home/arnab22/SPF-data/cr2_pantazis/rohf.def2svp/38__3d4d_2p3p/spin_corrected/data_bst_200.jld2" clusters init_fspace ints cluster_bases v
