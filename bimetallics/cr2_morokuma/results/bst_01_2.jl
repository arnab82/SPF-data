using QCBase
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2


@load "data_cmf_29_cr2.jld2"
@load "data_bst_5e3.jld2"
M = 100 

init_fspace = FockConfig([(6, 3), (4,4), (3,6)])

#cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3,3], init_fspace, max_roots=M, verbose=1);

clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters)
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

nroots=4


ebst, vbst = FermiCG.block_sparse_tucker(vbst, cluster_ops, clustered_ham,max_iter=50, thresh_var=4e-3, thresh_foi=1e-5, thresh_pt=4e-4,thresh_spin=4e-3)
@save "/home/arnab22/SPF-data/cr2_morokuma/29__3d4s4d_2s2p3p_3d4s4d/spin_corrected/data_bst_4e3.jld2" clusters init_fspace ints cluster_bases vbst
e2_bst = FermiCG.compute_pt2_energy(vbst, cluster_ops, clustered_ham, thresh_foi=1e-6);
ebst, vbst = FermiCG.block_sparse_tucker(vbst, cluster_ops, clustered_ham,max_iter=500, thresh_var=3e-3, thresh_foi=1e-5, thresh_pt=2e-4,thresh_spin=3e-3)
@save "/home/arnab22/SPF-data/cr2_morokuma/29__3d4s4d_2s2p3p_3d4s4d/spin_corrected/data_bst_3e3.jld2" clusters init_fspace ints cluster_bases  vbst
e2_bst = FermiCG.compute_pt2_energy(vbst, cluster_ops, clustered_ham, thresh_foi=1e-6);
ebst, vbst = FermiCG.block_sparse_tucker(vbst, cluster_ops, clustered_ham,max_iter=500, thresh_var=2e-3, thresh_foi=1e-5, thresh_pt=1e-4,thresh_spin=2e-3)
@save "/home/arnab22/SPF-data/cr2_morokuma/29__3d4s4d_2s2p3p_3d4s4d/spin_corrected/data_bst_2e3.jld2" clusters init_fspace ints cluster_bases  vbst
e2_bst = FermiCG.compute_pt2_energy(vbst, cluster_ops, clustered_ham, thresh_foi=1e-6);
ebst, vbst = FermiCG.block_sparse_tucker(vbst, cluster_ops, clustered_ham,max_iter=500, thresh_var=1e-3, thresh_foi=1e-5, thresh_pt=0.8e-4,thresh_spin=2e-3)
@save "/home/arnab22/SPF-data/cr2_morokuma/29__3d4s4d_2s2p3p_3d4s4d/spin_corrected/data_bst_1e3.jld2" clusters init_fspace ints cluster_bases  vbst
e2_bst = FermiCG.compute_pt2_energy(vbst, cluster_ops, clustered_ham, thresh_foi=1e-6);
