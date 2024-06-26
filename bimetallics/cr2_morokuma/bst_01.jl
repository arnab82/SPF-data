using QCBase
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2


@load "data_cmf_29_cr2.jld2"

M = 100 

init_fspace = FockConfig([(6, 3), (4,4), (3,6)])

cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3,3], init_fspace, max_roots=M, verbose=1);

clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters)
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

nroots=4

# start by defining P/Q spaces
p_spaces = Vector{ClusterSubspace}()

ssi = ClusterSubspace(clusters[1])
add_subspace!(ssi, (6,3), 1:1)
add_subspace!(ssi, (5,4), 1:1)
add_subspace!(ssi, (4,5), 1:1)
add_subspace!(ssi, (3,6), 1:1)
push!(p_spaces, ssi)

ssi = ClusterSubspace(clusters[2])
add_subspace!(ssi, (4,4), 1:1)
push!(p_spaces, ssi)

ssi = ClusterSubspace(clusters[3])
add_subspace!(ssi, (6,3), 1:1)
add_subspace!(ssi, (5,4), 1:1)
add_subspace!(ssi, (4,5), 1:1)
add_subspace!(ssi, (3,6), 1:1)
push!(p_spaces, ssi)


ci_vector = BSTstate(clusters, p_spaces, cluster_bases, R=4) 

na = 13
nb = 13
FermiCG.fill_p_space!(ci_vector, na, nb)
FermiCG.eye!(ci_vector)
ebst, vbst = FermiCG.ci_solve(ci_vector, cluster_ops, clustered_ham)

ebst, vbst = FermiCG.block_sparse_tucker(vbst, cluster_ops, clustered_ham,max_iter=50, thresh_var=1e-2, thresh_foi=1e-4, thresh_pt=1e-3,thresh_spin=1e-2)
@save "/home/arnab22/SPF-data/cr2_morokuma/29__3d4s4d_2s2p3p_3d4s4d/spin_corrected/data_bst_e2.jld2" clusters init_fspace ints cluster_bases  vbst
ebst, vbst = FermiCG.block_sparse_tucker(vbst, cluster_ops, clustered_ham,max_iter=50, thresh_var=8e-3, thresh_foi=5e-5, thresh_pt=1e-3,thresh_spin=8e-3)
@save "/home/arnab22/SPF-data/cr2_morokuma/29__3d4s4d_2s2p3p_3d4s4d/spin_corrected/data_bst_8e3.jld2" clusters init_fspace ints cluster_bases vbst

ebst, vbst = FermiCG.block_sparse_tucker(vbst, cluster_ops, clustered_ham,max_iter=50, thresh_var=6e-3, thresh_foi=1e-5, thresh_pt=8e-4,thresh_spin=8e-3)
@save "/home/arnab22/SPF-data/cr2_morokuma/29__3d4s4d_2s2p3p_3d4s4d/spin_corrected/data_bst_6e3.jld2" clusters init_fspace ints cluster_bases vbst
e2_bst = FermiCG.compute_pt2_energy(vbst, cluster_ops, clustered_ham, thresh_foi=1e-6);
ebst, vbst = FermiCG.block_sparse_tucker(vbst, cluster_ops, clustered_ham,max_iter=500, thresh_var=5e-3, thresh_foi=1e-5, thresh_pt=6e-4,thresh_spin=6e-3)
@save "/home/arnab22/SPF-data/cr2_morokuma/29__3d4s4d_2s2p3p_3d4s4d/spin_corrected/data_bst_5e3.jld2" clusters init_fspace ints cluster_bases  vbst
e2_bst = FermiCG.compute_pt2_energy(vbst, cluster_ops, clustered_ham, thresh_foi=1e-6);
