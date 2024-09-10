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
cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3,3,3,3], init_fspace, max_roots=M, verbose=1);

clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

nroots=4
# start by defining P/Q spaces
p_spaces = Vector{ClusterSubspace}()
ssi = ClusterSubspace(clusters[1])
add_subspace!(ssi, (2,5), 1:1)
add_subspace!(ssi, (3,4), 1:1)
add_subspace!(ssi, (4,3), 1:1)
add_subspace!(ssi, (5,2), 1:1)
push!(p_spaces, ssi)


ssi = ClusterSubspace(clusters[2])
add_subspace!(ssi, (2,5), 1:1)
add_subspace!(ssi, (3,4), 1:1)
add_subspace!(ssi, (4,3), 1:1)
add_subspace!(ssi, (5,2), 1:1)
push!(p_spaces, ssi)
ssi = ClusterSubspace(clusters[3])
add_subspace!(ssi, init_fspace[3], 1:1)
push!(p_spaces, ssi)
ssi = ClusterSubspace(clusters[4])
add_subspace!(ssi, init_fspace[4], 1:1)
push!(p_spaces, ssi)
ssi = ClusterSubspace(clusters[4])
add_subspace!(ssi, init_fspace[4], 1:1)
push!(p_spaces, ssi)
    


ci_vector = BSTstate(clusters, p_spaces, cluster_bases, R=4) 

na = sum([i[1] for i in init_fspace]) 
nb = sum([i[2] for i in init_fspace]) 

FermiCG.fill_p_space!(ci_vector, na, nb)
FermiCG.eye!(ci_vector)
e_ci, v0 = FermiCG.ci_solve(ci_vector, cluster_ops, clustered_ham);


bst_energy, v = FermiCG.block_sparse_tucker(v0, cluster_ops, clustered_ham,max_iter=2, thresh_var=9e-3, thresh_foi=5e-5,thresh_pt=1e-3,thresh_spin=9e-4);
@save "/home/arnab22/SPF-data/cr2_pantazis/rohf.def2svp/38__3d4d_2p3p/spin_corrected/data_bst_200.jld2" clusters init_fspace ints cluster_bases v
bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,max_iter=2, thresh_var=8e-3, thresh_foi=5e-5,thresh_pt=1e-3,thresh_spin=9e-4);
@save "/home/arnab22/SPF-data/cr2_pantazis/rohf.def2svp/38__3d4d_2p3p/spin_corrected/data_bst_200.jld2" clusters init_fspace ints cluster_bases v
bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,max_iter=2, thresh_var=8e-3, thresh_foi=5e-5,thresh_pt=1e-3,thresh_spin=9e-4);
@save "/home/arnab22/SPF-data/cr2_pantazis/rohf.def2svp/38__3d4d_2p3p/spin_corrected/data_bst_200.jld2" clusters init_fspace ints cluster_bases v
bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,max_iter=2, thresh_var=8e-3, thresh_foi=5e-5,thresh_pt=1e-3,thresh_spin=9e-4);
@save "/home/arnab22/SPF-data/cr2_pantazis/rohf.def2svp/38__3d4d_2p3p/spin_corrected/data_bst_200.jld2" clusters init_fspace ints cluster_bases v
bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,max_iter=2, thresh_var=8e-3, thresh_foi=5e-5,thresh_pt=1e-3,thresh_spin=9e-4);
@save "/home/arnab22/SPF-data/cr2_pantazis/rohf.def2svp/38__3d4d_2p3p/spin_corrected/data_bst_200.jld2" clusters init_fspace ints cluster_bases v
bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,max_iter=2, thresh_var=8e-3, thresh_foi=5e-5,thresh_pt=1e-3,thresh_spin=9e-4);
@save "/home/arnab22/SPF-data/cr2_pantazis/rohf.def2svp/38__3d4d_2p3p/spin_corrected/data_bst_200.jld2" clusters init_fspace ints cluster_bases v
bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,max_iter=2, thresh_var=8e-3, thresh_foi=5e-5,thresh_pt=1e-3,thresh_spin=9e-4);
@save "/home/arnab22/SPF-data/cr2_pantazis/rohf.def2svp/38__3d4d_2p3p/spin_corrected/data_bst_200.jld2" clusters init_fspace ints cluster_bases v
bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,max_iter=2, thresh_var=8e-3, thresh_foi=5e-5,thresh_pt=1e-3,thresh_spin=9e-4);
@save "/home/arnab22/SPF-data/cr2_pantazis/rohf.def2svp/38__3d4d_2p3p/spin_corrected/data_bst_200.jld2" clusters init_fspace ints cluster_bases v
bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,max_iter=2, thresh_var=8e-3, thresh_foi=5e-5,thresh_pt=1e-3,thresh_spin=9e-4);
@save "/home/arnab22/SPF-data/cr2_pantazis/rohf.def2svp/38__3d4d_2p3p/spin_corrected/data_bst_200.jld2" clusters init_fspace ints cluster_bases v
bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,max_iter=2, thresh_var=8e-3, thresh_foi=5e-5,thresh_pt=1e-3,thresh_spin=9e-4);
@save "/home/arnab22/SPF-data/cr2_pantazis/rohf.def2svp/38__3d4d_2p3p/spin_corrected/data_bst_200.jld2" clusters init_fspace ints cluster_bases v
bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,max_iter=50, thresh_var=8e-3, thresh_foi=5e-5,thresh_pt=1e-3,thresh_spin=9e-4);
@save "/home/arnab22/SPF-data/cr2_pantazis/rohf.def2svp/38__3d4d_2p3p/spin_corrected/data_bst_200.jld2" clusters init_fspace ints cluster_bases v
