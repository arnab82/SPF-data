using QCBase
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using Printf

@load "data_cmf_49_new.jld2"

M = 100


init_fspace =  [(5, 0), (4, 4), (0, 5), (12, 12), (0, 0)]
init_fspace = FermiCG.FockConfig(init_fspace)
@load "data_spt_08e4.jld2"
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

nroots=6


bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,max_iter=20, thresh_var=1e-2, thresh_foi=0.8e-4,thresh_pt=1e-3,thresh_spin=8e-4,ci_conv=4e-5,ci_max_iter=500,do_pt=true);
@save "/home/arnab22/FermiCG-data/bimetallics/fe2_morokuma/rohf.631gd/49/49_new/data_spt_0804_1e3.jld2" clusters init_fspace ints cluster_bases v
bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,max_iter=20, thresh_var=1e-2, thresh_foi=0.8e-4,thresh_pt=1e-3,thresh_spin=7e-4,ci_conv=4e-5,ci_max_iter=500,do_pt=true);
@save "/home/arnab22/FermiCG-data/bimetallics/fe2_morokuma/rohf.631gd/49/49_new/data_spt_07041e3.jld2" clusters init_fspace ints cluster_bases v
bst_energy, v = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,max_iter=20, thresh_var=1e-2, thresh_foi=0.5e-4,thresh_pt=1e-3,thresh_spin=6e-4,ci_conv=4e-5,ci_max_iter=500,do_pt=true);
@save "/home/arnab22/FermiCG-data/bimetallics/fe2_morokuma/rohf.631gd/49/49_new/data_spt_06041e3.jld2" clusters init_fspace ints cluster_bases v

