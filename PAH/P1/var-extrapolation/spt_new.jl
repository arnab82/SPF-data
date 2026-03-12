using FermiCG
using LinearAlgebra
using Printf
using JLD2
@load "data_benz_p1_cmf.jld2"
max_roots = 60
# Build Cluster basis
cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1_n, [3,3,3,3], FermiCG.FockConfig(init_fspace), max_roots=max_roots, verbose=1);

# Build ClusteredOperator
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);

# Build Cluster Operators
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);
nroots = 1

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1_n.a, d1_n.b)

v = FermiCG.BSTstate(clusters, FockConfig(init_fspace), cluster_bases)
thresh_spf=0.01
e_var1, v_var = FermiCG.block_sparse_tucker(v, cluster_ops, clustered_ham,
                                               max_iter    = 10,
                                               nbody       = 4,
                                               H0          = "Hcmf",
                                               thresh_var  = thresh_spf,
                                               thresh_foi  = thresh_spf/50,
                                               thresh_pt   = thresh_spf/2,
                                               ci_conv     = 5e-5,
                                               do_pt       = false,
                                               tol_tucker  = 1e-5,
                                               resolve_ss  = true,
                                               verbose     = 1)
@save "data_spt.jld2" e_var1 v_var cluster_bases

t1 = time()
alloc = @allocated begin
    σ = FermiCG.build_compressed_1st_order_state(v_var, cluster_ops, clustered_ham,
                                                 nbody=4, thresh=1e-9, compress_iter=false,compress_twice=false)
    σ = FermiCG.compress(σ, thresh=1e-9)
    H2 = FermiCG.orth_dot(σ, σ)
end
var1 = zeros(length(e_var1))
for r in 1:length(e_var1)
    var1[r] = H2[1] - e_var1[r] * e_var1[r]
end
@printf("Allocated memory: %.3f GB\n", alloc/1e9)
@printf("Variance, Energy: %1.5e, %1.8e\n", var1[1], e_var1[1])
t2= time()
println("Time taken: $(t2-t1) seconds")
t1=time()
alloc= @allocated σ2 = FermiCG.compute_spt_sigma_norm_blockwise(v_var, cluster_ops, clustered_ham;
                                       H0="Hcmf", nbody=4, thresh_foi=1e-8,
                                       max_number=nothing, opt_ref=true,
                                       ci_tol=1e-6, verbose=1)

var_like = σ2 .- e_var1.^2
println("Variance-like, Energy: $(var_like[1]), $(e_var1[1])")
@printf("Allocated memory: %.3f GB\n", alloc/1e9)
t2 = time()
println("Time taken: $(t2-t1) seconds")

t1=time()
alloc= @allocated σ2 = FermiCG.compute_spt_sigma_norm_blockwise_new(v_var, cluster_ops, clustered_ham;
                                       H0="Hcmf", nbody=4, thresh_foi=1e-8,
                                       max_number=nothing, opt_ref=true,
                                       ci_tol=1e-6, verbose=1)

var_like = σ2 .- e_var1.^2
println("Variance-like, Energy: $(var_like[1]), $(e_var1[1])")
@printf("Allocated memory: %.3f GB\n", alloc/1e9)
t2 = time()
println("Time taken: $(t2-t1) seconds")

thresh_spf=0.008
e_var2, v_var = FermiCG.block_sparse_tucker(v_var, cluster_ops, clustered_ham,
                                               max_iter    = 10,
                                               nbody       = 4,
                                               H0          = "Hcmf",
                                               thresh_var  = thresh_spf,
                                               thresh_foi  = thresh_spf/50,
                                               thresh_pt   = thresh_spf/2,
                                               ci_conv     = 5e-5,
                                               do_pt       = false,
                                               tol_tucker  = 1e-5,
                                               resolve_ss  = true,
                                               verbose     = 1)
@save "data_spt.jld2" e_var2 v_var cluster_bases
t1=time()
alloc= @allocated σ2 = FermiCG.compute_spt_sigma_norm_blockwise(v_var, cluster_ops, clustered_ham;
                                       H0="Hcmf", nbody=4, thresh_foi=1e-8,
                                       max_number=nothing, opt_ref=true,
                                       ci_tol=1e-6, verbose=1)

var_like2 = σ2 .- e_var2.^2
println("Variance-like, Energy: $(var_like2[1]), $(e_var2[1])")
@printf("Allocated memory: %.3f GB\n", alloc/1e9)
t2 = time()
println("Time taken: $(t2-t1) seconds")
t1=time()
alloc= @allocated σ2 = FermiCG.compute_spt_sigma_norm_blockwise_new(v_var, cluster_ops, clustered_ham;
                                       H0="Hcmf", nbody=4, thresh_foi=1e-8,
                                       max_number=nothing, opt_ref=true,
                                       ci_tol=1e-6, verbose=1)

var_like2 = σ2 .- e_var2.^2
println("Variance-like, Energy: $(var_like2[1]), $(e_var2[1])")
@printf("Allocated memory: %.3f GB\n", alloc/1e9)
t2 = time()
println("Time taken: $(t2-t1) seconds")

thresh_spf=0.006
e_var3, v_var = FermiCG.block_sparse_tucker(v_var, cluster_ops, clustered_ham,
                                               max_iter    = 20 ,
                                               nbody       = 4,
                                               H0          = "Hcmf",
                                               thresh_var  = thresh_spf,
                                               thresh_foi  = thresh_spf/50,
                                               thresh_pt   = thresh_spf/2,
                                               ci_conv     = 5e-5,
                                               do_pt       = false,
                                               tol_tucker  = 1e-5,
                                               resolve_ss  = true,
                                               verbose     = 1)
@save "data_spt.jld2" e_var3 v_var cluster_bases
t1=time()
alloc= @allocated σ2 = FermiCG.compute_spt_sigma_norm_blockwise(v_var, cluster_ops, clustered_ham;
                                       H0="Hcmf", nbody=4, thresh_foi=1e-8,
                                       max_number=nothing, opt_ref=true,
                                       ci_tol=1e-6, verbose=1)

var_like3 = σ2 .- e_var3.^2
println("Variance-like, Energy: $(var_like3[1]), $(e_var3[1])")
@printf("Allocated memory: %.3f GB\n", alloc/1e9)
t2 = time()
println("Time taken: $(t2-t1) seconds")

t1=time()
alloc= @allocated σ2 = FermiCG.compute_spt_sigma_norm_blockwise_new(v_var, cluster_ops, clustered_ham;
                                       H0="Hcmf", nbody=4, thresh_foi=1e-8,
                                       max_number=nothing, opt_ref=true,
                                       ci_tol=1e-6, verbose=1)

var_like3 = σ2 .- e_var3.^2
println("Variance-like, Energy: $(var_like3[1]), $(e_var3[1])")
@printf("Allocated memory: %.3f GB\n", alloc/1e9)
t2 = time()
println("Time taken: $(t2-t1) seconds")


thresh_spf=0.004
e_var4, v_var = FermiCG.block_sparse_tucker(v_var, cluster_ops, clustered_ham,
                                               max_iter    = 20,
                                               nbody       = 4,
                                               H0          = "Hcmf",
                                               thresh_var  = thresh_spf,
                                               thresh_foi  = thresh_spf/50,
                                               thresh_pt   = thresh_spf/2,
                                               ci_conv     = 5e-5,
                                               do_pt       = false,
                                               tol_tucker  = 1e-5,
                                               resolve_ss  = true,
                                               verbose     = 1)
@save "data_spt.jld2" e_var4 v_var cluster_bases

t1=time()
alloc= @allocated σ2 = FermiCG.compute_spt_sigma_norm_blockwise_new(v_var, cluster_ops, clustered_ham;
                                       H0="Hcmf", nbody=4, thresh_foi=1e-8,
                                       max_number=nothing, opt_ref=true,
                                       ci_tol=1e-6, verbose=1)

var_like4 = σ2 .- e_var4.^2
println("Variance-like, Energy: $(var_like4[1]), $(e_var4[1])")
@printf("Allocated memory: %.3f GB\n", alloc/1e9)
t2 = time()
println("Time taken: $(t2-t1) seconds")
t1=time()
alloc= @allocated σ2 = FermiCG.compute_spt_sigma_norm_blockwise_new(v_var, cluster_ops, clustered_ham;
                                       H0="Hcmf", nbody=4, thresh_foi=1e-8,
                                       max_number=nothing, opt_ref=true,
                                       ci_tol=1e-6, verbose=1)

var_like4 = σ2 .- e_var4.^2
println("Variance-like, Energy: $(var_like4[1]), $(e_var4[1])")
@printf("Allocated memory: %.3f GB\n", alloc/1e9)
t2 = time()
println("Time taken: $(t2-t1) seconds")
thresh_spf=0.002
e_var5, v_var = FermiCG.block_sparse_tucker(v_var, cluster_ops, clustered_ham,
                                               max_iter    = 20,
                                               nbody       = 4,
                                               H0          = "Hcmf",
                                               thresh_var  = thresh_spf,
                                               thresh_foi  = thresh_spf/50,
                                               thresh_pt   = thresh_spf/2,
                                               ci_conv     = 5e-5,
                                               do_pt       = false,
                                               tol_tucker  = 1e-5,
                                               resolve_ss  = true,
                                               verbose     = 1)
@save "data_spt.jld2" e_var5 v_var cluster_bases

t1=time()
alloc= @allocated σ2 = FermiCG.compute_spt_sigma_norm_blockwise(v_var, cluster_ops, clustered_ham;
                                       H0="Hcmf", nbody=4, thresh_foi=1e-8,
                                       max_number=nothing, opt_ref=true,
                                       ci_tol=1e-6, verbose=1)

var_like5 = σ2 .- e_var5.^2
println("Variance-like, Energy: $(var_like5[1]), $(e_var5[1])")
@printf("Allocated memory: %.3f GB\n", alloc/1e9)
t2 = time()
println("Time taken: $(t2-t1) seconds")
t1=time()
alloc= @allocated σ2 = FermiCG.compute_spt_sigma_norm_blockwise_new(v_var, cluster_ops, clustered_ham;
                                       H0="Hcmf", nbody=4, thresh_foi=1e-8,
                                       max_number=nothing, opt_ref=true,
                                       ci_tol=1e-6, verbose=1)

var_like5 = σ2 .- e_var5.^2
println("Variance-like, Energy: $(var_like5[1]), $(e_var5[1])")
@printf("Allocated memory: %.3f GB\n", alloc/1e9)
t2 = time()
println("Time taken: $(t2-t1) seconds")

thresh_spf=0.001
e_var7, v_var = FermiCG.block_sparse_tucker(v_var, cluster_ops, clustered_ham,
                                               max_iter    = 10,
                                               nbody       = 4,
                                               H0          = "Hcmf",
                                               thresh_var  = thresh_spf,
                                               thresh_foi  = thresh_spf/50,
                                               thresh_pt   = thresh_spf/2,
                                               ci_conv     = 5e-5,
                                               do_pt       = false,
                                               tol_tucker  = 1e-5,
                                               resolve_ss  = true,
                                               verbose     = 1)
@save "data_spt.jld2" e_var7 v_var cluster_bases
t1=time()
alloc= @allocated σ2 = FermiCG.compute_spt_sigma_norm_blockwise(v_var, cluster_ops, clustered_ham;
                                       H0="Hcmf", nbody=4, thresh_foi=1e-8,
                                       max_number=nothing, opt_ref=true,
                                       ci_tol=1e-6, verbose=1)

var_like7 = σ2 .- e_var7.^2
println("Variance-like, Energy: $(var_like7[1]), $(e_var7[1])")
@printf("Allocated memory: %.3f GB\n", alloc/1e9)
t2 = time()
println("Time taken: $(t2-t1) seconds")
t1=time()
alloc= @allocated σ2 = FermiCG.compute_spt_sigma_norm_blockwise_new(v_var, cluster_ops, clustered_ham;
                                       H0="Hcmf", nbody=4, thresh_foi=1e-8,
                                       max_number=nothing, opt_ref=true,
                                       ci_tol=1e-6, verbose=1)

var_like7 = σ2 .- e_var7.^2
println("Variance-like, Energy: $(var_like7[1]), $(e_var7[1])")
@printf("Allocated memory: %.3f GB\n", alloc/1e9)
t2 = time()
println("Time taken: $(t2-t1) seconds")
println("Variance, Energy")
@printf(" %1.5e, %1.5e\n", var_like2[1], e_var2[1])
@printf("%1.5e, %1.5e\n", var_like3[1], e_var3[1])
@printf(" %1.5e, %1.5e\n", var_like4[1], e_var4[1])
@printf(" %1.5e, %1.5e\n", var_like5[1], e_var5[1])
@printf(" %1.5e, %1.5e\n", var_like6[1], e_var6[1])
@printf(" %1.5e, %1.5e\n", var_like7[1], e_var7[1])