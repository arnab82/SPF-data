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



σ = FermiCG.build_compressed_1st_order_state(v_var, cluster_ops, clustered_ham, nbody=4, thresh=1e-9, compress_iter=false)
σ = FermiCG.compress(σ, thresh=1e-9)

# Add PT1 space
println()
v_var_=deepcopy(v_var)
time = @elapsed begin
    dim1 = length(v_var_)
    @timeit to "add" FermiCG.nonorth_add!(v_var_, σ)
    @timeit to "compress" v_var_ = FermiCG.compress(v_var_, thresh=1e-9)
    dim2 = length(v_var_)
    @printf(" %-50s", "Variational space increased from: ")
    @printf("%10i → %-10i (%-11s = %8.1e)\n", dim1, dim2, "thresh_pt", 1e-9)
    FermiCG.orthonormalize!(v_var_)
end
H2 = FermiCG.orth_dot(v_var_, v_var_)
var1 = zeros(length(e_var1))
for r in 1:length(e_var1)
    var1[r] = H2[1] - e_var1[r] * e_var1[r]
end
@printf("Variance, Energy: %1.5e, %1.8e\n", var1[1], e_var1[1])

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



σ = FermiCG.build_compressed_1st_order_state(v_var, cluster_ops, clustered_ham, nbody=4, thresh=1e-9, compress_iter=false)
σ = FermiCG.compress(σ, thresh=1e-9)
# Add PT1 space
println()
v_var_=deepcopy(v_var)
time = @elapsed begin
    dim1 = length(v_var_)
    @timeit to "add" FermiCG.nonorth_add!(v_var_, σ)
    @timeit to "compress" v_var_ = FermiCG.compress(v_var_, thresh=1e-9)
    dim2 = length(v_var_)
    @printf(" %-50s", "Variational space increased from: ")
    @printf("%10i → %-10i (%-11s = %8.1e)\n", dim1, dim2, "thresh_pt", 1e-9)
    FermiCG.orthonormalize!(v_var_)
end
H2 = FermiCG.orth_dot(v_var_, v_var_)
var2 = zeros(length(e_var2))
for r in 1:length(e_var2)
    var2[r] = H2[1] - e_var2[r] * e_var2[r]
end
@printf("Variance, Energy: %1.5e, %1.8e\n", var2[1], e_var2[1])



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



σ = FermiCG.build_compressed_1st_order_state(v_var, cluster_ops, clustered_ham, nbody=4, thresh=1e-9,compress_iter=false)
σ = FermiCG.compress(σ, thresh=1e-9)
# Add PT1 space
println()
v_var_=deepcopy(v_var)
time = @elapsed begin
    dim1 = length(v_var_)
    @timeit to "add" FermiCG.nonorth_add!(v_var_, σ)
    @timeit to "compress" v_var_ = FermiCG.compress(v_var_, thresh=1e-9)
    dim2 = length(v_var_)
    @printf(" %-50s", "Variational space increased from: ")
    @printf("%10i → %-10i (%-11s = %8.1e)\n", dim1, dim2, "thresh_pt", 1e-9)
    FermiCG.orthonormalize!(v_var_)
end
H2 = FermiCG.orth_dot(v_var_, v_var_)
var3 = zeros(length(e_var3))
for r in 1:length(e_var3)
    var3[r] = H2[1] - e_var3[r] * e_var3[r]
end
@printf("Variance, Energy: %1.5e, %1.8e\n", var3[1], e_var3[1])



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



σ = FermiCG.build_compressed_1st_order_state(v_var, cluster_ops, clustered_ham, nbody=4, thresh=1e-9,compress_iter=false)
σ = FermiCG.compress(σ, thresh=1e-9)
# Add PT1 space
println()
v_var_=deepcopy(v_var)
time = @elapsed begin
    dim1 = length(v_var_)
    @timeit to "add" FermiCG.nonorth_add!(v_var_, σ)
    @timeit to "compress" v_var_ = FermiCG.compress(v_var_, thresh=1e-9)
    dim2 = length(v_var_)
    @printf(" %-50s", "Variational space increased from: ")
    @printf("%10i → %-10i (%-11s = %8.1e)\n", dim1, dim2, "thresh_pt", 1e-9)
    FermiCG.orthonormalize!(v_var_)
end
H2 = FermiCG.orth_dot(v_var_, v_var_)
var4 = zeros(length(e_var4))
for r in 1:length(e_var4)
    var4[r] = H2[1] - e_var4[r] * e_var4[r]
end
@printf("Variance, Energy: %1.5e, %1.8e\n", var4[1], e_var4[1])

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



σ = FermiCG.build_compressed_1st_order_state(v_var, cluster_ops, clustered_ham, nbody=4, thresh=1e-9,compress_iter=false)
σ = FermiCG.compress(σ, thresh=1e-9)
# Add PT1 space
println()
v_var_=deepcopy(v_var)
time = @elapsed begin
    dim1 = length(v_var_)
    @timeit to "add" FermiCG.nonorth_add!(v_var_, σ)
    @timeit to "compress" v_var_ = FermiCG.compress(v_var_, thresh=1e-9)
    dim2 = length(v_var_)
    @printf(" %-50s", "Variational space increased from: ")
    @printf("%10i → %-10i (%-11s = %8.1e)\n", dim1, dim2, "thresh_pt", 1e-9)
    FermiCG.orthonormalize!(v_var_)
end
H2 = FermiCG.orth_dot(v_var_, v_var_)
var5 = zeros(length(e_var5))
for r in 1:length(e_var5)
    var5[r] = H2[1] - e_var5[r] * e_var5[r]
end
@printf("Variance, Energy: %1.5e, %1.8e\n", var5[1], e_var5[1])



thresh_spf=0.0015
e_var6, v_var = FermiCG.block_sparse_tucker(v_var, cluster_ops, clustered_ham,
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
@save "data_spt.jld2" e_var6 v_var cluster_bases



σ = FermiCG.build_compressed_1st_order_state(v_var, cluster_ops, clustered_ham, nbody=4, thresh=1e-9,compress_iter=false)
σ = FermiCG.compress(σ, thresh=1e-9)
# Add PT1 space
println()
v_var_=deepcopy(v_var)
time = @elapsed begin
    dim1 = length(v_var_)
    @timeit to "add" FermiCG.nonorth_add!(v_var_, σ)
    @timeit to "compress" v_var_ = FermiCG.compress(v_var_, thresh=1e-9)
    dim2 = length(v_var_)
    @printf(" %-50s", "Variational space increased from: ")
    @printf("%10i → %-10i (%-11s = %8.1e)\n", dim1, dim2, "thresh_pt", 1e-9)
    FermiCG.orthonormalize!(v_var_)
end
H2 = FermiCG.orth_dot(v_var_, v_var_)
var6 = zeros(length(e_var6))
for r in 1:length(e_var6)
    var6[r] = H2[1] - e_var6[r] * e_var6[r]
end
@printf("Variance, Energy: %1.5e, %1.8e\n", var6[1], e_var6[1])
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


σ = FermiCG.build_compressed_1st_order_state(v_var, cluster_ops, clustered_ham, nbody=4, thresh=1e-9,compress_iter=false)
σ = FermiCG.compress(σ, thresh=1e-9)
# Add PT1 space
println()
v_var_=deepcopy(v_var)
time = @elapsed begin
    dim1 = length(v_var_)
    @timeit to "add" FermiCG.nonorth_add!(v_var_, σ)
    @timeit to "compress" v_var_ = FermiCG.compress(v_var_, thresh=1e-9)
    dim2 = length(v_var_)
    @printf(" %-50s", "Variational space increased from: ")
    @printf("%10i → %-10i (%-11s = %8.1e)\n", dim1, dim2, "thresh_pt", 1e-9)
    FermiCG.orthonormalize!(v_var_)
end
H2 = FermiCG.orth_dot(v_var_, v_var_)
var7 = zeros(length(e_var7))
for r in 1:length(e_var7)
    var7[r] = H2[1] - e_var7[r] * e_var7[r]
end
@printf("Variance, Energy: %1.5e, %1.8e\n", var7[1], e_var7[1])


println("Variance, Energy")
@printf(" %1.5e, %1.5e\n", var1[1], e_var1[1])
@printf("%1.5e, %1.5e\n", var2[1], e_var2[1])
@printf(" %1.5e, %1.5e\n", var3[1], e_var3[1])
@printf(" %1.5e, %1.5e\n", var4[1], e_var4[1])
@printf(" %1.5e, %1.5e\n", var5[1], e_var5[1])
@printf(" %1.5e, %1.5e\n", var6[1], e_var6[1])
@printf(" %1.5e, %1.5e\n", var7[1], e_var7[1])