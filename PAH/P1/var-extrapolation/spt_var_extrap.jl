using FermiCG
using LinearAlgebra
using Printf
using JLD2
@load "data_benz_p1_cmf.jld2"
@load "data_spt.jld2"
max_roots = 60
# Build Cluster basis
#cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1_n, [3,3,3,3], FermiCG.FockConfig(init_fspace), max_roots=max_roots, verbose=1);

# Build ClusteredOperator
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);

# Build Cluster Operators
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);
nroots = 1

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1_n.a, d1_n.b)


v_var_=deepcopy(v_var)
e_var1 = FermiCG.compute_expectation_value(v_var, cluster_ops, clustered_ham)
σ = FermiCG.build_compressed_1st_order_state(v_var, cluster_ops, clustered_ham, nbody=4, thresh=1e-8)
σ = FermiCG.compress(σ, thresh=1e-9)
# Add PT1 space
println()
time = @elapsed begin
    #zero!(pt1_vec)
    dim1 = length(v_var)
    @timeit to "add" nonorth_add!(v_var, σ)
    @timeit to "compress" v_var = compress(v_var, thresh=1e-9)
    dim2 = length(v_var)
    @printf(" %-50s", "Variational space increased from: ")
    @printf("%10i → %-10i (%-11s = %8.1e)\n", dim1, dim2, "thresh_pt", 1e-9)
    orthonormalize!(v_var)
end
H2 = FermiCG.orth_dot(v_var, v_var)
var1 = zeros(length(e_var1))
for r in 1:length(e_var1)
    var1[r] = H2[1] - e_var1[r] * e_var1[r]
end
@printf("Variance, Energy: %1.5e, %1.5e\n", var1[1], e_var1[1])
v5=FermiCG.compress(v_var_, thresh=6e-4)
e_var2 = FermiCG.compute_expectation_value(v5, cluster_ops, clustered_ham)
σ = FermiCG.build_compressed_1st_order_state(v5, cluster_ops, clustered_ham, nbody=4, thresh=1e-8)
σ = FermiCG.compress(σ, thresh=1e-9)
@time = @elapsed begin
    #zero!(pt1_vec)
    dim1 = length(v5)
    @timeit to "add" nonorth_add!(v5, σ)
    @timeit to "compress" v5 = compress(v5, thresh=1e-9)
    dim2 = length(v5)
    @printf(" %-50s", "Variational space increased from: ")
    @printf("%10i → %-10i (%-11s = %8.1e)\n", dim1, dim2, "thresh_pt", 1e-9)
    orthonormalize!(v5)
end
H2 = FermiCG.orth_dot(v5,v5)
var2 = zeros(length(e_var2))
for r in 1:length(e_var2)
    var2[r] = H2[1] - e_var2[r] * e_var2[r]
end
@printf("Variance, Energy: %1.5e, %1.5e\n", var2[1], e_var2[1])

v4=FermiCG.compress(v_var_, thresh=1e-3)
e_var3 = FermiCG.compute_expectation_value(v4, cluster_ops, clustered_ham)
σ = FermiCG.build_compressed_1st_order_state(v4, cluster_ops, clustered_ham, nbody=4, thresh=1e-8)
σ = FermiCG.compress(σ, thresh=1e-9)
@time = @elapsed begin
    #zero!(pt1_vec)
    dim1 = length(v4)
    @timeit to "add" nonorth_add!(v4, σ)
    @timeit to "compress" v4 = compress(v4, thresh=1e-9)
    dim2 = length(v4)
    @printf(" %-50s", "Variational space increased from: ")
    @printf("%10i → %-10i (%-11s = %8.1e)\n", dim1, dim2, "thresh_pt", 1e-9)
    orthonormalize!(v4)
end
H2 = FermiCG.orth_dot(v4,v4)
var3 = zeros(length(e_var3))
for r in 1:length(e_var3)
    var3[r] = H2[1] - e_var3[r] * e_var3[r]
end
@printf("Variance, Energy: %1.5e, %1.5e\n", var3[1], e_var3[1])
v3=FermiCG.compress(v_var_, thresh=2e-3)
e_var4 = FermiCG.compute_expectation_value(v3, cluster_ops, clustered_ham)
σ = FermiCG.build_compressed_1st_order_state(v3, cluster_ops, clustered_ham, nbody=4, thresh=1e-8)
σ = FermiCG.compress(σ, thresh=1e-9)
@time = @elapsed begin
    #zero!(pt1_vec)
    dim1 = length(v3)
    @timeit to "add" nonorth_add!(v3, σ)
    @timeit to "compress" v3 = compress(v3, thresh=1e-9)
    dim2 = length(v3)
    @printf(" %-50s", "Variational space increased from: ")
    @printf("%10i → %-10i (%-11s = %8.1e)\n", dim1, dim2, "thresh_pt", 1e-9)
    orthonormalize!(v3)
end
H2 = FermiCG.orth_dot(v3,v3)
var4 = zeros(length(e_var4))
for r in 1:length(e_var4)
    var4[r] = H2[1] - e_var4[r] * e_var4[r]
end
@printf("Variance, Energy: %1.5e, %1.5e\n", var4[1], e_var4[1])
println("Variance, Energy")
@printf(" %1.5e, %1.5e\n", var1[1], e_var1[1])
@printf("%1.5e, %1.5e\n", var2[1], e_var2[1])
@printf(" %1.5e, %1.5e\n", var3[1], e_var3[1])
@printf(" %1.5e, %1.5e\n", var4[1], e_var4[1])