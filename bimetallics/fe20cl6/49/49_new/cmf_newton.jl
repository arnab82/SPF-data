using QCBase
using ClusterMeanField
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using Printf
using ActiveSpaceSolvers

C = npzread("mo_coeffs.npy")
h0 = npzread("ints_h0.npy")
h1 = npzread("ints_h1.npy")
h2 = npzread("ints_h2.npy")
ints = InCoreInts(h0, h1, h2)

Pa = npzread("Pa.npy")
Pb = npzread("Pb.npy")
@printf(" Input energy:    %12.8f\n", compute_energy(ints, RDM1(Pa, Pb)))


init_fspace =  [(5, 0), (4, 4), (5, 0), (12, 12), (0, 0)]
clusters    =  clusters    =  [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [11, 12, 13, 14, 15, 16, 17], [18, 19, 20, 21, 22, 23, 24, 25, 26, 27], [28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39], [40, 41, 42, 43, 44, 45, 46, 47, 48, 49]]

clusters = [MOCluster(i, collect(clusters[i])) for i = 1:length(clusters)]
display(clusters)

d1 = RDM1(n_orb(ints))

ansatze=[FCIAnsatz(length(ci),init_fspace[ci.idx][1],init_fspace[ci.idx][2]) for ci in clusters]

# Do CMF
e_cmf, U, d1 = ClusterMeanField.cmf_oo_newton(ints, clusters, init_fspace, d1, maxiter_oo = 400,
                                              tol_oo=1e-6, 
                                              tol_d1=1e-9, 
                                              tol_ci=1e-11,
                                              verbose=4, 
                                              sequential=true)

ints = orbital_rotation(ints, U)
C = C*U

npzwrite("Ccmf_49_new.npy", C)

@save "data_cmf_49_new.jld2" clusters init_fspace ints d1 e_cmf U 