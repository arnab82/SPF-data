using QCBase
using Printf
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2


function run_spt_scan()
    for geom in reverse(29:31)
        f_out = @sprintf("%03i.out",geom)
        f_err = @sprintf("%03i.err",geom)

        redirect_stdio(stdout=f_out, stderr=f_err) do
            file1 = "data_cmf_$(geom).jld2"
            println(" CMF File: ", file1)
            @eval @load $file1
            println(" Enuc: ", cmf_ints.h0)

            M = 100

            init_fspace = [(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
            ref_fspace = FockConfig(init_fspace)
            ecore = cmf_ints.h0

            cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(cmf_ints, clusters, d1, [1,1,1,1,1,1,1], ref_fspace, max_roots=M, verbose=4);

            clustered_ham = FermiCG.extract_ClusteredTerms(cmf_ints, clusters)
            cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, cmf_ints);

            FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, cmf_ints, d1.a, d1.b);

            nroots = 15 
            ci_vector=FermiCG.BSTstate(clusters, ref_fspace, cluster_bases, R=nroots)
            # ci_vector = FermiCG.add_spin_focksectors(ci_vector)
            FermiCG.add_single_excitons!(ci_vector, ref_fspace, nroots)
            display(ci_vector)
            FermiCG.eye!(ci_vector)   
            eguess, vguess = FermiCG.ci_solve(ci_vector, cluster_ops, clustered_ham)
            @time e2guess = FermiCG.compute_pt2_energy(vguess, cluster_ops, clustered_ham, thresh_foi=1e-8);
            e_var, v_var = block_sparse_tucker(vguess, cluster_ops, clustered_ham,
                                   max_iter    = 20,
                                   nbody       = 4,
                                   H0          = "Hcmf",
                                   thresh_var  = 1e-3,
                                   thresh_foi  = 1e-6,
                                   thresh_pt   = 2e-5,
                                   ci_conv     = 1e-5,
                                   ci_max_iter = 100,
                                   do_pt       = true,
                                   resolve_ss  = false,
                                   tol_tucker  = 1e-4,
                                   solver      = "davidson")
            @time e2 = FermiCG.compute_pt2_energy(v_var, cluster_ops, clustered_ham, thresh_foi=1e-6,prescreen   = true,compress_twice = true)
    
            e_var .+= cmf_ints.h0
            eguess .+= cmf_ints.h0
            e2guess .+= cmf_ints.h0
            e2 .+= cmf_ints.h0

            print("*EGuess: ")
            [@printf("%12.8f ", i) for i in eguess]
            println()

            print("*E2Guess: ")
            [@printf("%12.8f ", i) for i in e2guess]
            println()

            print("*E0: ")
            [@printf("%12.8f ", i) for i in e_var]
            println()

            print("*E2: ")
            [@printf("%12.8f ", i) for i in e2]
            println()

            @save @sprintf("data_spt_%03i.jld2",geom) vguess v_var  eguess e_var e2 e2guess 
        end
    end
end

run_spt_scan()


