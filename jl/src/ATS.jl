module ATSMIX
    using Distributions
    using Formatting
    using Memento
    using Random
    using RCall

    using Memento
    using Compat: @__MODULE__  # requires a minimum of Compat 0.26. Not required on Julia 0.7

    # Create our module level LOGGER (this will get precompiled)
    const LOGGER = Memento.config!("debug"; fmt="[{date} - {level} | {name}]: {msg}")

    # const LOGGER = getlogger(@__MODULE__)

    # Register the module level LOGGER at runtime so that folks can access the LOGGER via `get_LOGGER(MyModule)`
    # NOTE: If this line is not included then the precompiled `MyModule.LOGGER` won't be registered at runtime.
    # function __init__()
    #     Memento.register(LOGGER)
    # end

    Number = Union{Int64, Float64, Int}

    export emptyData, fit, fit_by_R

    struct Param
        st_arr::Vector
        en_arr::Vector
        mu_f::Number
        sigma_f::Number
        unif_log_lik::Number
        max_beta::Number
    end

    function createParam(
        st_arr::Vector, en_arr::Vector, 
        mu_f::Number, sigma_f::Number, 
        L::Number, max_beta::Number;
        single_end_mode::Bool = false
    )::Param
        unif_log_lik = lik_f0(L; do_log = true)
        if single_end_mode
            st_arr = en_arr .+ mu_f
            unif_log_lik = lik_f0_single(mu_f, sigma_f, L; do_log = true)
        end
        return Param(st_arr, en_arr, mu_f, sigma_f, unif_log_lik, max_beta)
    end

    mutable struct EMData
        ws::Union{AbstractArray, Nothing}
        alpha_arr::Union{Vector, Nothing}
        beta_arr::Union{Vector, Nothing}
        lb_arr::Union{Vector, Nothing}
        bic::Union{Number, Nothing}
        label
    end

    Base.show(io::IO, self::EMData) = begin
        res = Vector{String}()
        for i in [self.ws, self.alpha_arr, self.beta_arr, self.lb_arr, self.label]
            if isnothing(i)
                i = "NA"
            else
                i = join(map(string, i), ",")
            end
            push!(res, i)
        end

        if isnothing(self.bic)
            push!(res, "NA")
        else
            push!(res, string(self.bic))
        end

        print(io, join(res, "\t"))
    end

    function emptyData()::EMData
        return EMData(nothing, nothing, nothing, nothing, nothing, nothing)
    end


    function dnorm(weights::Union{Number, Vector}; mean::Number, sd::Number, do_log::Bool = false)::Union{Number, Vector}
        if do_log
            return logpdf(Normal(mean, sd), weights)
        end
        return pdf(Normal(mean, sd), weights)
    end

    # p(r|l)
    function lik_r_l(l_arr::Vector, r_arr::Vector; mu_f::Number, sigma_f::Number, do_log::Bool=false)::Vector
        return dnorm(l_arr .- r_arr; mean = mu_f, sd = sigma_f, do_log = do_log)
    end

    # p(l|alpha, beta)
    function lik_l_ab(l_arr::Vector, alpha::Number, beta::Number; do_log::Bool=false)::Vector
        return dnorm(l_arr; mean = alpha, sd = beta, do_log = do_log)
    end

    # p(l,r|alpha,beta)
    function lik_lr_ab(l_arr::Vector, r_arr::Vector, alpha::Number, beta::Number; params::Param, do_log::Bool=false)
        if do_log
        return lik_l_ab(l_arr, alpha, beta, do_log=true) .+ lik_r_l(l_arr, r_arr; mu_f=params.mu_f, sigma_f=params.sigma_f, do_log=true)
        else
        return ws .* lik_l_ab(l_arr,alpha, beta) .* lik_r_l(l_arr, r_arr; mu_f=params.mu_f, sigma_f=params.sigma_f)
        end
    end

    function cal_z_k(
        alpha_arr::Vector, beta_arr::Vector, ws::Vector, k::Int,
        log_zmat; params::Param)::AbstractArray

        K = length(ws) - 1

        if 0 < k <= K
            log_zmat[:, k] = log(ws[k]) .+ lik_lr_ab(
                params.st_arr, params.en_arr, alpha_arr[k], beta_arr[k]; 
                params=params, do_log = true
            )
        else
            log_zmat[:, K+1] .= log(ws[K + 1]) + params.unif_log_lik
        end
        return log_zmat
    end

    # Z is log likelihood
    function norm_z(Z::AbstractArray)::AbstractArray
        tmp_max_vec = maximum(Z, dims=2)
        Z = exp.(Z .- tmp_max_vec)
        return Z ./ sum(Z, dims=2)
    end

    # maximize ws given Z
    function maximize_ws(Z)
        res = sum(Z, dims = 1) / size(Z)[1]
        if length(size(res)) > 1
            res = res[1, :]
        end
        return res
    end

    function mstep(alpha_arr::Vector, beta_arr::Vector, ws::AbstractArray, Z::AbstractArray, k::Int; params::Param)::EMData
        K = length(ws) - 1 # last component is uniform component
        # println(string("size of Z in mstep: ", size(Z)))
        new_ws = maximize_ws(Z)

        new_alpha_arr = alpha_arr
        new_alpha_arr[k] = sum(Z[:,k] .* params.st_arr) / sum(Z[:,k])
        
        new_beta_arr = beta_arr;
        new_beta_arr[k] = sqrt.( sum(Z[:,k] .* (params.st_arr .- new_alpha_arr[k]) .^ 2 ) / sum(Z[:,k]) )

        if new_beta_arr[k] > params.max_beta
            new_beta_arr[k] = params.max_beta
        end

        return EMData(new_ws, new_alpha_arr, new_beta_arr, Vector(), NaN, nothing)
    end


    function elbo(log_zmat::AbstractArray, Z::AbstractArray)::Number
        LZ = deepcopy(Z)
        LZ[findall(!isequal(0), Z)] .= log.(Z[findall(!isequal(0), Z)])
        entropy = -1 .* Z .* LZ

        lb = exp_log_lik(log_zmat, Z) .+ sum(entropy)
        if isnan(lb)
            throw(string("lower bounder is na: ", lb))
        end
        return lb
    end

    function exp_log_lik(log_zmat::AbstractArray, Z::AbstractArray)::Number
        ZZ = Z .* log_zmat
        ZZ[findall(isequal(0), Z)] .= 0
        return sum(ZZ)
    end

    function lik_f0(L::Number; do_log::Bool=true)::Number
        if do_log
            return -2 * log(L)
        end
        return 1 / L / L
    end

    function lik_f0_single(mu_f::Number, sigma_f::Number, L::Number; do_log::Bool = true)::Number
        if do_log
            return -log(L) .+ dnorm(mu_f; mean=mu_f, sd=sigma_f, do_log=do_log)
        else
            return 1/L .* dnorm(mu_f; mean=mu_f, sd=sigma_f, do_log=do_log)
        end
    end

    function lik_r_s(r_arr::Vector, s::Number, do_log::Bool=false)::Vector
        res = [r_arr .<= s] ./ s
        if do_log
            return log.(res)
        else
            return res
        end
    end

    # generate random k such that each K is a group and no consecutive elements are the same
    function gen_k_arr(K::Int64, n::Int64)::Vector{Int64}

        if K <= 1
            return [K for _ in 1:n]
        else
            res = Vector()
            count_index = 0
            pool = [x for x in 1:K]
            last = nothing
            while count_index < n
                count_index += 1
                pool = shuffle(pool)

                if pool[1] == last
                    swap_with = rand(1:length(pool))
                    pool[1], pool[swap_with] = pool[swap_with], pool[1]
                end

                append!(res, pool)
                last = pool[length(pool)]
            end

            return res
        end
    end


    function cal_bic(log_zmat::AbstractArray, Z::AbstractArray)::Number
        N, K = size(Z)
        K = K - 1
        return -2 * exp_log_lik(log_zmat, Z) + (3 * K + 1) * log(N)
    end

    # replace query with closest values in ref_arr
    function replace_with_closest(ref_arr::Vector, query_arr::Vector)::Vector
        n = length(query_arr)
        res = zeros(n)

        for i = 1:n
            tmpind = findmin(abs.(ref_arr .- query_arr[i]))
            res[i] = ref_arr[tmpind]
        end
        return res 
    end

    # perform inference given alpha_arr and beta_arr
    function fixed_inference(alpha_arr::Vector, beta_arr::Vector, n_frag::Number, nround::Int; params::Param, verbose::Bool=false)::EMData

        lb, lb_arr = -Inf, Array{Number}(undef, nround, 1)
        N = n_frag
        K = length(alpha_arr)

        ws = [rand(Uniform()) + 1 for _ in 1:K+1]
        if K == 0
            ws[1] = ws[1] + 1
        else
            ws[1:K] = ws[1:K] .+ 1
        end
        ws = ws ./ sum(ws)

        k_arr = gen_k_arr(length(alpha_arr), nround)
        log_zmat = zeros( N, K + 1)

        for k = 1:K
            log_zmat = cal_z_k(alpha_arr, beta_arr, ws, k,log_zmat; params=params)
        end

        Z = nothing
        i = 1
        for i = 1:nround
            if verbose
                debug(LOGGER, Formatting.format(
                    FormatExpr("fixed_inference: iteration={}, lb={}"), i, lb
                ))
            end

            log_zmat = cal_z_k(alpha_arr, beta_arr, ws,k_arr[i], log_zmat, params=params)
            
            Z = norm_z(log_zmat)
            ws = maximize_ws(Z)

            lb_new = elbo(log_zmat, Z)
            lb_arr[i] = lb_new
            if lb_new == -Inf
                lb = -Inf
                break
            end
            if abs(lb_new - lb)  < abs(1e6 * lb)
                break
            else
                lb = lb_new
            end
        end
        if verbose
            debug(LOGGER, Formatting.format(
                FormatExpr("fixed_inference: Run all {} iterations. lb={}"), 
                i, lb
            ))
        end
        bic = NaN
        if !isnothing(Z)
            bic = cal_bic(log_zmat, Z)
        end
        
        # label = argmax(log_zmat, dims=2)
        if verbose
            debug(LOGGER, Formatting.format(
                FormatExpr("fixed_inference: bic={}; estimated ws: {}; estimated alpha:  {}; estimated beta: {}"), 
                bic, ws, alpha_arr, beta_arr
            ))
        end

        lb_arr = lb_arr[findall(x -> !isnan(x), lb_arr)]
        return EMData(ws, alpha_arr, beta_arr, lb_arr, bic, nothing)
    end

    # perform inference for K components
    # theta_arr => theta_win_mat
    function em_algo(
        ws::Vector, para_mat::EMData, nround::Int;
        params::Param, verbose::Bool=false)::EMData
        lb = -Inf
        lb_arr = [NaN for _ = 1:nround]
        K = length(ws) - 1

        # assume this is a NamedTuple
        alpha_arr = para_mat.alpha_arr
        beta_arr = para_mat.beta_arr

        k_arr = gen_k_arr(length(alpha_arr), nround)
        log_zmat = zeros(length(params.st_arr), K+1)

        for k = 1:(K+1)
            log_zmat = cal_z_k(alpha_arr, beta_arr, ws, k, log_zmat; params=params)
        end

        Z = nothing
        i = 1
        for i = 1:nround
            # println(string(i, ws))
            log_zmat = cal_z_k(alpha_arr, beta_arr, ws, k_arr[i], log_zmat; params=params)
            
            # println(string("size of log_zmat: ", size(log_zmat)))
            Z = norm_z(log_zmat)
            # println(string("size of Z: ", size(Z)))
            res = mstep(alpha_arr, beta_arr, ws, Z, k_arr[i]; params=params) # , log_zmat

            alpha_arr = res.alpha_arr
            beta_arr = res.beta_arr
            ws = res.ws

            lb_new = elbo(log_zmat, Z)
            lb_arr[i] = lb_new
            if lb_new == -Inf
                lb = -Inf
                break
            end
            
            if abs(lb_new - lb) < abs(1e-6 * lb)
                break
            else
                lb = lb_new
            end
        end

        if verbose
            debug(LOGGER, Formatting.format(
                FormatExpr("em_algo: Run all {} interactions. lb={}"), i, lb
            ))
        end

        bic = NaN
        if !isnothing(Z) 
            bic = cal_bic(log_zmat, Z)
        end
        
        if verbose
            debug(LOGGER, Formatting.format(
                FormatExpr("em_algo: bic={}; estimated ws: {}; estimated alpha:  {}; estimated beta: {}"), 
                bic, ws, alpha_arr, beta_arr
            ))
        end

        lb_arr = lb_arr[findall(x -> !isnan(x), lb_arr)]

        sorted_inds = sortperm(alpha_arr)
        alpha_arr = alpha_arr[sorted_inds]
        beta_arr = beta_arr[sorted_inds]
        ws[1:K] = ws[sorted_inds]

        return EMData(ws, alpha_arr, beta_arr, lb_arr, bic, nothing)
    end

    function init_para(n_ats::Int, st_arr::Vector, data::EMData)::EMData
        min_pos = minimum(st_arr)
        max_pos = maximum(st_arr)

        data.alpha_arr =  min_pos .+ (max_pos - min_pos) .* rand(n_ats)
        data.beta_arr = rand(Uniform(10, 70), n_ats)

        return data
    end

    function em_optim0(n_ats::Int; params::Param, n_trial::Int=20, verbose::Bool = false)

        lb_arr = [-Inf for _ = 1:n_trial]
        bic_arr = [-Inf for _ = 1:n_trial]
        res_list = Array{EMData}(undef, n_trial, 1)

        for i = 1:n_trial
            try
                ws = [rand(Uniform()) + 1 for _ in 1:(n_ats+1)]
                ws[n_ats + 1] = ws[n_ats + 1] - 1
                ws = ws / sum(ws)

                # initilize alpha_arr and beta_arr, considered pre_alpha_arr and pre_beta_arr
                para_mat = init_para(n_ats, params.st_arr, emptyData())

                if verbose
                    # println(para_mat)
                    debug(LOGGER, string("em_optim0: ", para_mat))
                end

                res_list[i] = em_algo(ws, para_mat, n_trial; params=params, verbose=verbose)
            catch e
                if verbose
                    println(e)
                    debug(LOGGER, Formatting.format(
                        FormatExpr("em_optim0: Error found  in {} trial. Next - {}"), i, e
                    ))
                end
                res_list[i] = emptyData()
            end

            if isnothing(res_list[i].alpha_arr)
                continue
            end

            lb_arr[i] =  res_list[i].lb_arr[length(res_list[i].lb_arr)]
            bic_arr[i] =  res_list[i].lb_arr[length(res_list[i].bic)]

            if verbose
                debug(LOGGER, Formatting.format(
                    FormatExpr("em_optim0: K = {}, {}_trial, n_trial_{}: ws: {}; alpha: {}; beta: {};lb = {}; bic = {}"), 
                    n_ats, i, n_trial, round.(res_list[i].ws),
                    res_list[i].alpha_arr, res_list[i].beta_arr, 
                    res_list[i].lb_arr[length(res_list[i].lb_arr)],
                    res_list[i].bic 
                ))
            end
        end

        min_ind = argmin(bic_arr)
        return res_list[min_ind]
    end

    function rm_component(res::EMData, min_ws::Number; n_frag::Number, params::Param, nround::Int=200, verbose::Bool = false)

        K = length(res.alpha_arr)
        res.alpha_arr = round.(Int, res.alpha_arr)
        res.beta_arr = round.(Int, res.beta_arr)

        rm_inds = findall(x -> x < min_ws, res.ws[1:K])

        if sum(rm_inds) == 0
            return res
        end

        alpha_arr = Vector()
        beta_arr = Vector()

        for i = 1:length(res.alpha_arr)
            if i in rm_inds
                continue
            end
            push!(alpha_arr, res.alpha_arr[i])
            push!(beta_arr, res.beta_arr[i])
        end

        return fixed_inference(alpha_arr, beta_arr, n_frag, nround, params=params, verbose=verbose)
    end


    function fit_by_R(
        # maximum number of ATS sites
        n_max_ats::Int, 
        # minimum number of ATS sites
        n_min_ats::Int, 

        # information for all DNA fragments, must be n_frag x 1 vector if specified
        # l, start location of each DNA fragment, from 5' to 3'on the 5'-UTR
        st_arr::Vector, 
        # r, end location of each DNA fragment
        en_arr::Vector; 

        # length of UTR region
        L::Int=nothing, 

        # fragment size information
        # fragment length mean
        mu_f::Int = 350, 
        # fragment length standard deviation
        sigma_f::Int = 50, 

        # pa site information
        # minimum weight of ATS site
        min_ws::Float64 = 0.01, 
        # maximum std for ATS site
        max_beta::Float64 = 30.0,
        # inference with fixed parameters
        fixed_inference_flag::Bool = false,

        #single end mode
        single_end_mode::Bool = false,
        verbose::Bool = false,
        error_log=nothing
    )::EMData

        try
            atsmix = R"""
                # ats mixture model inference
                atsmix = function(n_max_ats, # maximum number of ATS sites
                                n_min_ats=1, # minimum number of ATS sites
                                
                                # information for all DNA fragments, must be n_frag x 1 vector if specified
                                st_arr, # l, start location of each DNA fragment, from 5' to 3'on the 5'-UTR
                                en_arr, # r, end location of each DNA fragment

                                L=NULL,  # length of UTR region
                                
                                # fragment size information
                                mu_f=350, # fragment length mean
                                sigma_f=50, # fragment length standard deviation
                                
                                # pa site information
                                min_ws = 0.01, # minimum weight of ATS site
                                max_beta = 50, # maximum std for ATS site
                                
                                # inference with fixed parameters
                                fixed_inference_flag = FALSE,
                                
                                #single end mode
                                single_end_mode = FALSE, 
                                debug=T){
                
                cal_z_k = function(alpha_arr, beta_arr, ws, k, log_zmat){
                    K = length(ws) - 1  # last component is uniform component
                    if(k<=K){
                        log_zmat[,k] = log(ws[k]) + lik_lr_ab(st_arr, en_arr, alpha_arr[k], beta_arr[k], log=T)
                    }else{
                        log_zmat[,K+1] = log(ws[K+1]) + unif_log_lik
                    }
                    return(log_zmat)
                }
                
                # Z is log likelihood
                norm_z = function(Z){
                    tmp_max_vec = apply(Z,1,max)
                    Z = exp( Z - tmp_max_vec )
                    Z = Z/apply(Z,1,sum)
                    return(Z)
                }

                # maximize ws given Z
                maximize_ws = function(Z){
                    ws = apply(Z,2,sum)/dim(Z)[1]
                    return( ws )
                }
                
                mstep = function(alpha_arr, beta_arr, ws, Z, k){
                    
                    K = length(ws) - 1 # last component is uniform component
                    new_ws = maximize_ws(Z)
                    new_alpha_arr = alpha_arr
                    new_alpha_arr[k] = sum(Z[,k]*st_arr)/sum(Z[,k])
                    
                    new_beta_arr = beta_arr;
                    new_beta_arr[k] = sqrt( sum(Z[,k]*(st_arr-new_alpha_arr[k])^2)/sum(Z[,k]) )
                    if(new_beta_arr[k]>max_beta){
                        new_beta_arr[k] = max_beta
                    }
                    
                    return(list(alpha_arr=new_alpha_arr, beta_arr=new_beta_arr, ws=new_ws))
                }
                
                elbo = function(log_zmat, Z){
                    LZ = Z
                    LZ[Z!=0] = log(Z[Z!=0])
                    entropy = -1 * Z * LZ
                    
                    lb = exp_log_lik(log_zmat, Z) + sum(entropy)
                    if(is.na(lb)){
                        stop("lower bounder is na.")
                    }
                    return (lb)
                }
                
                # calculate the expected log joint likelihood
                exp_log_lik = function(log_zmat, Z){
                    ZZ = Z * log_zmat
                    ZZ[Z==0] = 0
                    
                    return(sum(ZZ))
                }
                
                # uniform component likelihood
                lik_f0 = function(log=F){
                    if(log){
                        return(-2*log(L))
                    }else{
                        return(1/L/L)
                    }
                }
                
                # uniform component likelihood for single end case
                lik_f0_single = function(log=F){
                    if(log){
                        return( -log(L)+dnorm(mu_f,mean=mu_f,sd=sigma_f,log=log) )
                    }else{
                        return( 1/L*dnorm(mu_f,mean=mu_f,sd=sigma_f,log=log) )
                    }
                }
                
                # p(l,r|alpha,beta)
                lik_lr_ab = function(l_arr, r_arr, alpha, beta, log=F){
                    if(log){
                        return(lik_l_ab(l_arr,alpha, beta, log=T)+lik_r_l(l_arr, r_arr, log=T))
                    }else{
                        return(ws * lik_l_ab(l_arr,alpha, beta) * lik_r_l(l_arr, r_arr))
                    }
                } 
                
                # p(l|alpha, beta)
                lik_l_ab = function(l_arr, alpha, beta, log=F){
                    return( dnorm(l_arr, mean=alpha, sd=beta, log=log)  )
                }
                
                # p(r|l)
                lik_r_l = function(l_arr, r_arr, log=F){
                    return( dnorm(l_arr-r_arr, mean=mu_f, sd=sigma_f, log=log)  )
                }
                
                # generate random k such that each K is a group and no consecutive elements are the same
                gen_k_arr = function(K, n){
                    if(K==0){
                        return(rep(0,n))
                    }
                    if(K==1){
                        return(rep(1,n))
                    }
                    if(K==2){
                        nn = ceiling(n/K)
                        res = rep(c(1,2),nn)
                        return(res[seq(n)])
                    }
                    nn = ceiling(n/K)
                    res = rep(0,nn*K)
                    res[seq(1,K)]=sample(K)
                    for(i in seq(2,nn)){
                        st = (i-1)*K
                        res[seq(st+1,st+K)]=sample(K)
                        if(res[st]==res[st+1]){
                            tmpind = sample(seq(2,K),1)
                            tmp = res[st+tmpind]
                            res[st+tmpind] = res[st+1]
                            res[st+1] = tmp
                        }
                    }
                    return(res[seq(n)])
                }
                
                cal_bic = function(log_zmat, Z){
                    N = dim(Z)[1]
                    K = dim(Z)[2]-1
                    
                    res = -2*exp_log_lik(log_zmat, Z) + (3*K+1) *log(N)  # the smaller bic, the better model
                    return(res)
                }
                
                # replace query with closest values in ref_arr
                replace_with_closest = function(ref_arr, query_arr){
                    n = length(query_arr)
                    res = rep(0,n)
                    for(i in seq(n)){
                        tmpind = which.min( abs(ref_arr-query_arr[i]) )
                        res[i] = ref_arr[tmpind]
                    }
                    return(res)
                }
                
                # perform inference given alpha_arr and beta_arr
                fixed_inference = function(alpha_arr, beta_arr, debug=F){
                    lb = -Inf
                    lb_arr = rep(NA,nround)
                    N = n_frag
                    K = length(alpha_arr)
                    
                    ws = runif(K+1)
                    ws[1:K] = ws[1:K]+1
                    ws = ws/sum(ws)

                    k_arr = gen_k_arr(length(alpha_arr), nround)
                    
                    log_zmat = matrix( 0, nrow=N, ncol=K+1 )
                    for(k in seq(K+1)){
                        log_zmat = cal_z_k(alpha_arr, beta_arr, ws, k, log_zmat)
                    }
                    
                    for(i in seq(nround)){
                        # estep
                        log_zmat = cal_z_k(alpha_arr, beta_arr, ws, k_arr[i], log_zmat)
                        Z = norm_z(log_zmat)
                    
                        # mstep, alpha and beta are fixed now
                        ws = maximize_ws(Z)  
                        
                        lb_new = elbo(log_zmat, Z)
                        lb_arr[i] = lb_new
                        
                        if(lb_new==-Inf){
                            lb = -Inf
                            break
                        }
                        
                        if( abs(lb_new-lb) < abs(1e-6*lb) )
                        {
                            break
                        }else{
                            lb = lb_new
                        }
                    }
                                        
                    bic = cal_bic(log_zmat, Z)
                    label = apply(log_zmat, 1, which.max)
                    
                    lb_arr = lb_arr[!is.na(lb_arr)]
                    return(list(ws=ws, alpha_arr=alpha_arr, beta_arr=beta_arr, lb_arr=lb_arr, bic=bic, label=label))
                }
                
                # perform inference for K components
                # theta_arr => theta_win_mat
                em_algo = function(ws, para_mat, debug=F){
                    lb = -Inf
                    lb_arr = rep(NA,nround)
                    N = n_frag
                    K = length(ws)-1
                    
                    alpha_arr = para_mat[,1]
                    beta_arr = para_mat[,2]
                    
                    k_arr = gen_k_arr(length(alpha_arr), nround)
                    log_zmat = matrix( 0, nrow=N, ncol=K+1 )
                    for(k in seq(K+1)){
                        log_zmat = cal_z_k(alpha_arr, beta_arr, ws, k, log_zmat)
                    }
                    
                    for(i in seq(nround)){
                        # estep
                        log_zmat = cal_z_k(alpha_arr, beta_arr, ws, k_arr[i], log_zmat)
                        temp = log_zmat
                        Z = norm_z(log_zmat)

                        res = mstep(alpha_arr, beta_arr, ws, Z, k_arr[i])

                        alpha_arr = res$alpha_arr
                        beta_arr = res$beta_arr
                        ws = res$ws
                        
                        lb_new = elbo(log_zmat, Z)
                        lb_arr[i] = lb_new
                        
                        if(lb_new==-Inf){
                            lb = -Inf
                            break
                        }
                        
                        if( abs(lb_new-lb) < abs(1e-6*lb) )
                        {
                            break
                        }else{
                            lb = lb_new
                        }
                    }
                    
                    bic = cal_bic(log_zmat, Z)              
                    lb_arr = lb_arr[!is.na(lb_arr)]

                    sorted_inds = order(alpha_arr)
                    alpha_arr = alpha_arr[sorted_inds]
                    beta_arr = beta_arr[sorted_inds]
                    ws[1:K] = ws[sorted_inds]
                    return(list(ws=ws, alpha_arr=alpha_arr, beta_arr=beta_arr, lb_arr=lb_arr, bic=bic))
                }
                
                init_para = function(n_ats, debug=F){
                    min_pos = min(st_arr)
                    max_pos = max(st_arr)
                    alpha_arr = min_pos + (max_pos-min_pos)*runif(n_ats)
                    beta_arr = runif(n_ats,10,70)
                    return(cbind(alpha_arr, beta_arr))
                }
                
                em_optim0 = function(n_ats, debug=F){
                    n_trial=20
                    
                    lb_arr = rep(-Inf,n_trial)
                    bic_arr = rep(Inf,n_trial)
                    res_list = vector("list",n_trial)
                    
                    for(i in seq(n_trial)){                    
                        res_list[[i]] = tryCatch({
                            ws = runif(n_ats+1)+1  # last component is for uniform component
                            ws[n_ats+1] = ws[n_ats+1]-1
                            ws = ws/sum(ws)
                            
                            # initilize alpha_arr and beta_arr, considered pre_alpha_arr and pre_beta_arr
                            para_mat = init_para(n_ats)
                            
                            em_algo(ws, para_mat, debug=debug)
                        },error=function(e){
                            print(e)
                            list(ws=NULL,alpha_arr=NULL,beta_arr=NULL,lb_arr=-Inf,bic=Inf)
                        })
                        if(is.null(res_list[[i]])){
                            next()
                        }
                        lb_arr[i] = tail(res_list[[i]]$lb_arr,1)
                        bic_arr[i] = res_list[[i]]$bic
                    }
                    
                    min_ind = which.min(bic_arr)
                    res = res_list[[min_ind]]
                                        
                    return(res)
                }
                
                # remove components with weight less than min_ws
                rm_component = function(res, min_ws){
                    K = length(res$alpha_arr)
                    res$alpha_arr = round(res$alpha_arr)
                    res$beta_arr = round(res$beta_arr)
                    rm_inds = which(res$ws[1:K]<min_ws)
                    
                    if(length(rm_inds)==0){
                        return(res)
                    }else{                    
                        alpha_arr = res$alpha_arr[-rm_inds]
                        beta_arr = res$beta_arr[-rm_inds]
                        res = fixed_inference(alpha_arr,beta_arr)
                        return(res)
                    }
                }
                
                n_frag = length(st_arr)
                stopifnot(n_frag==length(en_arr))
                nround=50

                if (typeof(st_arr) == "list") {
                    st_arr = unlist(st_arr)
                }

                if (typeof(en_arr) == "list") {
                    en_arr = unlist(en_arr)
                }

                lb_arr = rep(-Inf,n_max_ats)
                bic_arr = rep(Inf,n_max_ats)
                res_list = vector("list",n_max_ats)
                
                if(single_end_mode){
                    st_arr = en_arr + mu_f
                    unif_log_lik = lik_f0_single(log = T)
                }else{
                    unif_log_lik = lik_f0(log = T)
                }
                

                for(i in seq(n_max_ats,n_min_ats)){
                    res_list[[i]] = em_optim0(i, debug=debug)
                    lb_arr[i] = tail(res_list[[i]]$lb_arr,1)
                    bic_arr[i] = res_list[[i]]$bic
                }

                min_ind = which.min(bic_arr)
                res = res_list[[min_ind]]
                
                if(is.null(res$ws)){
                    return(res)
                }
                                
                # remove low weight component
                res = rm_component(res, min_ws)
                
                # calculate read assignment
                if(!'label' %in% names(res)){
                    N = n_frag
                    K = length(res$ws)-1
                    log_zmat = matrix( 0, nrow=N, ncol=K+1 )
                    for(k in seq(K+1)){
                        log_zmat = cal_z_k(res$alpha_arr, res$beta_arr, res$ws, k, log_zmat)
                    }
                    # Z = norm_z(log_zmat)
                    label = apply(log_zmat, 1, which.max)
                    res$label = label
                }
                
                if(single_end_mode){
                    res$beta_arr = res$beta_arr + sigma_f 
                }
                
                return(res)
            }


            atsmix(n_max_ats=$n_max_ats, 
                    n_min_ats=$n_min_ats,
                    st_arr=$st_arr, 
                    en_arr=$en_arr, 
                    L=$L, 
                    mu_f=$mu_f,
                    sigma_f=$sigma_f,
                    min_ws = $min_ws,
                    max_beta = $max_beta,
                    fixed_inference_flag = FALSE,
                    single_end_mode = $single_end_mode
            )
            """

            res = rcopy(atsmix)
            res = collect(res)
            # println(res)
            ws = res[1][2]
            if isa(ws, Number)
                ws = [ws]
            end

            alpha_arr = res[2][2]
            if isa(alpha_arr, Number)
                alpha_arr = [alpha_arr]
            end

            beta_arr = res[3][2]
            if isa(beta_arr, Number)
                beta_arr = [beta_arr]
            end

            lb_arr = res[4][2]
            if isa(lb_arr, Number)
                lb_arr = [lb_arr]
            end

            label = nothing
            if length(res) >= 6
                label = res[6][2]
            end

            if isnothing(alpha_arr)
                return emptyData()
            end

            return EMData(
                ws, alpha_arr, 
                beta_arr, lb_arr,
                res[5][2], label
            )
        catch e
            
            if isnothing(error_log)
                error_log = "error.log"
            end

            open(error_log, "a+") do w
                write(w, string(
                    string(e), "\t",
                    join(map(string, st_arr), ","), "\t",
                    join(map(string, en_arr), ","), "\t",
                    L, "\n"
                ))
                close(w)
            end

            warn(LOGGER, e)
            return emptyData()
        end
    end

    function fit(
        # maximum number of ATS sites
        n_max_ats::Int, 
        # minimum number of ATS sites
        n_min_ats::Int, 

        # information for all DNA fragments, must be n_frag x 1 vector if specified
        # l, start location of each DNA fragment, from 5' to 3'on the 5'-UTR
        st_arr::Vector, 
        # r, end location of each DNA fragment
        en_arr::Vector; 

        # length of UTR region
        L::Int=nothing, 

        # fragment size information
        # fragment length mean
        mu_f::Int = 350, 
        # fragment length standard deviation
        sigma_f::Int = 30, 

        # pa site information
        # minimum weight of ATS site
        min_ws::Float64 = 0.01, 
        # maximum std for ATS site
        max_beta::Float64 = 30.0,

        # inference with fixed parameters
        fixed_inference_flag::Bool = false,

        #single end mode
        single_end_mode::Bool = false,
        verbose::Bool = false,

        nround::Int = 50,
        using_R::Bool = false,
        error_log = nothing
    )

        try
            n_frag = length(st_arr)

            if n_frag != length(en_arr)
                throw("the length of st_arr != en_arr")
            end
    
            nround=50
    
            lb_arr = [-Inf for _ = 1:n_max_ats]
            bic_arr = [Inf for _ = 1:n_max_ats]
            res = nothing
            
            params = createParam(
                st_arr, en_arr, mu_f, sigma_f, L, max_beta,
                single_end_mode=single_end_mode
            )
            
            for i = n_max_ats:-1:n_min_ats
                temp = em_optim0(i; params=params, verbose=verbose)
    
                if isnothing(res)
                    res = temp
                elseif !isnothing(temp.bic) && !isnan(temp.bic) && !isinf(temp.bic)
                    if temp.bic < res.bic
                        res = temp
                    end
                end
            end
            
            if isnothing(res.ws) 
                debug(LOGGER, "fit: Inference failed. No results available.")
                return res
            end
                
            # remove low weight component
            res = rm_component(res, min_ws; n_frag=n_frag, verbose=verbose, params=params)
            
            # calculate read assignment
            if isnothing(res.label)
                N = n_frag
                K = length(res.ws) - 1
                log_zmat = zeros(N, K + 1)
                for k = 1:(K + 1)
                    log_zmat = cal_z_k(res.alpha_arr, res.beta_arr, res.ws, k, log_zmat; params=params)
                end
                # Z = norm_z(log_zmat)
                label = argmax.(log_zmat)
                res.label = label
            end
            
            if single_end_mode
                res.beta_arr = res.beta_arr .+ sigma_f 
            end 
            return res
        catch e
            if isnothing(error_log)
                error_log = "error.log"
            end

            open(error_log, "a+") do w
                write(w, string(
                    string(e), "\t",
                    join(map(string, st_arr), ","), "\t",
                    join(map(string, en_arr), ","), "\t",
                    L, "\n"
                ))
                close(w)
            end

            warn(LOGGER, e)
        end
        return emptyData()
    end
end

#=
Test section
=#

using Distributions
using Random
using RCall

function main()
    mu_f = 350
    sigma_f = 50
    L = 1500

    n_frag = 1000
    # f_len_arr = round.(Int, rand(Normal(mu_f, sigma_f), n_frag))

    # alpha_arr = sort([500, 800, 1000])
    # n_ats = length(alpha_arr)

    # beta_arr = [[10, 20, 30][Random.randperm(3)[1]] for _ = 1:n_ats]

    # unif_ws = 0.05
    # ws = 1 .+ rand(n_ats)
    # ws = ws ./ sum(ws)
    # ws = (1 .- unif_ws) .* ws
      
    # boarder = round.(Int, n_frag .* cumsum(ws))
    # seg_st = [1; boarder]
    # seg_en = [boarder .- 1; n_frag]
    
    # label_arr = zeros(n_frag)
    # st_arr = zeros(n_frag)
    # en_arr = zeros(n_frag)
    # for i = 1:(n_ats+1)
      
    #   tmpinds = seg_st[i]:seg_en[i] 
    #   tmpn = length(tmpinds)
      
    #   label_arr[tmpinds] .= i
      
    #   if i<=n_ats # ATS component
    #     st_arr[tmpinds] = round.(Int, rand(Normal(alpha_arr[i], beta_arr[i]), tmpn))
    #     en_arr[tmpinds] = st_arr[tmpinds] - f_len_arr[tmpinds]
    #   else # uniform component
    #     st_arr[tmpinds] = [collect(1:L)[Random.randperm(L)[1]] for _ = 1:tmpn]
    #     en_arr[tmpinds] = [collect(1:L)[Random.randperm(L)[1]] for _ = 1:tmpn]
    #   end
    # end

    # println(en_arr)
    en_arr = rcopy(R""" readRDS('test_en_arr.rds') """)
    L = 1500

    res = ATSMIX.fit_by_R(5, 1, [NaN for _ = 1:n_frag] , en_arr; L = L, mu_f=350, min_ws = 0.01, fixed_inference_flag = false, single_end_mode = true, verbose = true)

    

    println(res)
    
end

# main()