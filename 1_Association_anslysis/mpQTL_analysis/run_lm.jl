using DataFrames, MixedModels, RData, StatsBase, Distributed, CSV, RCall, GLM

function lm1(data)
    try
        rhs = Expr(:call, :+, names(data)[2:end]...)
        fit1 = lm(@eval(@formula(y ~ $(rhs))), data)
        res1 = coeftable(fit1)
        res1 = DataFrame(Est = res1.cols[1][2], P = res1.cols[4][2])
        res1.cancers = ["panCancer"]
        return(res1)
    catch
        return(DataFrame(Est = Float64[], P = Float64[], cancers = String[]))
    end
end

function lm2(data)
    try
        rhs = Expr(:call, :+, names(data)[4:end]...)
        fit2 = lm(@eval(@formula(y ~ cancer + x&cancer + $(rhs))), data)
        res2 = coeftable(fit2)
        res2 = DataFrame(Est = res2.cols[1][end-13:end], P = res2.cols[4][end-13:end])
        res2.cancers = ["BRCA","COAD","ESCA","ESCC","GBM","KIRC","LIHC","LUAD","OV","PRAD","STAD","THCA","UCEC"]
        return(res2)
    catch
        return(DataFrame(Est = Float64[], P = Float64[], cancers = String[]))
    end
end

function run_lm(sig, snp, geno, MP, batch)
    println(snp)
    data = hcat(DataFrame(y = MP[Symbol(sig)], x = geno[snp], cancer = batch.cancer))
    res1 = lm1(data)
    res2 = lm2(data)
    res = vcat(res1, res2)
    res.sig = sig
    res.snp = snp
    return(res)
end

function main()
    if length(ARGS) == 0 || ARGS[1] == "-h" || ARGS[1] == "--help"
        println("usage: run_lm.jl chrom sig out_dir")
        return true
    elseif length(ARGS) == 4
        chrom = ARGS[1]
        sig = ARGS[2]
        out_dir = ARGS[4]
    else
        error("the ARGS should be <chrom> <sig> <out_dir>")
    end

    R"""
    load(paste0("Projects/TCGA/hg19/workflow/09com/julia/input_geno/", $chrom, ".geno.RData"))
    """

    @rget geno
    @rget MP
    @rget batch

    MP = map(x->replace(x, NaN=>missing), eachcol(MP))
    res = map(x -> run_lm(sig, x, geno, MP, batch), names(geno))
    res = vcat(res...)
    res.snp = [replace(string(x), "x" => "") for x in res.snp]
    res = res[:, [:cancers,:sig,:snp,:Est,:P]]
    CSV.write(joinpath(out_dir, string(chrom, ".", sig, ".txt")), res, delim = "\t")
end

main()