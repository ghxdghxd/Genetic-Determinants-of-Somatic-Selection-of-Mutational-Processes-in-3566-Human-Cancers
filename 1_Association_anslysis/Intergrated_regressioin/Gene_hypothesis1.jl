using DataFrames, MixedModels, RData, StatsBase, CSVFiles, RCall, DataFramesMeta, FixedEffectModels

function get_data(mat, genes, row_num, mem, MP, miss, trun, large, snv, scna, methy, fusion, mRNA)
    gene = genes[row_num]
    sig = Symbol(mat.sig[row_num])
    cancer = mat.cancer[row_num]
    feature = mat.feature[row_num]
    if feature == "snv"
        feature_dat = snv[!, gene]
    elseif feature == "scna"
        feature_dat = scna[!, gene]
    elseif feature == "methy"
        feature_dat = methy[!, gene]
    elseif feature == "miss"
        feature_dat = miss[!, gene]
    elseif feature == "trun"
        feature_dat = trun[!, gene]
    elseif feature == "large"
        feature_dat = large[!, gene]
    elseif feature == "fusion"
        feature_dat = fusion[!, gene]
    end
    dat = DataFrame(sig = MP[!,sig], gene = mRNA[!, gene], feature = feature_dat, cancer = mem.cancer)
    if cancer != "panCan"
        dat = dat[dat.cancer .== cancer, :]
    end
    try
        res = reg(dat, @formula(sig ~ (gene ~ feature) + fe(cancer)))
        res1 = coeftable(res)
        index = res1.rownms .== "gene"
        res_dat = DataFrame(gene = mat.gene[row_num], sig = sig, cancer = cancer, feature = feature, Weak = res.p_kp, Est = res1.mat[index,1],
            P = res1.mat[index,4], low95 = res1.mat[index,5], up95 = res1.mat[index,6], r2 = res.r2, adjr2 = res.adjr2)
        return res_dat
    catch
        return(DataFrame(gene = Symbol[], sig = Symbol[], cancer = Symbol[], feature = Symbol[], Weak = Float64[], Est = Float64[],
            P = Float64[], low95 = Float64[], up95 = Float64[], r2 = Float64[], adjr2 = Float64[]))
    end
end

function main()
    if length(ARGS) == 0 || ARGS[1] == "-h" || ARGS[1] == "--help"
        println("usage: Gene_hypothesis1.jl <all_RData> <mat_file> <out_dir> <out_name>")
        return true
    elseif length(ARGS) == 8
        all_RData = ARGS[1]
        mat_file = ARGS[2]
        out_dir = ARGS[3]
        out_name = ARGS[4]
    else
        error("the ARGS should be <all_RData> <mat_file> <out_dir> <out_name>")
    end

    mat = DataFrame(load(File(format"TSV", mat_file)))
    @rput all_RData

    R"""
    load($all_RData)
    load("/share/data4/TCGA/TCGA_The_Immune_Landscape_of_Cancer/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.log2.RData")
    load("/share/data4/TCGA/TCGA_fusion/fusion_mat.RData")
    int = intersect(rownames(snv), rownames(mRNA))
    mem = mem[int, ]
    mRNA = mRNA[int, ]
    MP = MP[int, ]
    snv = snv[int, ]
    scna = cnvT3[int, ]
    methy = methyT0[int, ]
    miss = miss[int, ]
    miss_split = miss_split[int, ]
    trun = trun[int, ]
    trun_split = trun_split[int, ]
    large = large[int, ]
    large_split = large_split[int, ]
    fusion_mat = fusion_mat[match(int, rownames(fusion_mat)),]
    rownames(fusion_mat) = int
    fusion_mat[is.na(fusion_mat)] = 0
    """

    @rget mem
    @rget MP
    @rget mRNA
    @rget snv
    @rget scna
    @rget methy
    @rget miss
    @rget trun
    @rget large
    @rget miss_split
    @rget trun_split
    @rget large_split
    @rget fusion_mat

    MP = mapcols(x -> replace(x, NaN=>missing), MP)
    mem[:cancer] = categorical(mem[:cancer])
    pan_mat = mat[mat.cancer .== "panCan", :]
    pan_genes = map(x -> Symbol(replace(x, "."=>"_")), pan_mat.gene)
    specific_mat = mat[mat.cancer .!= "panCan", :]
    specific_genes = map(x -> Symbol(replace(x, "."=>"_")), specific_mat.gene)
    res_pan_all = []
    if nrow(pan_mat) > 1
        for x in 1:nrow(pan_mat)
            println(x)
            push!(res_pan_all, get_data(pan_mat, pan_genes, x, mem, MP, miss, trun, large, snv, scna, methy, fusion_mat, mRNA))
        end
        res_pan_all = vcat(res_pan_all...)
    end
    res_specific_all = []
    if nrow(specific_mat) > 1
        for x in 1:nrow(specific_mat)
            println(x)
            push!(res_specific_all, get_data(specific_mat, specific_genes, x, mem, MP, miss_split, trun_split, large_split, snv, scna, methy, fusion_mat, mRNA))
        end
        res_specific_all = vcat(res_specific_all...)
        res_pan_all = vcat(res_pan_all, res_specific_all)
    end
    save(File(format"TSV", joinpath(out_dir, string(out_name, ".txt.gz"))), res_pan_all; quotechar = "", escapechar = "")
end

main()

