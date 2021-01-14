using DataFrames, MixedModels, RData, StatsBase, CSVFiles, RCall, DataFramesMeta, FixedEffectModels

function get_data(mat, genes, row_num, mem, MP, miss, trun, large, snv, scna, methy, fusion, peer_factors, immu, sig, cell)
    gene = genes[row_num]
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
    dat = DataFrame(sig = MP[!,sig], cell = immu[!,cell], feature = feature_dat, cancer = mem.cancer)
    rhs = 1
    if peer_factors != 0
        dat = hcat(dat, peer_factors)
        rhs = Expr(:call, :+, names(peer_factors)...)
    end
    if cancer != "panCan"
        dat = dat[dat.cancer .== cancer, :]
    end
    try
        if rhs == 1
            res = reg(dat, @formula(sig ~ (cell ~ feature) + fe(cancer)))
        else
            res = reg(dat, @eval @formula(sig ~ (cell ~ feature) + $(rhs)))
        end
        res1 = coeftable(res)
        index = res1.rownms .== "cell"
        res_dat = DataFrame(gene = mat.gene[row_num], sig = sig, cancer = cancer, feature = feature, cell = cell, Weak = res.p_kp, Est = res1.mat[index,1],
            P = res1.mat[index,4], low95 = res1.mat[index,5], up95 = res1.mat[index,6], r2 = res.adjr2)
        return res_dat
    catch
        return(DataFrame(gene = Symbol[], sig = Symbol[], cancer = Symbol[], feature = Symbol[], cell = cell, Weak = Float64[], Est = Float64[],
            P = Float64[], low95 = Float64[], up95 = Float64[], r2 = Float64[]))
    end
end

function main()
    if length(ARGS) == 0 || ARGS[1] == "-h" || ARGS[1] == "--help"
        println("usage: Gene_hypothesis2 <all_RData> <mat_file> <out_dir> <out_name>")
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
    @rput peer_num
    @rput peer_type
    @rput snv_type
    @rput split_BRCA

    R"""
    load($all_RData)

    immu <- read.table("Projects/signature_immune/The_Immune_Landscape_of_Cancer_Figure1.csv", header=T,sep="\t",stringsAsFactors=F)
    rownames(immu) = immu$TCGA.Participant.Barcode
    immu = immu[, -c(1:4, 33:36)]
    sample_int = intersect(rownames(immu), rownames(snv))
    load("/share/data4/TCGA/TCGA_fusion/fusion_mat.RData")
    int = intersect(sample_int, rownames(immu))
    mem = mem[int, ]
    immu = immu[int, ]
    immu$SNV.Neoantigens.rate = immu$SNV.Neoantigens/(mem$varSum/38)
    immu$Indel.Neoantigens.rate = immu$Indel.Neoantigens/(mem$varSum/38)
    for(i in grep("Cell", names(which(apply(immu, 2, function(x){sd(as.numeric(x), na.rm = T)}) > 5)), value = T, invert = T)){
        immu[, i] = log(immu[, i] + 1)
        colnames(immu)[which(colnames(immu)==i)] = paste0(colnames(immu)[which(colnames(immu)==i)], "_log")
    }
    MP = MP[int, ]
    if($snv_type == "count"){
        snv = snv_count
    }
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
    @rget immu
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
            for y in names(MP)
                for c in names(immu)
                    push!(res_pan_all, get_data(pan_mat, pan_genes, x, mem, MP, miss, trun, large, snv, scna, methy, fusion_mat, peer_factors, immu, y, c))
                end
            end
        end
    end
    res_pan_all = vcat(res_pan_all...)
    res_specific_all = []
    if nrow(specific_mat) > 1
        for x in 1:nrow(specific_mat)
            println(x)
            for y in names(MP)
                for c in names(immu)
                    push!(res_specific_all, get_data(specific_mat, specific_genes, x, mem, MP, miss_split, trun_split, large_split, snv, scna, methy, fusion_mat, peer_factors, immu, y, c))
                end
            end
        end
    end
    res_specific_all = vcat(res_specific_all...)
    res_all = vcat(res_pan_all, res_specific_all)
    save(File(format"TSV", joinpath(out_dir, string(out_name, ".txt.gz"))), res_all; quotechar = "", escapechar = "")
end

main()

