using DataFrames, MixedModels, RData, StatsBase, Distributed, CSV, RCall, DataFramesMeta, FixedEffectModels, GZip
# addprocs(60)
@everywhere using DataFrames, MixedModels, RData, StatsBase, Distributed, CSV, RCall, DataFramesMeta, FixedEffectModels, GZip

@everywhere function ivreg(cancer, sig, immu_cell, snp, data)
    # println(gene)
    dat = @linq data[!, ["cancer", sig, immu_cell, snp]] |>
        rename([Symbol(sig) => :sig, Symbol(immu_cell) =>:immu_cell, Symbol(snp)=>:snp])
    try
        if cancer == "panCan"
            # res = reg(dat, @model(sig ~ (immu_cell ~ snp), fe = cancer))
            res = reg(dat, @formula(sig ~ (immu_cell ~ snp) + fe(cancer)))
        else
            # res = reg(dat, @model(sig ~ (gene ~ snp)))
            res = reg(dat, @formula(sig ~ (immu_cell ~ snp)))
        end
        res1 = coeftable(res)
        index = res1.rownms .== "immu_cell"
        res_dat = DataFrame(immu_cell = immu_cell, Weak = res.p_kp, Est = res1.mat[index,1], P = res1.mat[index,4], low95 = res1.mat[index,5], up95 = res1.mat[index,6])
        return res_dat
    catch
        return(DataFrame(immu_cell = Symbol[], Weak = Float64[], Est = Float64[], P = Float64[], low95 = Float64[], up95 = Float64[]))
    end
end

@everywhere function get_data(mat, row_num, data, immu_cell)
    out = map(x -> ivreg(mat.cancer[row_num], mat.sig[row_num], x, mat.ID[row_num], data), immu_cell)
    out = vcat(out...)
    out[!, "cancer"] .= mat.cancer[row_num]
    out[!, "sig"] .= mat.sig[row_num]
    out[!, "snp"] .= mat.ID[row_num]
    return(out)
end

function main()
    if length(ARGS) == 0 || ARGS[1] == "-h" || ARGS[1] == "--help"
        println("usage: mpQTLs_hypothesis2.jl chrom cancer out_dir mat_file")
        return true
    elseif length(ARGS) == 4
        chrom = ARGS[1] # string
        cancer = ARGS[2]
        out_dir = ARGS[3]
        mat_file = ARGS[4]
    else
        error("the ARGS should be <chrom> <cancer> <out_dir> <mat_file>")
    end
    # mat = CSV.read("Projects/TCGA/hg19/workflow/09com/julia/sig.mat.txt.1", delim="\t")
    mat = CSV.read(mat_file, delim = "\t")
    mat.chrom = [split(x, r"[p|q]")[1] for x in mat.cytoband]
    # mat2 = @linq mat |> by([:cancer,:chrom], num = length(:ID))
    mat = @where(mat, :chrom .== string(chrom), :cancer .== cancer)
    # cancers = @where(mat2, :chrom .== chrom)[:cancer]

    @rput mat
    @rput chrom
    @rput cancer

    R"""
    if(cancer == "panCan"){
        geno_file = paste0("Projects/TCGA/hg19/workflow/09com/julia/input_geno_0.9/panCan/", chrom, ".geno.RData")
    }else{
        geno_file = list.files(paste0("Projects/TCGA/hg19/workflow/09com/julia/input_geno_0.9/", cancer),
            pattern=paste0("*chr", chrom, ".COM.geno.RData"), full.names = T)
    }
    load(geno_file)
    immu <- read.table("Projects/signature_immune/The_Immune_Landscape_of_Cancer_Figure1.csv", header=T,sep="\t",stringsAsFactors=F)
    rownames(immu) = immu$TCGA.Participant.Barcode
    immu = immu[, -c(1:4, 33:36)]

    load("Projects/TCGA/hg19/workflow/01maf/mem.RData")
    int = intersect(rownames(geno), rownames(immu))
    mem = mem[int, ]
    mem$cancer[grep('BRCA', mem$cancer)] = "BRCA"
    logSigMat = logSigMat[int, ]
    immu = immu[int, ]
    immu$SNV.Neoantigens.rate = immu$SNV.Neoantigens/(mem$varSum/38)
    immu$Indel.Neoantigens.rate = immu$Indel.Neoantigens/(mem$varSum/38)
    for(i in grep("Cell", names(which(apply(immu, 2, function(x){sd(as.numeric(x), na.rm = T)}) > 5)), value = T, invert = T)){
        immu[, i] = log(immu[, i] + 1)
        colnames(immu)[which(colnames(immu)==i)] = paste0(colnames(immu)[which(colnames(immu)==i)], "_log")
    }
    rsid = unique(mat$ID)
    geno = as.data.frame(geno[int, rsid])
    colnames(geno) = rsid
    """

    @rget mem
    @rget logSigMat
    @rget immu
    @rget geno

    logSigMat = mapcols(x -> replace(x, NaN=>missing), logSigMat)
    data = hcat(mem[!, ["sample","cancer"]], logSigMat, immu, geno)
    data[!, "cancer"] = categorical(data[!, "cancer"])

    res_all = []
    for x in 1:nrow(mat)
        # push!(res_all, @spawn get_data(mat, x, data, genes_mat))
        push!(res_all, get_data(mat, x, data, names(immu)))
    end
    # res_all = pmap(x -> get_data(mat, x, data, genes_mat), 1:nrow(mat))
    # res_all = map(fetch, res_all)
    res_all = vcat(res_all...)
    GZip.open(joinpath(out_dir, string(cancer, ".", chrom, ".IV.txt.gz")), "w") do io
        CSV.write(io, res_all, delim = "\t")
    end
end

main()


# for x in 1:nrow(mat)
#     println(x)
#     push!(res_all, get_data(mat, x, data, genes_mat))
# end