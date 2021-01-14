using DataFrames, MixedModels, RData, StatsBase, Distributed, CSV, RCall, DataFramesMeta, FixedEffectModels, GZip
# addprocs(60)
@everywhere using DataFrames, MixedModels, RData, StatsBase, Distributed, CSV, RCall, DataFramesMeta, FixedEffectModels, GZip

@everywhere function ivreg(cancer, sig, gene, snp, data)
    # println(gene)
    dat = @linq data[[:cancer, Symbol(sig), gene, Symbol(snp)]] |>
        rename([Symbol(sig) => :sig, gene =>:gene, Symbol(snp)=>:snp])
    try
        if cancer == "panCan"
            res = reg(dat, @model(sig ~ (gene ~ snp), fe = cancer))
        else
            res = reg(dat, @model(sig ~ (gene ~ snp)))
        end
        res1 = coeftable(res)
        index = res1.rownms .== "gene"
        res_dat = DataFrame(gene = gene, Weak = res.p_kp, Est = res1.mat[index,1], P = res1.mat[index,4], low95 = res1.mat[index,5], up95 = res1.mat[index,6])
        return res_dat
    catch
        return(DataFrame(gene = Symbol[], Weak = Float64[], Est = Float64[], P = Float64[], low95 = Float64[], up95 = Float64[]))
    end
end

@everywhere function get_data(mat, row_num, data, genes_mat)
    gene_list = genes_mat.Column2[genes_mat.Column1 .== mat.ID[row_num]]
    out = map(x -> ivreg(mat.cancer[row_num], mat.sig[row_num], x, mat.ID[row_num], data), gene_list)
    out = vcat(out...)
    out.cancer = mat.cancer[row_num]
    out.sig = mat.sig[row_num]
    out.snp = mat.ID[row_num]
    return(out)
end

function main()
    if length(ARGS) == 0 || ARGS[1] == "-h" || ARGS[1] == "--help"
        println("usage: mpQTLs_hypothesis2.jl chrom out_dir mat_file gene_file")
        return true
    elseif length(ARGS) == 5
        chrom = ARGS[1] # string
        cancer = ARGS[2]
        out_dir = ARGS[3]
        mat_file = ARGS[4]
        gene_file = ARGS[5]
    else
        error("the ARGS should be <chrom> <cancer> <out_dir> <mat_file> <gene_file>")
    end
    # mat = CSV.read("Projects/TCGA/hg19/workflow/09com/julia/sig.mat.txt.1", delim="\t")
    mat = CSV.read(mat_file, delim = "\t")
    mat.chrom = [split(x, r"[p|q]")[1] for x in mat.cytoband]
    # mat2 = @linq mat |> by([:cancer,:chrom], num = length(:ID))
    mat = @where(mat, :chrom .== string(chrom), :cancer .== cancer)
    # cancers = @where(mat2, :chrom .== chrom)[:cancer]
    genes_mat = CSV.read(gene_file, delim = "\t", datarow = 1)

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
    # load("/share/data4/TCGA/TCGA_The_Immune_Landscape_of_Cancer/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.RData")
    load("Projects/TCGA/hg19/workflow/04CNV_methy_mRNA/imm/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.log2.RData")
    load("Projects/TCGA/hg19/workflow/01maf/mem.RData")
    rownames(mRNA) = substr(rownames(mRNA), 1, 12)
    int = intersect(rownames(geno), rownames(mRNA))
    mem = mem[int, ]
    mem$cancer[grep('BRCA', mem$cancer)] = "BRCA"
    logSigMat = logSigMat[int, ]
    mRNA = mRNA[int, ]
    rsid = unique(mat$ID)
    geno = as.data.frame(geno[int, rsid])
    colnames(geno) = rsid
    """

    @rget mem
    @rget logSigMat
    @rget mRNA
    @rget geno

    logSigMat = mapcols(x -> replace(x, NaN=>missing), logSigMat)
    data = hcat(mem[[:sample,:cancer]], logSigMat, mRNA, geno)
    data[:cancer] = categorical(data[:cancer])
    genes_mat.Column2 = Symbol.(genes_mat.Column2)
    genes_mat = genes_mat[map(x -> x in names(mRNA), genes_mat.Column2), :]

    int_ID = intersect(genes_mat.Column1, mat.ID)

    genes_mat = genes_mat[map(x -> x in int_ID, genes_mat.Column1), :]
    mat = mat[map(x -> x in int_ID, mat.ID), :]

    res_all = []
    for x in 1:nrow(mat)
        # push!(res_all, @spawn get_data(mat, x, data, genes_mat))
        push!(res_all, get_data(mat, x, data, genes_mat))
    end
    res_all = pmap(x -> get_data(mat, x, data, genes_mat), 1:nrow(mat))
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