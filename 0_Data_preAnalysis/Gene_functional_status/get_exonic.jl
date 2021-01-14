# run_exonic.jl
using CSVFiles, Grep, DataFrames

anno_head = DataFrame(load(File(format"TSV", "anno.head.txt"), header_exists = false))
samples = DataFrame(load(File(format"TSV", "samples.txt"), header_exists = false))

function filter_variant(chrom, out_dir)
    anno = DataFrame(load(File(format"TSV", string(chrom, ".hg19_multianno.txt.gz")), header_exists = false, skiplines_begin = 1, colparsers = Dict(40=>String)))
    names!(anno, Symbol.(anno_head[1,:]))
    exon_var = ["nonsynonymous SNV", "frameshift deletion", "nonframeshift insertion", "stoploss", "frameshift insertion", "nonframeshift deletion", "stopgain"]
    anno = anno[map(x -> x in exon_var, anno.ExonicFunc_refGene), :]
    anno = anno[:, Symbol.(["Chr", "Start", "End", "Ref", "Alt", "Gene_refGene", "ExonicFunc_refGene", "AAChange_refGene", "avsnp150", "ExAC_ALL",
                        "ExAC_AFR", "ExAC_AMR", "ExAC_EAS", "ExAC_FIN", "ExAC_NFE", "ExAC_OTH", "ExAC_SAS", "POS", "REF", "ALT"])]
    anno = anno[map(x -> !occursin(r"^OR", x), anno.Gene_refGene), :]
    anno = anno[map(x -> !occursin(r"^HLA", x), anno.Gene_refGene), :]
    anno_SNPs = unique(map(x -> join(x, ":"), zip(anno.Chr, anno.POS, anno.REF, anno.ALT)))

    vcf = DataFrame(load(File(format"TSV", string(chrom, ".vcf.gz")), header_exists = false))
    vcf = vcf[:,[1:2;4:5]]
    names!(vcf, Symbol.(["CHROM", "POS", "REF", "ALT"]))
    vcf_SNPs = map(x -> join(x, ":"), zip(vcf.CHROM, vcf.POS, vcf.REF, vcf.ALT))
    index = map(x -> x in anno_SNPs, vcf_SNPs)
    vcf = vcf[index, :]

    geno = DataFrame(load(File(format"CSV", string(chrom, ".geno.gz")), header_exists = false))
    geno = geno[index, 1]
    geno = split.(geno[:,1], "\t", limit = 5)
    vcf.geno = map(x -> x[5], geno)
    rename!(vcf, Dict(:geno => Symbol(join(samples[:,1], "\t"))))
    save(File(format"CSV", joinpath(out_dir, string(chrom, ".exonic.anno.gz"))), anno; delim='\t', quotechar = "", escapechar = "")
    save(File(format"CSV", joinpath(out_dir, string(chrom, ".exonic.geno.gz"))), vcf; delim='\t', quotechar = "", escapechar = "")
end

if length(ARGS) == 0 || ARGS[1] == "-h" || ARGS[1] == "--help"
    println("usage: run_lme.jl chrom out_dir")
    return true
elseif length(ARGS) == 2
    chrom = ARGS[1]
    out_dir = ARGS[2]
else
    error("the ARGS should be <chrom> <out_dir>")
end

filter_variant(chrom, out_dir)
