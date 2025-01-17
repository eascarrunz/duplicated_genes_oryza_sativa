using Pkg; Pkg.activate(".")
# using CSV, DataFrames
using ProgressMeter
using Dates

is_fasta_description(line) = first(line) == '>'

function beginswith(needle, haystack)
    length(needle) > length(haystack) && return false
    
    for (c1, c2) in zip(needle, haystack)
        c1 ≠ c2 && return false
    end

    return true
end

function get_fasta_id(line)
    i = 2
    id = ""

    while (i ≤ length(line))
        isspace(line[i]) && break
        id *= line[i]
        i += 1
    end

    return id
end

function get_seq_from_fasta(file)
    line = readline(file)
    @assert is_fasta_description(line)
    seq = ""

    eof(file) && return seq

    line = readline(file)
    while ! is_fasta_description(line) && isspace(first(line))
        seq *= line
        eof(file) && break
        line = readline(file)
    end

    return seq
end


"""
Return a dictionary with the line numbers of each seqid in a fasta file
"""
function fasta_seqid_dict(fasta_file)
    d = Dict{String,Int}()
    i = 1
    open(fasta_file) do file
        line = readline(file)
        while ! eof(file)
            if is_fasta_description(line)
                seqid = get_fasta_id(line)
                d[seqid] = i
            end
            line = readline(file)
            i += 1
        end
    end

    return d
end

function extract_sequence(file, linenum)
    seq = ""

    open(file) do io
        line = ""
        for _ in 1:linenum
            line = readline(io)
            @assert ! eof(io)
        end
        @assert is_fasta_description(line)

        line = readline(io)
        @assert ! is_fasta_description(line)

        while ! is_fasta_description(line) && (! isspace(first(line)))
            eof(io) && break
            seq *= line
            line = readline(io)
        end
    end

    return seq
end

function write_fasta(file, seqid, seq, mode)
    b = 0
    open(file, mode) do io
        b += write(io, '>', seqid, '\n')
        b += write(io, seq, '\n')
    end

    return b
end

function write_yn00_control_file(ctlfile, seqfile, outfile)
    b = 0
    open(ctlfile, "w") do io
        b += write(io, "seqfile = ", seqfile, '\n')
        b += write(io, "outfile = ", outfile, '\n')
        b += write(io, "verbose = 1", '\n')
        b += write(io, "icode = 0", '\n')
        b += write(io, "commonf3x4 = 0", '\n')
    end

    return b
end

function read_yn00_output(file)
    result_table_header = "seq. seq.     S       N        t   kappa   omega     dN +- SE    dS +- SE"

    line = ""
    open(file) do io
        while ! eof(io)
            line = readline(io)
            if beginswith(result_table_header, line)
                line = readline(io)
                line = isempty(line) ? readline(io) : line
                break
            end
        end
    end

    data_indices = [3:8..., 10, 11, 13]
    data = join(split(line)[data_indices], '\t')

    return data
end


function get_gene_combination_from_file(file, linelimit)
    return function get_new(c::Channel)
        open(file) do io
            line = readline(io) # Read table headers
            i = 0
            while (! eof(io)) && (i < linelimit)
                line = readline(io)
                i += 1
                family, seqid1, seqid2 = split(line, '\t')
                put!(c, (family=family, seqid1=seqid1, seqid2=seqid2))
            end
        end
    end
end

struct JobContext
    id::Int
    jdir::String    # Work directory for job
    cds_fasta_file::String
    pep_fasta_file::String
    aln_file::String
    phy_file::String
    outfile::String
    logfile::String
    logio::IO
    cmd_clustal::Base.AbstractCmd
    cmd_pal2nal::Base.AbstractCmd
    cmd_yn00::Base.AbstractCmd
    skipped::Int

    function JobContext(id)
        jdir = mktempdir()
        cds_fasta_file, pep_fasta_file, aln_file, phy_file, outfile, logfile = (mktemp(jdir)[1] for _ in 1:6)
        logio = open(logfile, "a")
        cmd_clustal = pipeline(`clustalw2 -quiet -align -infile=$(pep_fasta_file) -outfile=$(aln_file)`, stdout=devnull, stderr=logio)
        cmd_pal2nal = pipeline(`./src/pal2nal.pl $(aln_file) $(cds_fasta_file) -output paml`, stdout=phy_file, stderr=logio)
        cmd_yn00 = pipeline(setenv(`yn00`, dir=jdir), stdout=logio, stderr=logio)

        return new(id, jdir, cds_fasta_file, pep_fasta_file, aln_file, phy_file, outfile, logfile, logio, cmd_clustal, cmd_pal2nal, cmd_yn00)
    end
end


function drive_job_wrapper(input_chnl, job_context_chnl, batch_ctx)
    family, seqid1, seqid2 = take!(input_chnl)
    job_ctx = take!(job_context_chnl)

    ctx_skipped[job_ctx.id] += drive_job(family, seqid1, seqid2, job_ctx, batch_ctx)

    next!(batch_ctx.meter)
    put!(job_context_chnl, job_ctx)
end


function drive_job(family, seqid1, seqid2, job_ctx, batch_ctx)
    println(job_ctx.logio, "Processing $(seqid1) and $(seqid2)------------------------------------")
    
    if ! haskey(batch_ctx.line_nums.cds, seqid1) || ! haskey(batch_ctx.line_nums.cds, seqid2)
        println(job_ctx.logio, "Missing key error")
        return 1
    end

    linenum_cds1 = batch_ctx.line_nums.cds[seqid1]
    linenum_pep1 = batch_ctx.line_nums.pep[seqid1]
    linenum_cds2 = batch_ctx.line_nums.cds[seqid2]
    linenum_pep2 = batch_ctx.line_nums.pep[seqid2]

    # Temporary cds fasta file
    seq = extract_sequence(batch_ctx.all_cds_fasta, linenum_cds1)
    b = write_fasta(job_ctx.cds_fasta_file, seqid1, seq, "w")
    @assert b > 0
    seq = extract_sequence(batch_ctx.all_cds_fasta, linenum_cds2)
    b = write_fasta(job_ctx.cds_fasta_file, seqid2, seq, "a")
    @assert b > 0
    
    # Temporary pep fasta file
    seq = extract_sequence(batch_ctx.all_pep_fasta, linenum_pep1)
    b = write_fasta(job_ctx.pep_fasta_file, seqid1, seq, "w")
    @assert b > 0
    seq = extract_sequence(batch_ctx.all_pep_fasta, linenum_pep2)
    b = write_fasta(job_ctx.pep_fasta_file, seqid2, seq, "a")
    @assert b > 0

    run(job_ctx.cmd_clustal)

    try
        run(job_ctx.cmd_pal2nal)
    catch
        println(job_ctx.logio, "\n ^^    Error in pal2nal with $(seqid1) and $(seqid2)    ^^")
        # @lock batch_ctx.loglock println(logio, "\n ^^    Error in pal2nal with $(seqid1) and $(seqid2)    ^^")
        return 1
    end
     
    try
        run(job_ctx.cmd_yn00)
    catch
        open(job_ctx.phy_file) do f
            while ! eof(f)
                write(job_ctx.logio, "PHY FILE BELOW: ---------------\n")
                line = readline(f, keep=true)
                write(job_ctx.logio, line)
                write(job_ctx.logio, "--------------- END OF PHY FILE")
            end
        end
        println(job_ctx.logio, "\n ^^    Error in yn00 with $(seqid1) and $(seqid2)    ^^")
        # @lock batch_ctx.loglock println(logio, "\n ^^    Error in yn00 with $(seqid1) and $(seqid2)    ^^")
        return 1
    end

    result_string = read_yn00_output(job_ctx.outfile)
    @lock batch_ctx.outlock write(batch_ctx.outio, family, '\t', seqid1, '\t', seqid2, '\t', result_string, '\n')

    println(job_ctx.logio, "\n")

    return 0
end

function main()
    t0 = Dates.now()
    nt = Threads.nthreads()
    
    cds_fasta = "data/Oryza_sativa.IRGSP-1.0.cds.all.fa"
    pep_fasta = "data/Oryza_sativa.IRGSP-1.0.pep.all.fa"
    cds_line_nums = fasta_seqid_dict(cds_fasta)
    pep_line_nums = fasta_seqid_dict(pep_fasta)
    
    combinations_file = "output/gene_combinations.tsv"

    maxlines = countlines(combinations_file)
    linelimit = isempty(ARGS) ? maxlines : parse(Int, ARGS[1])
    linelimit = min(maxlines, linelimit)

    input_chnl = Channel(get_gene_combination_from_file(combinations_file, linelimit), 256)
    outfile = "output/seq_divergence.tsv"
    global ctx_skipped = zeros(Int, nt)
    final_logfile = "output/seq_divergence.log"

    open(final_logfile, "w") do logio
        println(logio, "============= GENOMIC SEQUENCE DIVERGENCE ANALYSIS =============")
        println(logio, "Started at $(t0)")
        println(logio, "Logging errors only")
        println(logio, "----------------------------------------------------------------")
    end
    
    BATCH_CONTEXT = (
        all_cds_fasta = cds_fasta,
        all_pep_fasta = pep_fasta,
        line_nums = (cds=cds_line_nums, pep=pep_line_nums),
        outio = open(outfile, "w"),
        outlock = ReentrantLock(),
        meter = Progress(linelimit; dt=1.0)
        )
    
    job_context_pool = [JobContext(i) for i in 1:nt]
    job_context_chnl = Channel{JobContext}(nt)
    for ctx in job_context_pool
        put!(job_context_chnl, ctx)
        write_yn00_control_file(joinpath(ctx.jdir, "yn00.ctl"), ctx.phy_file, ctx.outfile)
    end

    println(BATCH_CONTEXT.outio, "family\tseqid1\tseqid2\tS\tN\tt\tκ\tω\tdN\tSE_dN\tdS\tSE_dS")
    
    Threads.@threads for _ in 1:nt
        while ! isempty(input_chnl)
            drive_job_wrapper(input_chnl, job_context_chnl, BATCH_CONTEXT)
        end
    end
    
    close(BATCH_CONTEXT.outio)
    skipped = sum(ctx_skipped)
    
    t1 = Dates.now()
    open(final_logfile, "a") do io
        for ctx in job_context_pool
            close(ctx.logio)
            open(ctx.logfile) do cio
                while ! eof(cio)
                    line = readline(cio, keep=true)
                    write(io, line)
                end
            end
        end
    println(io, "------------- GENOMIC SEQUENCE DIVERGENCE ANALYSIS -------------")
    println(io, "Finished at $(t1) - Duration: $(t1 - t0)")
    println(io, "Comparisons skipped: $(skipped) ($(round(100 * skipped/linelimit; digits=1)) %)")
    println(io, "================================================================")
    end

    finish!(BATCH_CONTEXT.meter)

    # println("Skipped: $(skipped)")

    return 0
end

main()
