using VGT
using Test

@testset "interrupt detection covers direct and threaded exceptions" begin
    @test VGT.is_interrupt_exception(InterruptException())
    @test VGT.is_interrupt_exception(CompositeException([InterruptException()]))
    @test !VGT.is_interrupt_exception(ErrorException("ordinary failure"))

    task = Threads.@spawn throw(InterruptException())
    err = try
        fetch(task)
        nothing
    catch err
        err
    end

    @test err isa TaskFailedException
    @test VGT.is_interrupt_exception(err)
end

@testset "interrupt cleanup removes partial chunks and only new empty output dirs" begin
    mktempdir() do dir
        output_dir = joinpath(dir, "run-output")
        mkpath(output_dir)
        chunk = joinpath(output_dir, "chunk_2d_1.jld2")
        tmp_chunk = joinpath(output_dir, "chunk_2d_2.jld2.tmp.123.1")
        unrelated = joinpath(output_dir, "notes.txt")
        write(chunk, "partial")
        write(tmp_chunk, "partial")
        write(unrelated, "keep")

        cleanup = VGT.cleanup_interrupted_run!(output_dir; remove_created_output_dir=true)

        @test !isfile(chunk)
        @test !isfile(tmp_chunk)
        @test isfile(unrelated)
        @test !cleanup.removed_output_dir
        @test sort(basename.(cleanup.removed_chunks)) == ["chunk_2d_1.jld2", "chunk_2d_2.jld2.tmp.123.1"]
    end

    mktempdir() do dir
        output_dir = joinpath(dir, "created-run")
        mkpath(output_dir)
        write(joinpath(output_dir, "chunk_2d_1.jld2"), "partial")

        cleanup = VGT.cleanup_interrupted_run!(output_dir; remove_created_output_dir=true)

        @test cleanup.removed_output_dir
        @test !isdir(output_dir)
    end

    mktempdir() do dir
        output_dir = joinpath(dir, "existing-run")
        mkpath(output_dir)
        write(joinpath(output_dir, "chunk_2d_1.jld2"), "partial")

        cleanup = VGT.cleanup_interrupted_run!(output_dir; remove_created_output_dir=false)

        @test !cleanup.removed_output_dir
        @test isdir(output_dir)
        @test isempty(readdir(output_dir))
    end
end
