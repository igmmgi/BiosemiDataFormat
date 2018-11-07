using BioSemiBDF
using Test

@testset "Newtest17-256_read_bdf" begin

    bdf_filename = joinpath(dirname(@__FILE__), "Newtest17-256.bdf")
    dat = read_bdf(bdf_filename)

    # some header details
    @test isequal(dat.header["num_bytes_header"], 4608)
    @test isequal(dat.header["num_channels"], 17)
    @test isequal(dat.header["num_data_records"], 60)
    @test isequal(dat.header["sample_rate"][1], 256)

    # some data
    @test isequal(size(dat.data), (16, 15360))

    # channel labels
    @test isequal(dat.header["channel_labels"][1], dat.labels[1])
    @test isequal(dat.header["channel_labels"][16], dat.labels[16])

    # triggers
    @test isequal(dat.triggers["idx"][1], 415)
    @test isequal(dat.triggers["val"][1], 255)
    @test isequal(dat.triggers["count"][255], 40)

end

@testset "Newtest17-256_write_bdf" begin

    bdf_filename = joinpath(dirname(@__FILE__), "Newtest17-256.bdf")
    dat1 = read_bdf(bdf_filename)
    dat1.header["filename"] = "Newtest17-256_write.bdf"
    write_bdf(dat1)

    bdf_filename = joinpath(dirname(@__FILE__), "Newtest17-256_write.bdf")
    dat2 = read_bdf(bdf_filename)

    # some header details
    @test isequal(dat1.data, dat2.data)
    @test isequal(dat1.labels, dat2.labels)
    @test isequal(dat1.status, dat2.status)
    @test isequal(dat1.time, dat2.time)
    @test isequal(dat1.triggers, dat2.triggers)

    rm(bdf_filename)

end

@testset "Newtest17-2048_write_bdf" begin

    bdf_filename = joinpath(dirname(@__FILE__), "Newtest17-2048.bdf")
    dat1 = read_bdf(bdf_filename)
    dat1.header["filename"] = "Newtest17-2048_write.bdf"
    write_bdf(dat1)

    bdf_filename = joinpath(dirname(@__FILE__), "Newtest17-2048_write.bdf")
    dat2 = read_bdf(bdf_filename)

    # some header details
    @test isequal(dat1.data, dat2.data)
    @test isequal(dat1.labels, dat2.labels)
    @test isequal(dat1.status, dat2.status)
    @test isequal(dat1.time, dat2.time)
    @test isequal(dat1.triggers, dat2.triggers)

    rm(bdf_filename)

end


@testset "Newtest17-2048_read_bdf" begin

    bdf_filename = joinpath(dirname(@__FILE__), "Newtest17-2048.bdf")
    dat = read_bdf(bdf_filename)

    # some header details
    @test isequal(dat.header["num_bytes_header"], 4608)
    @test isequal(dat.header["num_channels"], 17)
    @test isequal(dat.header["num_data_records"], 60)
    @test isequal(dat.header["sample_rate"][1], 2048)

    # some data
    @test isequal(size(dat.data), (16, 122880))

    # channel labels
    @test isequal(dat.header["channel_labels"][1], dat.labels[1])
    @test isequal(dat.header["channel_labels"][16], dat.labels[16])

    # triggers
    @test isequal(dat.triggers["idx"][1], 3353)
    @test isequal(dat.triggers["val"][1], 255)
    @test isequal(dat.triggers["count"][255], 39)

end


@testset "Newtest17-256_merge_bdf" begin

    bdf_filename = joinpath(dirname(@__FILE__), "Newtest17-256.bdf")
    dat1 = read_bdf(bdf_filename)
    dat2 = read_bdf(bdf_filename)

    dat3 = merge_bdf([dat1, dat2], "merged_datafile")

    # some header details
    @test isequal(dat3.header["num_bytes_header"], 4608)
    @test isequal(dat3.header["num_channels"], 17)
    @test isequal(dat3.header["num_data_records"], 120)
    @test isequal(dat3.header["sample_rate"][1], 256)

    # some data
    @test isequal(size(dat3.data), (16, 15360*2))

    # channel labels
    @test isequal(dat3.header["channel_labels"][1], dat3.labels[1])
    @test isequal(dat3.header["channel_labels"][16], dat3.labels[16])

    # triggers
    @test isequal(dat3.triggers["idx"][1], 415)
    @test isequal(dat3.triggers["val"][1], 255)
    @test isequal(dat3.triggers["count"][255], 80)

end

@testset "Newtest17-2048_merge_bdf" begin

    bdf_filename = joinpath(dirname(@__FILE__), "Newtest17-2048.bdf")
    dat1 = read_bdf(bdf_filename)
    dat2 = read_bdf(bdf_filename)

    dat3 = merge_bdf([dat1, dat2], "merged_datafile")

    # some header details
    @test isequal(dat3.header["num_bytes_header"], 4608)
    @test isequal(dat3.header["num_channels"], 17)
    @test isequal(dat3.header["num_data_records"], 120)
    @test isequal(dat3.header["sample_rate"][1], 2048)

    # some data
    @test isequal(size(dat3.data), (16, 122880*2))

    # channel labels
    @test isequal(dat3.header["channel_labels"][1], dat3.labels[1])
    @test isequal(dat3.header["channel_labels"][16], dat3.labels[16])

    # triggers
    @test isequal(dat3.triggers["idx"][1], 3353)
    @test isequal(dat3.triggers["val"][1], 255)
    @test isequal(dat3.triggers["count"][255], 78)

end
