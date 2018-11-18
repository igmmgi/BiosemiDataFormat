using BioSemiBDF
using Test

@testset "Newtest17-256_read_bdf" begin

    bdf_filename = joinpath(dirname(@__FILE__), "Newtest17-256.bdf")
    dat = read_bdf(bdf_filename)
    
    @test isequal(dat.header["num_bytes_header"], 4608)
    @test isequal(dat.header["num_channels"], 17)
    @test isequal(dat.header["num_data_records"], 60)
    @test isequal(dat.header["sample_rate"][1], 256)
    @test isequal(size(dat.data), (16, 15360))
    @test isequal(dat.triggers["idx"][1], 415)
    @test isequal(dat.triggers["val"][1], 255)
    @test isequal(dat.triggers["count"][255], 40)

end

@testset "Newtest17-256_read_bdf_channels" begin

    bdf_filename = joinpath(dirname(@__FILE__), "Newtest17-256.bdf")
    dat1 = read_bdf(bdf_filename)
    dat2 = read_bdf(bdf_filename, channels=[1, 3, 5])

    @test isequal(dat1.data[1,:], dat2.data[1,:])
    @test isequal(dat1.data[5,:], dat2.data[3,:])
    @test isequal(dat2.header["num_channels"], 4)
    @test isequal(dat2.header["num_data_records"], 60)
    @test isequal(dat2.header["sample_rate"][1], 256)
    @test isequal(size(dat2.data), (3, 15360))
    @test isequal(dat2.triggers["idx"][1], 415)
    @test isequal(dat2.triggers["val"][1], 255)
    @test isequal(dat2.triggers["count"][255], 40)

end

@testset "Newtest17-256_write_bdf" begin

    bdf_filename = joinpath(dirname(@__FILE__), "Newtest17-256.bdf")
    dat1 = read_bdf(bdf_filename)
    dat1.header["filename"] = "Newtest17-256_write.bdf"
    write_bdf(dat1)

    bdf_filename = joinpath(dirname(@__FILE__), "Newtest17-256_write.bdf")
    dat2 = read_bdf(bdf_filename)

    @test isequal(dat1.data, dat2.data)
    @test isequal(dat1.status, dat2.status)
    @test isequal(dat1.time, dat2.time)
    @test isequal(dat1.triggers, dat2.triggers)

    rm(bdf_filename)

end

@testset "Newtest17-2048_read_bdf" begin

    bdf_filename = joinpath(dirname(@__FILE__), "Newtest17-2048.bdf")
    dat = read_bdf(bdf_filename)

    @test isequal(dat.header["num_bytes_header"], 4608)
    @test isequal(dat.header["num_channels"], 17)
    @test isequal(dat.header["num_data_records"], 60)
    @test isequal(dat.header["sample_rate"][1], 2048)
    @test isequal(size(dat.data), (16, 122880))
    @test isequal(dat.triggers["idx"][1], 3353)
    @test isequal(dat.triggers["val"][1], 255)
    @test isequal(dat.triggers["count"][255], 39)

end

@testset "Newtest17-2048_read_bdf_channels" begin

    bdf_filename = joinpath(dirname(@__FILE__), "Newtest17-2048.bdf")
    dat1 = read_bdf(bdf_filename)
    dat2 = read_bdf(bdf_filename, channels=["A1", "A3", "A5"])

    @test isequal(dat1.data[1,:], dat2.data[1,:])
    @test isequal(dat1.data[5,:], dat2.data[3,:])
    @test isequal(dat2.header["num_channels"], 4)
    @test isequal(dat2.header["num_data_records"], 60)
    @test isequal(dat2.header["sample_rate"][1], 2048)
    @test isequal(size(dat2.data), (3, 122880))
    @test isequal(dat2.triggers["idx"][1], 3353)
    @test isequal(dat2.triggers["val"][1], 255)
    @test isequal(dat2.triggers["count"][255], 39)

end

@testset "Newtest17-2048_write_bdf" begin

    bdf_filename = joinpath(dirname(@__FILE__), "Newtest17-2048.bdf")
    dat1 = read_bdf(bdf_filename)
    dat1.header["filename"] = "Newtest17-2048_write.bdf"
    write_bdf(dat1)

    bdf_filename = joinpath(dirname(@__FILE__), "Newtest17-2048_write.bdf")
    dat2 = read_bdf(bdf_filename)

    @test isequal(dat1.data, dat2.data)
    @test isequal(dat1.status, dat2.status)
    @test isequal(dat1.time, dat2.time)
    @test isequal(dat1.triggers, dat2.triggers)

    rm(bdf_filename)

end

@testset "Newtest17-256_merge_bdf" begin

    bdf_filename = joinpath(dirname(@__FILE__), "Newtest17-256.bdf")
    dat1 = read_bdf(bdf_filename)
    dat2 = read_bdf(bdf_filename)
    dat3 = merge_bdf([dat1, dat2], "merged_datafile")

    @test isequal(dat3.header["num_bytes_header"], 4608)
    @test isequal(dat3.header["num_channels"], 17)
    @test isequal(dat3.header["num_data_records"], 120)
    @test isequal(dat3.header["sample_rate"][1], 256)
    @test isequal(size(dat3.data), (16, 15360*2))
    @test isequal(dat3.triggers["idx"][1], 415)
    @test isequal(dat3.triggers["val"][1], 255)
    @test isequal(dat3.triggers["count"][255], 80)

end

@testset "Newtest17-2048_merge_bdf" begin

    bdf_filename = joinpath(dirname(@__FILE__), "Newtest17-2048.bdf")
    dat1 = read_bdf(bdf_filename)
    dat2 = read_bdf(bdf_filename)
    dat3 = merge_bdf([dat1, dat2], "merged_datafile")

    @test isequal(dat3.header["num_bytes_header"], 4608)
    @test isequal(dat3.header["num_channels"], 17)
    @test isequal(dat3.header["num_data_records"], 120)
    @test isequal(dat3.header["sample_rate"][1], 2048)
    @test isequal(size(dat3.data), (16, 122880*2))
    @test isequal(dat3.triggers["idx"][1], 3353)
    @test isequal(dat3.triggers["val"][1], 255)
    @test isequal(dat3.triggers["count"][255], 78)

end

@testset "Newtest17-256_select_channels_bdf" begin

    bdf_filename = joinpath(dirname(@__FILE__), "Newtest17-256.bdf")
    dat1 = read_bdf(bdf_filename)
    dat2 = select_channels_bdf(dat1, [1, 3, 5])

    @test isequal(dat1.data[1,:], dat2.data[1,:])
    @test isequal(dat1.data[5,:], dat2.data[3,:])
    @test isequal(dat2.header["num_channels"], 4)
    @test isequal(dat2.header["num_data_records"], 60)
    @test isequal(dat2.header["sample_rate"][1], 256)
    @test isequal(size(dat2.data), (3, 15360))
    @test isequal(dat2.triggers["idx"][1], 415)
    @test isequal(dat2.triggers["val"][1], 255)
    @test isequal(dat2.triggers["count"][255], 40)

end

@testset "Newtest17-2048_select_channels_bdf" begin

    bdf_filename = joinpath(dirname(@__FILE__), "Newtest17-2048.bdf")
    dat1 = read_bdf(bdf_filename)
    dat2 = select_channels_bdf(dat1, ["A1", "A3", "A5"])

    @test isequal(dat1.data[1,:], dat2.data[1,:])
    @test isequal(dat1.data[5,:], dat2.data[3,:])
    @test isequal(dat2.header["num_channels"], 4)
    @test isequal(dat2.header["num_data_records"], 60)
    @test isequal(dat2.header["sample_rate"][1], 2048)
    @test isequal(size(dat2.data), (3, 122880))
    @test isequal(dat2.triggers["idx"][1], 3353)
    @test isequal(dat2.triggers["val"][1], 255)
    @test isequal(dat2.triggers["count"][255], 39)

end

@testset "Newtest17-256_crop_bdf" begin

    bdf_filename = joinpath(dirname(@__FILE__), "Newtest17-256.bdf")
    dat1 = read_bdf(bdf_filename)
    dat2 = crop_bdf(dat1, "records", [10 20], "Newtest17-256_cropped.bdf")

    @test isequal(size(dat1.data,  1), size(dat2.data, 1))
    @test !isequal(size(dat1.data, 2), size(dat2.data, 2))
    @test isequal(dat2.header["num_data_records"], 11)

end

@testset "Newtest17-2048_crop_bdf" begin

    bdf_filename = joinpath(dirname(@__FILE__), "Newtest17-2048.bdf")
    dat1 = read_bdf(bdf_filename)
    dat2 = crop_bdf(dat1, "records", [10 20], "Newtest17-2048_cropped.bdf")

    @test isequal(size(dat1.data,  1), size(dat2.data, 1))
    @test !isequal(size(dat1.data, 2), size(dat2.data, 2))
    @test isequal(dat2.header["num_data_records"], 11)

end


@testset "Newtest17-256_downsample_bdf" begin

    bdf_filename = joinpath(dirname(@__FILE__), "Newtest17-256.bdf")
    dat1 = read_bdf(bdf_filename)
    dat2 = downsample_bdf(dat1, 2, "Newtest17-128.bdf")

    @test isequal(dat2.header["sample_rate"][1], 128)

end

@testset "Newtest17-2048_downsample_bdf" begin

    bdf_filename = joinpath(dirname(@__FILE__), "Newtest17-2048.bdf")
    dat1 = read_bdf(bdf_filename)
    dat2 = downsample_bdf(dat1, 2, "Newtest17-1024.bdf")

    @test isequal(dat2.header["sample_rate"][1], 1024)

end

@testset "channel_idx" begin

    x = BioSemiBDF.channel_idx(["A1", "A2", "A3"], ["A1", "A3"])
    @test isequal(x, [1, 3])

    x = BioSemiBDF.channel_idx(["A1", "A2", "A3", "A4"], ["A1", "A4"])
    @test isequal(x, [1, 4])

    x = BioSemiBDF.channel_idx(["A1", "A2", "A3", "A4"], ["A1", "A4", "A4"])
    @test isequal(x, [1, 4])

    @test_throws(ErrorException, BioSemiBDF.channel_idx(["A1"], ["zzz"]))

    x = BioSemiBDF.channel_idx(["A1", "A2", "A3"], [1, 3])
    @test isequal(x, [1, 3])

    x = BioSemiBDF.channel_idx(["A1", "A2", "A3", "A4"], [1, 4])
    @test isequal(x, [1, 4])

    x = BioSemiBDF.channel_idx(["A1", "A2", "A3", "A4"], [1, 4, 4])
    @test isequal(x, [1, 4])

    @test_throws(ErrorException, BioSemiBDF.channel_idx(["A1"], [2]))
    @test_throws(ErrorException, BioSemiBDF.channel_idx(["A1"], [1, 2]))

end
