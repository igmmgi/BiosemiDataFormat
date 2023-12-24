using BioSemiBDF
using Test

@testset "read_bdf" begin

  # header only
  bdf_file = joinpath(dirname(@__FILE__), "Newtest17-256.bdf")
  hdr = read_bdf(bdf_file, header_only=true)
  @test isa(hdr, BioSemiBDF.BioSemiHeader)
  @test isequal(hdr.num_bytes_header, 18 * 256)
  @test isequal(hdr.num_channels, 17)
  @test isequal(hdr.num_data_records, 60)
  @test isequal(hdr.sample_rate[1], 256)

  bdf_file = joinpath(dirname(@__FILE__), "Newtest17-2048.bdf")
  hdr = read_bdf(bdf_file, header_only=true)
  @test isa(hdr, BioSemiBDF.BioSemiHeader)
  @test isequal(hdr.num_bytes_header, 18 * 256)
  @test isequal(hdr.num_channels, 17)
  @test isequal(hdr.num_data_records, 60)
  @test isequal(hdr.sample_rate[1], 2048)

  # full files
  bdf_file = joinpath(dirname(@__FILE__), "Newtest17-256.bdf")
  dat = read_bdf(bdf_file)

  @test isequal(dat.header.num_bytes_header, 18 * 256)
  @test isequal(dat.header.num_channels, 17)
  @test isequal(dat.header.num_data_records, 60)
  @test isequal(dat.header.sample_rate[1], 256)
  @test isequal(size(dat.data), (15360, 16))
  @test isequal(dat.triggers.idx[1], 415)
  @test isequal(dat.triggers.val[1], 255)
  @test isequal(dat.triggers.count[255], 40)

  bdf_file = joinpath(dirname(@__FILE__), "Newtest17-2048.bdf")
  dat = read_bdf(bdf_file)

  @test isequal(dat.header.num_bytes_header, 18 * 256)
  @test isequal(dat.header.num_channels, 17)
  @test isequal(dat.header.num_data_records, 60)
  @test isequal(dat.header.sample_rate[1], 2048)
  @test isequal(size(dat.data), (122880, 16))
  @test isequal(dat.triggers.idx[1], 3353)
  @test isequal(dat.triggers.val[1], 255)
  @test isequal(dat.triggers.count[255], 39)

  # specific channels
  bdf_file1 = joinpath(dirname(@__FILE__), "Newtest17-256.bdf")
  bdf_file2 = joinpath(dirname(@__FILE__), "Newtest17-256.bdf")
  dat1 = read_bdf(bdf_file1)
  dat2 = read_bdf(bdf_file1, channels=[1, 3, 5])

  @test isequal(dat1.data[:, 1], dat2.data[:, 1])
  @test isequal(dat1.data[:, 5], dat2.data[:, 3])
  @test isequal(dat2.header.num_bytes_header, 5 * 256)
  @test isequal(dat2.header.num_channels, 4)
  @test isequal(dat2.header.num_data_records, 60)
  @test isequal(dat2.header.sample_rate[1], 256)
  @test isequal(size(dat2.data), (15360, 3))
  @test isequal(dat2.triggers.idx[1], 415)
  @test isequal(dat2.triggers.val[1], 255)
  @test isequal(dat2.triggers.count[255], 40)

  bdf_file1 = joinpath(dirname(@__FILE__), "Newtest17-2048.bdf")
  bdf_file2 = joinpath(dirname(@__FILE__), "Newtest17-2048.bdf")
  dat1 = read_bdf(bdf_file1)
  dat2 = read_bdf(bdf_file2, channels=[1, 3, 5])

  @test isequal(dat1.data[:, 1], dat2.data[:, 1])
  @test isequal(dat1.data[:, 5], dat2.data[:, 3])
  @test isequal(dat2.header.num_bytes_header, 5 * 256)
  @test isequal(dat2.header.num_channels, 4)
  @test isequal(dat2.header.num_data_records, 60)
  @test isequal(dat2.header.sample_rate[1], 2048)
  @test isequal(size(dat2.data), (122880, 3))
  @test isequal(dat2.triggers.idx[1], 3353)
  @test isequal(dat2.triggers.val[1], 255)
  @test isequal(dat2.triggers.count[255], 39)

  bdf_file1 = joinpath(dirname(@__FILE__), "Newtest17-2048.bdf")
  bdf_file2 = joinpath(dirname(@__FILE__), "Newtest17-2048.bdf")
  dat1 = read_bdf(bdf_file1)
  dat2 = read_bdf(bdf_file2, channels=["A1", "A3", "A5"])

  @test isequal(dat1.data[:, 1], dat2.data[:, 1])
  @test isequal(dat1.data[:, 5], dat2.data[:, 3])
  @test isequal(dat2.header.num_bytes_header, 5 * 256)
  @test isequal(dat2.header.num_channels, 4)
  @test isequal(dat2.header.num_data_records, 60)
  @test isequal(dat2.header.sample_rate[1], 2048)
  @test isequal(size(dat2.data), (122880, 3))
  @test isequal(dat2.triggers.idx[1], 3353)
  @test isequal(dat2.triggers.val[1], 255)
  @test isequal(dat2.triggers.count[255], 39)

  bdf_file1 = joinpath(dirname(@__FILE__), "Newtest17-2048.bdf")
  bdf_file2 = joinpath(dirname(@__FILE__), "Newtest17-2048.bdf")
  dat1 = read_bdf(bdf_file1)
  dat2 = read_bdf(bdf_file2, channels=["A1", "A3", "A5"])

  @test isequal(dat1.data[:, 1], dat2.data[:, 1])
  @test isequal(dat1.data[:, 5], dat2.data[:, 3])
  @test isequal(dat2.header.num_bytes_header, 5 * 256)
  @test isequal(dat2.header.num_channels, 4)
  @test isequal(dat2.header.num_data_records, 60)
  @test isequal(dat2.header.sample_rate[1], 2048)
  @test isequal(size(dat2.data), (122880, 3))
  @test isequal(dat2.triggers.idx[1], 3353)
  @test isequal(dat2.triggers.val[1], 255)
  @test isequal(dat2.triggers.count[255], 39)

end


@testset "read_bdf_write_bdf" begin

  bdf_file1 = joinpath(dirname(@__FILE__), "Newtest17-256.bdf")
  dat1 = read_bdf(bdf_file1)
  dat1.filename = "Newtest17-256_write.bdf"
  write_bdf(dat1)
  bdf_file2 = joinpath(dirname(@__FILE__), "Newtest17-256_write.bdf")
  dat2 = read_bdf(bdf_file2)

  @test isequal(dat1.data, dat2.data)
  @test isequal(dat1.status, dat2.status)
  @test isequal(dat1.time, dat2.time)
  @test isequal(dat1.triggers.count, dat2.triggers.count)
  @test isequal(dat1.triggers.idx, dat2.triggers.idx)
  @test isequal(dat1.triggers.raw, dat2.triggers.raw)
  @test isequal(dat1.triggers.time, dat2.triggers.time)
  @test isequal(dat1.triggers.val, dat2.triggers.val)

  rm(bdf_file2)


  bdf_file1 = joinpath(dirname(@__FILE__), "Newtest17-2048.bdf")
  dat1 = read_bdf(bdf_file1)
  dat1.filename = "Newtest17-2048_write.bdf"
  write_bdf(dat1)

  bdf_file2 = joinpath(dirname(@__FILE__), "Newtest17-2048_write.bdf")
  dat2 = read_bdf(bdf_file2)

  @test isequal(dat1.data, dat2.data)
  @test isequal(dat1.status, dat2.status)
  @test isequal(dat1.time, dat2.time)
  @test isequal(dat1.triggers.count, dat2.triggers.count)
  @test isequal(dat1.triggers.idx, dat2.triggers.idx)
  @test isequal(dat1.triggers.raw, dat2.triggers.raw)
  @test isequal(dat1.triggers.time, dat2.triggers.time)
  @test isequal(dat1.triggers.val, dat2.triggers.val)

  rm(bdf_file2)

end


@testset "delete_channels_bdf" begin

  bdf_file = joinpath(dirname(@__FILE__), "Newtest17-256.bdf")
  dat = read_bdf(bdf_file)
  @test isequal(dat.header.num_channels, 17)
  dat = delete_channels_bdf(dat, [1])
  @test isequal(dat.header.num_channels, 16)

  bdf_file = joinpath(dirname(@__FILE__), "Newtest17-2048.bdf")
  dat = read_bdf(bdf_file)
  @test isequal(dat.header.num_channels, 17)
  dat = delete_channels_bdf(dat, [1])
  @test isequal(dat.header.num_channels, 16)

  bdf_file = joinpath(dirname(@__FILE__), "Newtest17-256.bdf")
  dat = read_bdf(bdf_file)
  @test isequal(dat.header.num_channels, 17)
  dat = delete_channels_bdf(dat, 1)
  @test isequal(dat.header.num_channels, 16)

  bdf_file = joinpath(dirname(@__FILE__), "Newtest17-2048.bdf")
  dat = read_bdf(bdf_file)
  @test isequal(dat.header.num_channels, 17)
  dat = delete_channels_bdf(dat, 1)
  @test isequal(dat.header.num_channels, 16)

end

@testset "merge_bdf" begin

  bdf_file = joinpath(dirname(@__FILE__), "Newtest17-256.bdf")
  dat1 = read_bdf(bdf_file)
  dat2 = read_bdf(bdf_file)
  dat3 = merge_bdf([dat1, dat2])

  @test isequal(dat3.header.num_bytes_header, 18 * 256)
  @test isequal(dat3.header.num_channels, 17)
  @test isequal(dat3.header.num_data_records, 120)
  @test isequal(dat3.header.sample_rate[1], 256)
  @test isequal(size(dat3.data), (15360 * 2, 16))
  @test isequal(dat3.triggers.idx[1], 415)
  @test isequal(dat3.triggers.val[1], 255)
  @test isequal(dat3.triggers.count[255], 80)

  bdf_file = joinpath(dirname(@__FILE__), "Newtest17-2048.bdf")
  dat1 = read_bdf(bdf_file)
  dat2 = read_bdf(bdf_file)
  dat3 = merge_bdf([dat1, dat2])

  @test isequal(dat3.header.num_bytes_header, 18 * 256)
  @test isequal(dat3.header.num_channels, 17)
  @test isequal(dat3.header.num_data_records, 120)
  @test isequal(dat3.header.sample_rate[1], 2048)
  @test isequal(size(dat3.data), (122880 * 2, 16))
  @test isequal(dat3.triggers.idx[1], 3353)
  @test isequal(dat3.triggers.val[1], 255)
  @test isequal(dat3.triggers.count[255], 79)


end


@testset "select_channels_bdf" begin

  bdf_file = joinpath(dirname(@__FILE__), "Newtest17-256.bdf")
  dat1 = read_bdf(bdf_file)
  dat2 = select_channels_bdf(dat1, [1, 3, 5])

  @test isequal(dat1.data[:, 1], dat2.data[:, 1])
  @test isequal(dat1.data[:, 5], dat2.data[:, 3])
  @test isequal(dat2.header.num_bytes_header, 5 * 256)
  @test isequal(dat2.header.num_channels, 4)
  @test isequal(dat2.header.num_data_records, 60)
  @test isequal(dat2.header.sample_rate[1], 256)
  @test isequal(size(dat2.data), (15360, 3))
  @test isequal(dat2.triggers.idx[1], 415)
  @test isequal(dat2.triggers.val[1], 255)
  @test isequal(dat2.triggers.count[255], 40)

  bdf_file = joinpath(dirname(@__FILE__), "Newtest17-2048.bdf")
  dat1 = read_bdf(bdf_file)
  dat2 = select_channels_bdf(dat1, ["A1", "A3", "A5"])

  @test isequal(dat1.data[:, 1], dat2.data[:, 1])
  @test isequal(dat1.data[:, 5], dat2.data[:, 3])
  @test isequal(dat2.header.num_bytes_header, 5 * 256)
  @test isequal(dat2.header.num_channels, 4)
  @test isequal(dat2.header.num_data_records, 60)
  @test isequal(dat2.header.sample_rate[1], 2048)
  @test isequal(size(dat2.data), (122880, 3))
  @test isequal(dat2.triggers.idx[1], 3353)
  @test isequal(dat2.triggers.val[1], 255)
  @test isequal(dat2.triggers.count[255], 39)



end



@testset "crop_bdf" begin

  bdf_file = joinpath(dirname(@__FILE__), "Newtest17-256.bdf")
  dat1 = read_bdf(bdf_file)
  dat2 = crop_bdf(dat1, "records", [10 20])

  @test !isequal(size(dat1.data, 1), size(dat2.data, 1))
  @test isequal(size(dat1.data, 2), size(dat2.data, 2))
  @test isequal(dat2.header.num_data_records, 11)

  bdf_file = joinpath(dirname(@__FILE__), "Newtest17-2048.bdf")
  dat1 = read_bdf(bdf_file)
  dat2 = crop_bdf(dat1, "records", [10 20])

  @test !isequal(size(dat1.data, 1), size(dat2.data, 1))
  @test isequal(size(dat1.data, 2), size(dat2.data, 2))
  @test isequal(dat2.header.num_data_records, 11)



end


@testset "downsample_bdf" begin

  bdf_file = joinpath(dirname(@__FILE__), "Newtest17-256.bdf")
  dat1 = read_bdf(bdf_file)
  dat2 = downsample_bdf(dat1, 2)

  @test isequal(dat2.header.sample_rate[1], 128)

  bdf_file = joinpath(dirname(@__FILE__), "Newtest17-2048.bdf")
  dat1 = read_bdf(bdf_file)
  dat2 = downsample_bdf(dat1, 2)

  @test isequal(dat2.header.sample_rate[1], 1024)



end

@testset "channel_idx" begin

  x = BioSemiBDF.channel_idx(["A1", "A2", "A3"], ["A1", "A3"])
  @test isequal(x, [1, 3])

  x = BioSemiBDF.channel_idx(["A1", "A2", "A3"], "A1")
  @test isequal(x, [1, 3])

  x = BioSemiBDF.channel_idx(["A1", "A2", "A3"], "A3")
  @test isequal(x, [3])

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
