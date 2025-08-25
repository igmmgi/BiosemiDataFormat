using BiosemiDataFormat
using Test

testfile(name) = joinpath(@__DIR__, name)

const cases = [
  ("Newtest17-256.bdf", 256, 60, (idx1=415, val1=255, cnt255=40), (rows=15360, chans=16)),
  ("Newtest17-2048.bdf", 2048, 60, (idx1=3353, val1=255, cnt255=39), (rows=122880, chans=16)),
]

@testset "BiosemiDataFormat" begin

  @testset "read_bdf (header + data + selections)" begin
    for (fname, sr, nrecs, trig, sz) in cases
      @testset "read_bdf $(sr)Hz" begin
        bdf = testfile(fname)

        # header-only
        hdr = read_bdf(bdf, header_only=true)
        @test isa(hdr, BiosemiDataFormat.BiosemiHeader)
        @test hdr.num_bytes_header == 18 * 256
        @test hdr.num_channels == 17
        @test hdr.num_data_records == nrecs
        @test hdr.sample_rate[1] == sr

        # full file
        dat = read_bdf(bdf)
        @test dat.header.num_bytes_header == 18 * 256
        @test dat.header.num_channels == 17
        @test dat.header.num_data_records == nrecs
        @test dat.header.sample_rate[1] == sr
        @test size(dat.data) == (sz.rows, sz.chans)
        @test dat.triggers.idx[1] == trig.idx1
        @test dat.triggers.val[1] == trig.val1
        @test dat.triggers.count[trig.val1] == trig.cnt255

        # selection by indices
        dat2 = read_bdf(bdf, channels=[1, 3, 5])
        @test dat.data[:, 1] == dat2.data[:, 1]
        @test dat.data[:, 5] == dat2.data[:, 3]
        @test dat2.header.num_bytes_header == 5 * 256
        @test dat2.header.num_channels == 4
        @test dat2.header.num_data_records == nrecs
        @test dat2.header.sample_rate[1] == sr
        @test size(dat2.data) == (sz.rows, 3)
        @test dat2.triggers.idx[1] == trig.idx1
        @test dat2.triggers.val[1] == trig.val1
        @test dat2.triggers.count[trig.val1] == trig.cnt255

        # selection by labels (only assert known labels on 2048Hz fixture)
        if sr == 2048
          dat3 = read_bdf(bdf, channels=["A1", "A3", "A5"])
          @test dat.data[:, 1] == dat3.data[:, 1]
          @test dat.data[:, 5] == dat3.data[:, 3]
          @test dat3.header.num_bytes_header == 5 * 256
          @test dat3.header.num_channels == 4
          @test dat3.header.num_data_records == nrecs
          @test dat3.header.sample_rate[1] == sr
          @test size(dat3.data) == (sz.rows, 3)
          @test dat3.triggers.idx[1] == trig.idx1
          @test dat3.triggers.val[1] == trig.val1
          @test dat3.triggers.count[trig.val1] == trig.cnt255
        end
      end
    end
  end

  @testset "write_bdf roundtrip" begin
    for (fname, sr, nrecs, trig, sz) in cases
      mktempdir() do tmp
        src = testfile(fname)
        dat1 = read_bdf(src)
        out = joinpath(tmp, "out.bdf")
        write_bdf(dat1, out)
        dat2 = read_bdf(out)

        @test dat1.data == dat2.data
        @test dat1.status == dat2.status
        @test dat1.time == dat2.time
        @test dat1.triggers.count == dat2.triggers.count
        @test dat1.triggers.idx == dat2.triggers.idx
        @test dat1.triggers.raw == dat2.triggers.raw
        @test dat1.triggers.time == dat2.triggers.time
        @test dat1.triggers.val == dat2.triggers.val
      end
    end
  end

  @testset "delete_channels_bdf" begin
    for (fname, sr, nrecs, trig, sz) in cases
      bdf = testfile(fname)
      dat = read_bdf(bdf)
      @test dat.header.num_channels == 17
      dat = delete_channels_bdf(dat, [1])
      @test dat.header.num_channels == 16
    end
  end

  @testset "merge_bdf extras" begin
    for (fname, sr, nrecs, trig, sz) in cases
      bdf = testfile(fname)
      d1 = read_bdf(bdf)
      d2 = read_bdf(bdf)
      m = merge_bdf([d1, d2])

      @test m.header.num_bytes_header == 18 * 256
      @test m.header.num_channels == 17
      @test m.header.num_data_records == nrecs * 2
      @test m.header.sample_rate[1] == sr
      @test size(m.data) == (sz.rows * 2, sz.chans)
      @test m.triggers.idx[1] == trig.idx1
      @test m.triggers.val[1] == trig.val1
      # trigger count for primary value: known behavior differs per fixture
      if sr == 256
        @test m.triggers.count[trig.val1] == trig.cnt255 * 2
      else
        @test m.triggers.count[trig.val1] == 79
      end
    end
  end

  @testset "channel label/metadata integrity" begin
    bdf = testfile("Newtest17-256.bdf")
    dat = read_bdf(bdf)
    sel = select_channels_bdf(dat, [1, 3, 5])
    @test sel.header.num_channels == 4 # includes status
    @test length(sel.header.channel_labels) == sel.header.num_channels
    @test sel.header.num_bytes_header == (sel.header.num_channels + 1) * 256
  end

  @testset "crop_bdf" begin
    for (fname, sr, nrecs, trig, sz) in cases
      bdf = testfile(fname)
      dat = read_bdf(bdf)
      cropped = crop_bdf(dat, "records", [10 20])
      @test size(dat.data, 1) != size(cropped.data, 1)
      @test size(dat.data, 2) == size(cropped.data, 2)
      @test cropped.header.num_data_records == 11

      # trigger-based crop between first and last distinct trigger
      uniq = sort!(collect(keys(dat.triggers.count)))
      if length(uniq) >= 2
        cropped2 = crop_bdf(dat, "triggers", [uniq[1] uniq[end]])
        @test size(cropped2.data, 2) == size(dat.data, 2)
        @test cropped2.header.num_data_records > 0
        @test first(cropped2.time) == 0
      end
    end
  end

  @testset "downsample_bdf details" begin
    for (fname, sr, nrecs, trig, sz) in cases
      dat = read_bdf(testfile(fname))
      dec = 2
      ds = downsample_bdf(dat, dec)
      @test ds.header.sample_rate[1] == div(dat.header.sample_rate[1], dec)
      @test ds.header.num_samples[1] == div(dat.header.num_samples[1], dec)
      @test size(ds.data, 1) == div(size(dat.data, 1), dec)
      @test length(ds.time) == size(ds.data, 1)
      @test length(ds.triggers.raw) == size(ds.data, 1)
      @test ds.triggers.idx[1] == round(Int, dat.triggers.idx[1] / dec)
    end
  end

  @testset "read triggers-only" begin
    for (fname, sr, nrecs, trig, sz) in cases
      dat = read_bdf(testfile(fname), channels=[-1])
      @test size(dat.data, 2) == 0
      @test length(dat.triggers.raw) == length(dat.time)
    end
  end

  @testset "errors" begin
    @test_throws ErrorException read_bdf(testfile("does_not_exist.bdf"))
    dat = read_bdf(testfile("Newtest17-256.bdf"))
    @test_throws ErrorException downsample_bdf(dat, 3)
    @test_throws ErrorException crop_bdf(dat, "invalid", [1 2])
    @test_throws ErrorException crop_bdf(dat, "records", [1])
    @test_throws ErrorException crop_bdf(dat, "records", [0 20])
    @test_throws ErrorException crop_bdf(dat, "records", [10 61])
    @test_throws ErrorException crop_bdf(dat, "triggers", [999 1000])
  end

  @testset "merge_bdf error paths" begin
    bdf = testfile("Newtest17-256.bdf")
    d1 = read_bdf(bdf)
    d2 = read_bdf(bdf)
    d3 = read_bdf(bdf)

    # Force header mismatches
    d2.header.num_channels = 16
    @test_throws ErrorException merge_bdf([d1, d2])

    d2.header.num_channels = 17
    d2.header.channel_labels = d2.header.channel_labels[1:end-1]
    @test_throws ErrorException merge_bdf([d1, d2])

    d2.header.channel_labels = d1.header.channel_labels
    d2.header.sample_rate = [512]
    @test_throws ErrorException merge_bdf([d1, d2])
  end

  @testset "bang variants mutate in place" begin
    for (fname, sr, nrecs, trig, sz) in cases
      bdf = testfile(fname)
      dat = read_bdf(bdf)
      orig_channels = dat.header.num_channels
      orig_bytes = dat.header.num_bytes_header
      orig_labels = copy(dat.header.channel_labels)

      # Test select_channels_bdf!
      select_channels_bdf!(dat, [1, 3, 5])
      @test dat.header.num_channels == 4
      @test dat.header.num_bytes_header == 5 * 256
      @test length(dat.header.channel_labels) == 4
      @test size(dat.data, 2) == 3

      # Test delete_channels_bdf!
      delete_channels_bdf!(dat, [1])
      @test dat.header.num_channels == 3
      @test dat.header.num_bytes_header == 4 * 256
      @test length(dat.header.channel_labels) == 3
      @test size(dat.data, 2) == 2
    end
  end

  @testset "write_bdf default filename" begin
    for (fname, sr, nrecs, trig, sz) in cases
      mktempdir() do tmp
        src = testfile(fname)
        dat = read_bdf(src)
        out = joinpath(tmp, "out.bdf")
        dat.filename = out
        write_bdf(dat)  # No filename arg
        @test isfile(out)
        dat2 = read_bdf(out)
        @test dat.data == dat2.data
      end
    end
  end

  @testset "Base.show methods" begin
    bdf = testfile("Newtest17-256.bdf")
    hdr = read_bdf(bdf, header_only=true)
    dat = read_bdf(bdf)
    trig = dat.triggers

    hdr_str = sprint(show, hdr)
    @test occursin("Number of Channels: 16", hdr_str)
    @test occursin("Sample Rate: 256", hdr_str)

    trig_str = sprint(show, trig)
    @test occursin("Triggers (Value => Count):", trig_str)

    dat_str = sprint(show, dat)
    @test occursin("Filename:", dat_str)
    @test occursin("Data Size:", dat_str)
  end

  @testset "time boundaries and consistency" begin
    for (fname, sr, nrecs, trig, sz) in cases
      dat = read_bdf(testfile(fname))
      @test length(dat.time) == size(dat.data, 1)
      @test last(dat.time) == (nrecs - 1/sr)
      @test first(dat.time) == 0

      # After cropping
      cropped = crop_bdf(dat, "records", [10 20])
      @test length(cropped.time) == size(cropped.data, 1)
      @test cropped.header.num_data_records == 11
    end
  end

  @testset "mixed channel selection including status" begin
    for (fname, sr, nrecs, trig, sz) in cases
      dat = read_bdf(testfile(fname), channels=[1, -1])
      @test dat.header.num_channels == 2
      @test size(dat.data, 2) == 1
      @test length(dat.triggers.raw) == size(dat.data, 1)
      @test dat.header.num_bytes_header == 3 * 256
    end
  end

  @testset "downsample additional cases" begin
    for (fname, sr, nrecs, trig, sz) in cases
      dat = read_bdf(testfile(fname))
      dec = 4
      ds = downsample_bdf(dat, dec)
      @test ds.header.sample_rate[1] == div(dat.header.sample_rate[1], dec)
      @test size(ds.data, 1) == div(size(dat.data, 1), dec)

      # Trigger indices within bounds
      valid_indices = filter(x -> 1 <= x <= size(ds.data, 1), ds.triggers.idx)
      @test count(x -> x != 0, ds.triggers.raw) == length(valid_indices)
      @test all(1 .<= ds.triggers.idx .<= size(ds.data, 1))
    end
  end

  @testset "header parsing robustness" begin
    bdf = testfile("Newtest17-256.bdf")
    dat = read_bdf(bdf)
    sel = select_channels_bdf(dat, [1, 3, 5])

    # Header field lengths track channel selections
    @test length(sel.header.channel_labels) == sel.header.num_channels
    @test length(sel.header.transducer_type) == sel.header.num_channels
    @test length(sel.header.channel_unit) == sel.header.num_channels
    @test length(sel.header.physical_min) == sel.header.num_channels
    @test length(sel.header.physical_max) == sel.header.num_channels
    @test length(sel.header.digital_min) == sel.header.num_channels
    @test length(sel.header.digital_max) == sel.header.num_channels
    @test length(sel.header.pre_filter) == sel.header.num_channels
    @test length(sel.header.num_samples) == sel.header.num_channels
    @test length(sel.header.reserved) == sel.header.num_channels
    @test length(sel.header.sample_rate) == sel.header.num_channels
    @test length(sel.header.scale_factor) == sel.header.num_channels
  end

  @testset "channel_index" begin
    x = BiosemiDataFormat.channel_index(["A1", "A2", "A3"], ["A1", "A3"])
    @test x == [1, 3]

    x = BiosemiDataFormat.channel_index(["A1", "A2", "A3"], "A1")
    @test x == [1, 3]

    x = BiosemiDataFormat.channel_index(["A1", "A2", "A3"], "A3")
    @test x == [3]

    x = BiosemiDataFormat.channel_index(["A1", "A2", "A3", "A4"], ["A1", "A4"])
    @test x == [1, 4]

    x = BiosemiDataFormat.channel_index(["A1", "A2", "A3", "A4"], ["A1", "A4", "A4"])
    @test x == [1, 4]

    @test_throws ErrorException BiosemiDataFormat.channel_index(["A1"], ["zzz"])

    x = BiosemiDataFormat.channel_index(["A1", "A2", "A3"], [1, 3])
    @test x == [1, 3]

    x = BiosemiDataFormat.channel_index(["A1", "A2", "A3", "A4"], [1, 4])
    @test x == [1, 4]

    x = BiosemiDataFormat.channel_index(["A1", "A2", "A3", "A4"], [1, 4, 4])
    @test x == [1, 4]

    @test_throws ErrorException BiosemiDataFormat.channel_index(["A1"], [2])
    @test_throws ErrorException BiosemiDataFormat.channel_index(["A1"], [1, 2])
  end

end
