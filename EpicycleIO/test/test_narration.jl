# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# The caption strip: what `narrate` publishes is what the dashboard reads, so it is checked at the
# file rather than in a browser.

using Test
using EpicycleIO
using EpicycleIO: narrate, clear_narration   # public, not exported
using JSON

const _IO = EpicycleIO
_IO.auto_open!(false)

_caption() = (_IO.flush_publishes(); JSON.parse(read(joinpath(_IO.assets_dir(), "narration.json"), String)))

@testset "narration" begin

@testset "a caption carries every field the strip draws" begin
    narrate("Propagating a quarter of a day.";
            title = "A spacecraft in orbit", subtitle = "Act 1",
            code = "propagate!(prop, sat, StopAt(sat, PropDurationDays(), 0.25))",
            facts = ["inclination" => "51.6 deg", "period" => 92.6],
            status = "running")
    c = _caption()

    @test c["title"]    == "A spacecraft in orbit"
    @test c["subtitle"] == "Act 1"
    @test c["text"]     == "Propagating a quarter of a day."
    @test occursin("propagate!", c["code"])
    @test c["status"]   == "running"
    # Values are shown as given, so a number arrives as its string.
    @test c["facts"] == [["inclination", "51.6 deg"], ["period", "92.6"]]
end

@testset "each caption supersedes the last, which is how the page knows it changed" begin
    narrate("first")
    s1 = _caption()["seq"]
    narrate("second")
    c = _caption()
    @test c["seq"] > s1
    @test c["text"] == "second"
end

@testset "clearing leaves nothing for the strip to draw" begin
    narrate("something"; title = "Title")
    clear_narration()
    c = _caption()
    @test !haskey(c, "title") && !haskey(c, "text")
end

@testset "a caption for some panels stays with them after the script moves on" begin
    # This is what lets a tab opened for one act keep that act's words: the Act 3 tab must not
    # show Act 6's caption when the presenter goes back to it.
    narrate("coasting to apoapsis"; title = "Act 3", panels = ["Act 3 · Transfer", "Act 3 · Radius"])
    narrate("touching down";        title = "Act 6", panels = ["Act 6 · Descent"])
    _IO.flush_publishes()

    read_caption(name) = JSON.parse(read(joinpath(_IO.assets_dir(), name), String))
    act3 = read_caption(_IO.narration_file(["Act 3 · Radius", "Act 3 · Transfer"]))   # order does not matter
    @test act3["title"] == "Act 3"
    @test act3["text"]  == "coasting to apoapsis"

    # The unscoped dashboard follows the latest caption.
    @test _caption()["title"] == "Act 6"

    # The file name is the one the dashboard derives from its scope: panel ids, sorted, joined.
    @test _IO.narration_file(["Act 3 · Transfer", "Act 3 · Radius"]) ==
          "narration_Act_3___Radius~Act_3___Transfer.json"
    @test _IO.narration_file(String[]) == "narration.json"

    # Clearing takes the kept captions away too.
    clear_narration()
    @test !isfile(joinpath(_IO.assets_dir(), _IO.narration_file(["Act 3 · Transfer", "Act 3 · Radius"])))
end

@testset "facts must be label => value pairs, and say so" begin
    err = try
        narrate("x"; facts = ["not a pair"])
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("\"label\" => \"value\"", err.msg)
end

end # testset
