# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# The code adapted from Tempo.jl reproduces Tempo.jl.
#
# Truth: Tempo.jl v1.3.1's own outputs, captured by calling it before it was removed from
# AstroEpochs' dependencies (see THIRD_PARTY_NOTICES.md). The offsets and calendar routines were
# brought in unchanged, so the comparison is exact: any difference is a transcription error, not a
# tolerance. Physical correctness is tested elsewhere, against Astropy, in
# test_scaleconversions.jl and test_formatconversions.jl; this file only says the move lost nothing.
#
# The points: offsets from 1900 to 2100 including the 1977 TCG/TCB epoch; dates across leap years,
# the century rules and month ends; day fractions at their edges; and Julian dates split between
# the two parts in each of the ways the calendar routine allows.

using Test
using AstroEpochs
using AstroEpochs: offset_tt2tai, offset_tai2tt, offset_tt2tdb, offset_tdb2tt, offset_tcg2tt,
                   offset_tt2tcg, offset_tcb2tdb, offset_tdb2tcb, cal2jd, jd2cal, hms2fd, fd2hms,
                   _builtin_leap_table

const _OFFSET_SECONDS = [-3.15576e9, -1.5e9, -7.25803167816e8, -1.0e6, -1.0, 0.0, 1.0, 1.0e6, 2.5e8, 8.5e8, 1.5e9, 3.15576e9]
const _OFFSETS = Dict(
    :offset_tt2tai => [-32.184, -32.184, -32.184, -32.184, -32.184, -32.184, -32.184, -32.184, -32.184, -32.184, -32.184, -32.184],
    :offset_tai2tt => [32.184, 32.184, 32.184, 32.184, 32.184, 32.184, 32.184, 32.184, 32.184, 32.184, 32.184, 32.184],
    :offset_tt2tdb => [-4.480752804288902e-5, 0.00038302803026703274, -5.906411547295813e-5, -0.0004039976869954267, -7.273711127905638e-5, -7.273677619130569e-5, -7.2736441103552e-5, 0.0002615398804270446, -0.0008561343776250286, -0.0007418862478001235, -0.00024484077380588793, -0.00010064500530621369],
    :offset_tdb2tt => [4.480752802801252e-5, -0.00038302803038799566, 5.90641154531909e-5, 0.0004039976868640724, 7.273711125468333e-5, 7.273677616693264e-5, 7.273644107918043e-5, -0.00026153988034043594, 0.0008561343773797322, 0.0007418862475780317, 0.0002448407738845113, 0.0001006450052725561],
    :offset_tcg2tt => [1.6935074176585845, 0.5395602344314004, -0.0, -0.5051363566551995, -0.5058332849716705, -0.5058332856685995, -0.5058332863655285, -0.5065302146819995, -0.6800655390185995, -1.0982229470585996, -1.5512268057685994, -2.7051739889957833],
    :offset_tt2tcg => [-1.6935074188388388, -0.5395602348074356, 0.0, 0.5051363570072437, 0.5058332853242004, 0.5058332860211293, 0.5058332867180584, 0.506530215035015, 0.6800655394925569, 1.0982229478239829, 1.5512268068496944, 2.7051739908810974],
    :offset_tcb2tdb => [37.67696103687951, 12.004074926242707, -0.0, -11.238216396077295, -11.253721578252097, -11.253721593757295, -11.253721609262492, -11.269226791437294, -15.130021013757295, -24.433139621757295, -34.5115181137573, -60.1844042243941],
    :offset_tdb2tcb => [-37.676961621068244, -12.004075112368264, 0.0, 11.238216570328063, 11.253721752743276, 11.253721768248475, 11.253721783753672, 11.269226966168885, 15.130021248351264, 24.43314000059796, 34.51151864886521, 60.184405157565195],
)
const _CAL2JD = [
    (1583, 1, 1) => (2451545, -152306),
    (1600, 2, 29) => (2451545, -146038),
    (1700, 2, 28) => (2451545, -109514),
    (1858, 11, 17) => (2451545, -51544),
    (1900, 3, 1) => (2451545, -36465),
    (1972, 1, 1) => (2451545, -10227),
    (1999, 12, 31) => (2451545, -1),
    (2000, 1, 1) => (2451545, 0),
    (2000, 2, 29) => (2451545, 59),
    (2016, 12, 31) => (2451545, 6209),
    (2017, 1, 1) => (2451545, 6210),
    (2024, 6, 1) => (2451545, 8918),
    (2100, 2, 28) => (2451545, 36583),
    (2400, 12, 31) => (2451545, 146462),
]
const _HMS2FD = [
    (0, 0, 0.0) => 0.0,
    (0, 0, 59.999) => 0.0006944328703703704,
    (6, 30, 15.25) => 0.27100983796296296,
    (12, 0, 0.0) => 0.5,
    (23, 59, 59.999999) => 0.999999999988426,
]
const _FD2HMS = [
    0.0 => (0, 0, 0.0),
    1.0e-9 => (0, 0, 8.64e-5),
    0.25 => (6, 0, 0.0),
    0.5 => (12, 0, 0.0),
    0.7071067811865476 => (16, 58, 14.025894517712004),
    0.999999999 => (23, 59, 59.999913599996944),
]
const _JD2CAL = [
    (2.451545e6, 0.0) => (2000, 1, 1, 0.5),
    (2.451545e6, 365.5) => (2001, 1, 1, 0.0),
    (2.45191e6, 0.5) => (2001, 1, 1, 0.0),
    (2.4000005e6, 50123.2) => (1996, 2, 10, 0.19999999999708962),
    (2.4501235e6, 0.2) => (1996, 2, 10, 0.2),
    (2.457754e6, 0.4999999999) => (2016, 12, 31, 0.9999999999),
    (2.459143e6, -0.25) => (2020, 10, 20, 0.25),
    (2.4604625e6, 0.499999999) => (2024, 6, 1, 0.499999999),
    (2.2991605e6, 0.0) => (1582, 10, 15, 0.0),
    (2.4880695e6, 0.75) => (2100, 1, 1, 0.75),
]
const _TEMPO_LEAP_STARTS = [-10227.5, -10045.5, -9861.5, -9496.5, -9131.5, -8766.5, -8400.5, -8035.5, -7670.5, -7305.5, -6758.5, -6393.5, -6028.5, -5297.5, -4383.5, -3652.5, -3287.5, -2740.5, -2375.5, -2010.5, -1461.5, -914.5, -365.5, 2191.5, 3287.5, 4564.5, 5659.5, 6209.5]
const _TEMPO_LEAP_DELTA  = [10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0]

@testset "Tempo parity — scale offsets" begin
    for (name, expected) in _OFFSETS
        f = getfield(AstroEpochs, name)
        @test [f(s) for s in _OFFSET_SECONDS] == expected
    end
end

@testset "Tempo parity — calendar" begin
    for ((y, m, d), expected) in _CAL2JD
        @test cal2jd(y, m, d) == expected
    end
    for ((h, m, s), expected) in _HMS2FD
        @test hms2fd(h, m, s) == expected
    end
    for (fd, expected) in _FD2HMS
        @test fd2hms(fd) == expected
    end
    for ((a, b), expected) in _JD2CAL
        @test jd2cal(a, b) == expected
    end
end

@testset "Tempo parity — the built-in leap-second table" begin
    table = _builtin_leap_table()
    @test table.starts == _TEMPO_LEAP_STARTS
    @test table.delta  == _TEMPO_LEAP_DELTA
end

@testset "Tempo parity — the domain errors are kept" begin
    @test_throws DomainError cal2jd(1582, 12, 31)
    @test_throws DomainError cal2jd(2023, 2, 29)
    @test_throws DomainError cal2jd(2024, 13, 1)
    @test_throws DomainError hms2fd(24, 0, 0.0)
    @test_throws DomainError hms2fd(0, 0, 60.0)
    @test_throws DomainError fd2hms(1.5)
    @test_throws DomainError jd2cal(-1.0e5, 0.0)
end
