# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Where the leap-second table comes from when the network is not simply there: the stored list,
# a fresh download, an expired copy, or the built-in table, and what each says. Every test works
# in a temporary folder with URLs that fail on purpose, so none of it touches the network or the
# stored list the session uses.
#
# The URLs:
#   refused     a port nothing listens on, so the connection is refused at once
#   file://     a list written here, standing in for a successful download
#   silent      a local server that accepts the connection and never answers, which without a
#               timeout would block the first UTC conversion indefinitely

using Test
using AstroEpochs
using AstroEpochs: leap_second_table, _fetch_leap_second_table, _refresh_leap_seconds!,
                   _try_parse
using Sockets

const _NTP_JD = 2415020.5                                   # 1900-01-01, the list's time origin
_ntp(jd) = round(Int, (jd - _NTP_JD) * 86400)
_today_jd() = time() / 86400 + 2440587.5

# A two-entry list in the IANA format, expiring at `expires_jd`.
_list(expires_jd) = """
    #\tA test list
    #@\t$(_ntp(expires_jd))
    2272060800\t10\t# 1 Jan 1972
    3692217600\t37\t# 1 Jan 2017
    """

_file_url(path) = "file:///" * replace(abspath(path), "\\" => "/")
const _REFUSED = "http://127.0.0.1:1/leap-seconds.list"
_expires_jd(table) = table.expires + 2451545.0

@testset "leap-second list: where the table comes from" begin
    mktempdir() do dir
        path = joinpath(dir, "leap-seconds.list")
        now = _today_jd()

        @testset "a current stored list is used, and nothing is downloaded" begin
            write(path, _list(now + 100))
            table = @test_logs _fetch_leap_second_table(path; url = _REFUSED)
            @test _expires_jd(table) ≈ now + 100 atol = 1e-5
        end

        @testset "an expired list and no network: the expired copy, with a warning" begin
            write(path, _list(now - 100))
            table = @test_logs (:warn, r"expired copy") _fetch_leap_second_table(path; url = _REFUSED)
            @test _expires_jd(table) ≈ now - 100 atol = 1e-5
            @test table.delta[end] == 37.0
        end

        @testset "no list and no network: the built-in table, with a warning" begin
            rm(path)
            table = @test_logs (:warn, r"built-in table") _fetch_leap_second_table(path; url = _REFUSED)
            @test table.expires == -Inf
            @test table.delta[end] == 37.0
            @test !isfile(path)
        end

        @testset "an expired list and a download: the fresh list, stored for next time" begin
            write(path, _list(now - 100))
            fresh = joinpath(dir, "fresh.list")
            write(fresh, _list(now + 200))
            table = @test_logs _fetch_leap_second_table(path; url = _file_url(fresh))
            @test _expires_jd(table) ≈ now + 200 atol = 1e-5
            @test _try_parse(read(path, String)).expires == table.expires
            @test !isfile(path * ".download")
        end

        @testset "a download that is not a list: the stored copy is kept, nothing left behind" begin
            write(path, _list(now - 100))
            stored = read(path)
            bad = joinpath(dir, "bad.list")
            write(bad, "<html>not a leap-second list</html>\n")
            table = @test_logs (:warn, r"expired copy") _fetch_leap_second_table(path; url = _file_url(bad))
            @test _expires_jd(table) ≈ now - 100 atol = 1e-5
            @test read(path) == stored
            @test !isfile(path * ".download")
        end

        @testset "a server that never answers: the timeout ends the wait" begin
            port, server = listenany(ip"127.0.0.1", 49152)
            held = TCPSocket[]
            acceptor = @async try
                while true
                    push!(held, accept(server))             # accept, and never reply
                end
            catch
            end
            try
                write(path, _list(now - 100))
                silent = "http://127.0.0.1:$(port)/leap-seconds.list"
                elapsed = @elapsed begin
                    table = @test_logs (:warn, r"expired copy") _fetch_leap_second_table(
                        path; url = silent, timeout = 1.0)
                end
                @test _expires_jd(table) ≈ now - 100 atol = 1e-5
                @test elapsed < 10.0
                @test !isfile(path * ".download")
            finally
                close(server)
                foreach(close, held)
            end
        end

        @testset "a refresh that fails throws, and the table in use is unchanged" begin
            in_use = leap_second_table()
            @test_throws Exception _refresh_leap_seconds!(path, _REFUSED, 1.0)
            @test leap_second_table() === in_use
        end

        @testset "a refresh that succeeds is used from then on" begin
            in_use = leap_second_table()
            fresh = joinpath(dir, "refreshed.list")
            write(fresh, _list(now + 300))
            try
                table = _refresh_leap_seconds!(path, _file_url(fresh), 1.0)
                @test leap_second_table() === table
                @test _expires_jd(table) ≈ now + 300 atol = 1e-5
                @test _try_parse(read(path, String)).expires == table.expires
            finally
                AstroEpochs._LEAP[] = in_use              # the session's own table back
            end
            @test leap_second_table() === in_use
        end
    end
end
