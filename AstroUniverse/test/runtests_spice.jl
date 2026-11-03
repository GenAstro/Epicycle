# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

using Test
using AstroUniverse
using SPICE
using Logging
using Scratch: with_scratch_directory

# Managing kernels needs no ephemeris and no network, and these tests leave the real kernel cache
# alone. They run in a temporary cache holding a copy of the default kernels, and their extra kernel
# is a small file written here and fetched from a file URL.
_EXTRA = "extra_constants.tpc"
AstroUniverse.require_spice_kernels()
_DEFAULT_PATHS = [joinpath(get_spice_directory(), k) for k in AstroUniverse.DEFAULT_KERNELS]

mktempdir() do dir
_source = joinpath(dir, _EXTRA)
write(_source, "KPL/PCK\n\n\\begindata\n   EPICYCLE_TEST_VALUE = ( 1.0 )\n\\begintext\n")
_url = "file:///" * lstrip(replace(_source, '\\' => '/'), '/')
_scratch = mkpath(joinpath(dir, "scratch"))

with_scratch_directory(_scratch) do
unload_all_spice_kernels()
for path in _DEFAULT_PATHS
    cp(path, joinpath(get_spice_directory(), basename(path)))
end
for k in AstroUniverse.DEFAULT_KERNELS
    load_spice_kernel(k)
end

@testset "SPICE Kernel Management" begin
    
    @testset "Default Kernels on Init" begin
        # Verify default kernels are loaded after module initialization
        @test SPICE.ktotal("ALL") >= 2
        
        # Check that storage directory exists
        cache_dir = get_spice_directory()
        @test isdir(cache_dir)
        
        # Verify default kernel files exist on disk
        @test isfile(joinpath(cache_dir, "naif0012.tls"))
        @test isfile(joinpath(cache_dir, "epicycle_de440_1950-2100.bsp"))
        
        # Verify kernels are actually cached (loaded in SPICE memory)
        loaded_kernels = String[]
        for i in 1:SPICE.ktotal("ALL")
            result = SPICE.kdata(i, "ALL")
            if result !== nothing
                push!(loaded_kernels, basename(result[1]))
            end
        end
        for filename in AstroUniverse.DEFAULT_KERNELS
            @test filename in loaded_kernels
        end
    end
    
    @testset "Download Operations" begin
        download_spice_kernel(_EXTRA, _url)
        
        cache_dir = get_spice_directory()
        kernel_path = joinpath(cache_dir, _EXTRA)
        
        # Verify file was downloaded
        @test isfile(kernel_path)
        
        @test read(kernel_path) == read(_source)

        # A file already in the cache is kept, so this does not go back to the source.
        rm(_source)
        @test_nowarn download_spice_kernel(_EXTRA, _url)
    end
    
    @testset "get_spice_directory()" begin
        dir = get_spice_directory()
        @test isdir(dir)
        @test contains(dir, "spice_kernels")
    end
    
    @testset "cached_spice_kernels()" begin
        @test "naif0012.tls" in cached_spice_kernels()

        # Verify default kernels are in the storage directory
        storage_dir = get_spice_directory()
        downloaded_files = readdir(storage_dir)
        @test "naif0012.tls" in downloaded_files
        @test "epicycle_de440_1950-2100.bsp" in downloaded_files
        @test _EXTRA in downloaded_files  # From earlier test
    end
    
    @testset "loaded_spice_kernels()" begin
        # Ensure the extra kernel is NOT loaded to start with clean state
        try
            unload_spice_kernel(_EXTRA)
        catch
            # Ignore if not already loaded
        end
        
        # After loading a kernel, verify it appears in the kernel pool
        download_spice_kernel(_EXTRA, _url)
        initial_count = SPICE.ktotal("ALL")
        load_spice_kernel(_EXTRA)
        @test SPICE.ktotal("ALL") == initial_count + 1
        
        # Verify the extra kernel is in the cached (loaded) list
        cached_kernels = String[]
        for i in 1:SPICE.ktotal("ALL")
            result = SPICE.kdata(i, "ALL")
            if result !== nothing
                push!(cached_kernels, basename(result[1]))
            end
        end
        @test _EXTRA in cached_kernels
        for filename in AstroUniverse.DEFAULT_KERNELS
            @test filename in cached_kernels
        end
        
        # Clean up
        unload_spice_kernel(_EXTRA)
    end
    
    @testset "Load/Unload Operations" begin
        # Ensure the extra kernel is NOT loaded to start with clean state
        try
            unload_spice_kernel(_EXTRA)
        catch
            # Ignore if not already loaded
        end
        
        # Count kernels before operations
        initial_count = SPICE.ktotal("ALL")
        
        # Load the extra kernel
        @test_nowarn load_spice_kernel(_EXTRA)
        @test SPICE.ktotal("ALL") == initial_count + 1
        
        # Verify the extra kernel is actually in the kernel pool
        loaded_kernels = String[]
        for i in 1:SPICE.ktotal("ALL")
            result = SPICE.kdata(i, "ALL")
            if result !== nothing
                push!(loaded_kernels, basename(result[1]))
            end
        end
        @test _EXTRA in loaded_kernels
        
        # Try to load again - should warn and not increase count
        @test_logs (:warn, r"already loaded") load_spice_kernel(_EXTRA)
        @test SPICE.ktotal("ALL") == initial_count + 1
        
        # Unload the extra kernel
        @test_nowarn unload_spice_kernel(_EXTRA)
        @test SPICE.ktotal("ALL") == initial_count
        
        # Verify the extra kernel is no longer in the kernel pool
        loaded_kernels = String[]
        for i in 1:SPICE.ktotal("ALL")
            result = SPICE.kdata(i, "ALL")
            if result !== nothing
                push!(loaded_kernels, basename(result[1]))
            end
        end
        @test !(_EXTRA in loaded_kernels)
        
        # Try to unload again - should warn
        @test_logs (:warn, r"not currently loaded") unload_spice_kernel(_EXTRA)
        @test SPICE.ktotal("ALL") == initial_count
        
        # Load it again for later tests
        load_spice_kernel(_EXTRA)
    end
    
    @testset "unload_all_spice_kernels()" begin
        # Ensure some kernels are loaded
        @test SPICE.ktotal("ALL") > 0
        
        # Clear all kernels
        @test_nowarn unload_all_spice_kernels()
        @test SPICE.ktotal("ALL") == 0
        
        # Put back exactly what AstroUniverse.__init__ loaded, read from DEFAULT_KERNELS so
        # this cannot drift from it again. It previously reloaded a hand-written pair —
        # naif0012.tls and de440.bsp — which is two of the four and the wrong ephemeris:
        # de440.bsp is not epicycle_de440_1950-2100.bsp and covers a different span. Every
        # suite running after this one in the same process then failed on
        # "Insufficient ephemeris data", and the `== 2` below asserted the broken state was
        # correct.
        for k in AstroUniverse.DEFAULT_KERNELS
            load_spice_kernel(k)
        end
        @test SPICE.ktotal("ALL") == length(AstroUniverse.DEFAULT_KERNELS)
    end
    
    @testset "Error Handling" begin
        # Loading non-existent kernel should error
        @test_throws ErrorException load_spice_kernel("nonexistent.bsp")
        
        # Unloading non-existent file (not in cache directory) should error
        @test_throws ErrorException unload_spice_kernel("nonexistent.bsp")
        
        # Error messages should be informative
        try
            load_spice_kernel("missing.bsp")
            @test false  # Should not reach here
        catch e
            @test occursin("not found in cache", e.msg)
            @test occursin("download_spice_kernel", e.msg)
        end
        
        try
            unload_spice_kernel("missing.bsp")
            @test false  # Should not reach here
        catch e
            @test occursin("not found in cache", e.msg)
            @test occursin("Cannot unload", e.msg)
        end
    end
    
    @testset "Kernel Swapping Workflow" begin
        # Ensure clean state - unload the extra kernel from previous test, ensure epicycle_de440_1950-2100.bsp is loaded
        try
            unload_spice_kernel(_EXTRA)
        catch
        end
        
        # Make sure epicycle_de440_1950-2100.bsp is loaded
        try
            load_spice_kernel("epicycle_de440_1950-2100.bsp")
        catch
        end
        initial_count = SPICE.ktotal("ALL")
        
        # Unload epicycle_de440_1950-2100.bsp
        unload_spice_kernel("epicycle_de440_1950-2100.bsp")
        @test SPICE.ktotal("ALL") == initial_count - 1
        
        # Verify epicycle_de440_1950-2100.bsp is no longer in the pool
        loaded_kernels = String[]
        for i in 1:SPICE.ktotal("ALL")
            result = SPICE.kdata(i, "ALL")
            if result !== nothing
                push!(loaded_kernels, basename(result[1]))
            end
        end
        @test !("epicycle_de440_1950-2100.bsp" in loaded_kernels)
        
        # Load the extra kernel instead
        load_spice_kernel(_EXTRA)
        @test SPICE.ktotal("ALL") == initial_count
        
        # Verify the extra kernel is in the pool
        loaded_kernels = String[]
        for i in 1:SPICE.ktotal("ALL")
            result = SPICE.kdata(i, "ALL")
            if result !== nothing
                push!(loaded_kernels, basename(result[1]))
            end
        end
        @test _EXTRA in loaded_kernels
        @test !("epicycle_de440_1950-2100.bsp" in loaded_kernels)
        
        # Swap back
        unload_spice_kernel(_EXTRA)
        load_spice_kernel("epicycle_de440_1950-2100.bsp")
        @test SPICE.ktotal("ALL") == initial_count
        
        # Verify final state: epicycle_de440_1950-2100.bsp is in pool, the extra kernel is not
        loaded_kernels = String[]
        for i in 1:SPICE.ktotal("ALL")
            result = SPICE.kdata(i, "ALL")
            if result !== nothing
                push!(loaded_kernels, basename(result[1]))
            end
        end
        @test "epicycle_de440_1950-2100.bsp" in loaded_kernels
        @test !(_EXTRA in loaded_kernels)
    end
    
    @testset "Custom Configuration Workflow" begin
        # Clear all kernels
        unload_all_spice_kernels()
        @test SPICE.ktotal("ALL") == 0
        
        # Load selective set
        load_spice_kernel("naif0012.tls")
        @test SPICE.ktotal("ALL") == 1
        
        load_spice_kernel(_EXTRA)
        @test SPICE.ktotal("ALL") == 2
        
        # Restore what AstroUniverse.__init__ loaded, read from DEFAULT_KERNELS. This is the
        # last testset in the file, so the state it leaves is what every later suite in the
        # process inherits — test_all_packages.jl runs all thirteen in one process. Reloading
        # a hand-written pair here left two kernels of four, and de440.bsp in place of the
        # merged epicycle_de440_1950-2100.bsp, which carries the satellite ephemerides. Plain
        # DE440 has planetary barycentres only, so body 499 (Mars) and 699 (Saturn) go missing
        # and every later frame test fails on "Insufficient ephemeris data".
        unload_all_spice_kernels()
        for k in AstroUniverse.DEFAULT_KERNELS
            load_spice_kernel(k)
        end
        @test SPICE.ktotal("ALL") == length(AstroUniverse.DEFAULT_KERNELS)
    end

end

unload_all_spice_kernels()   # releases the temporary copies, so they can be deleted
end
end

# Every later suite in the process reads the real cache again.
for k in AstroUniverse.DEFAULT_KERNELS
    load_spice_kernel(k)
end
