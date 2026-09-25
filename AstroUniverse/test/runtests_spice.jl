# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

using Test
using AstroUniverse
using SPICE
using Logging

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
        # Download de440s.bsp for testing
        download_spice_kernel("de440s.bsp",
            "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp")
        
        cache_dir = get_spice_directory()
        kernel_path = joinpath(cache_dir, "de440s.bsp")
        
        # Verify file was downloaded
        @test isfile(kernel_path)
        
        # Verify file size is reasonable (de440s.bsp should be ~31 MB)
        @test filesize(kernel_path) > 1_000_000  # At least 1 MB
        
        # Test idempotency - re-downloading should not error
        @test_nowarn download_spice_kernel("de440s.bsp",
            "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp")
    end
    
    @testset "get_spice_directory()" begin
        dir = get_spice_directory()
        @test isdir(dir)
        @test contains(dir, "spice_kernels")
    end
    
    @testset "list_downloaded_spice_kernels()" begin
        # Should not error and should display downloaded kernel info
        @test_nowarn list_downloaded_spice_kernels()
        
        # Verify default kernels are in the storage directory
        storage_dir = get_spice_directory()
        downloaded_files = readdir(storage_dir)
        @test "naif0012.tls" in downloaded_files
        @test "epicycle_de440_1950-2100.bsp" in downloaded_files
        @test "de440s.bsp" in downloaded_files  # From earlier test
    end
    
    @testset "list_cached_spice_kernels()" begin
        # Should not error - shows what's loaded in SPICE memory
        @test_nowarn list_cached_spice_kernels()
        
        # Ensure de440s.bsp is NOT loaded to start with clean state
        try
            unload_spice_kernel("de440s.bsp")
        catch
            # Ignore if not already loaded
        end
        
        # After loading a kernel, verify it appears in the kernel pool
        download_spice_kernel("de440s.bsp", "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp")
        initial_count = SPICE.ktotal("ALL")
        load_spice_kernel("de440s.bsp")
        @test SPICE.ktotal("ALL") == initial_count + 1
        
        # Verify de440s.bsp is in the cached (loaded) list
        cached_kernels = String[]
        for i in 1:SPICE.ktotal("ALL")
            result = SPICE.kdata(i, "ALL")
            if result !== nothing
                push!(cached_kernels, basename(result[1]))
            end
        end
        @test "de440s.bsp" in cached_kernels
        for filename in AstroUniverse.DEFAULT_KERNELS
            @test filename in cached_kernels
        end
        
        # Clean up
        unload_spice_kernel("de440s.bsp")
    end
    
    @testset "Load/Unload Operations" begin
        # Ensure de440s.bsp is NOT loaded to start with clean state
        try
            unload_spice_kernel("de440s.bsp")
        catch
            # Ignore if not already loaded
        end
        
        # Count kernels before operations
        initial_count = SPICE.ktotal("ALL")
        
        # Load de440s.bsp
        @test_nowarn load_spice_kernel("de440s.bsp")
        @test SPICE.ktotal("ALL") == initial_count + 1
        
        # Verify de440s.bsp is actually in the kernel pool
        loaded_kernels = String[]
        for i in 1:SPICE.ktotal("ALL")
            result = SPICE.kdata(i, "ALL")
            if result !== nothing
                push!(loaded_kernels, basename(result[1]))
            end
        end
        @test "de440s.bsp" in loaded_kernels
        
        # Try to load again - should warn and not increase count
        @test_logs (:warn, r"already loaded") load_spice_kernel("de440s.bsp")
        @test SPICE.ktotal("ALL") == initial_count + 1
        
        # Unload de440s.bsp
        @test_nowarn unload_spice_kernel("de440s.bsp")
        @test SPICE.ktotal("ALL") == initial_count
        
        # Verify de440s.bsp is no longer in the kernel pool
        loaded_kernels = String[]
        for i in 1:SPICE.ktotal("ALL")
            result = SPICE.kdata(i, "ALL")
            if result !== nothing
                push!(loaded_kernels, basename(result[1]))
            end
        end
        @test !("de440s.bsp" in loaded_kernels)
        
        # Try to unload again - should warn
        @test_logs (:warn, r"not currently loaded") unload_spice_kernel("de440s.bsp")
        @test SPICE.ktotal("ALL") == initial_count
        
        # Load it again for later tests
        load_spice_kernel("de440s.bsp")
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
        # Ensure clean state - unload de440s.bsp from previous test, ensure epicycle_de440_1950-2100.bsp is loaded
        try
            unload_spice_kernel("de440s.bsp")
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
        
        # Load de440s.bsp instead
        load_spice_kernel("de440s.bsp")
        @test SPICE.ktotal("ALL") == initial_count
        
        # Verify de440s.bsp is in the pool
        loaded_kernels = String[]
        for i in 1:SPICE.ktotal("ALL")
            result = SPICE.kdata(i, "ALL")
            if result !== nothing
                push!(loaded_kernels, basename(result[1]))
            end
        end
        @test "de440s.bsp" in loaded_kernels
        @test !("epicycle_de440_1950-2100.bsp" in loaded_kernels)
        
        # Swap back
        unload_spice_kernel("de440s.bsp")
        load_spice_kernel("epicycle_de440_1950-2100.bsp")
        @test SPICE.ktotal("ALL") == initial_count
        
        # Verify final state: epicycle_de440_1950-2100.bsp is in pool, de440s.bsp is not
        loaded_kernels = String[]
        for i in 1:SPICE.ktotal("ALL")
            result = SPICE.kdata(i, "ALL")
            if result !== nothing
                push!(loaded_kernels, basename(result[1]))
            end
        end
        @test "epicycle_de440_1950-2100.bsp" in loaded_kernels
        @test !("de440s.bsp" in loaded_kernels)
    end
    
    @testset "Custom Configuration Workflow" begin
        # Clear all kernels
        unload_all_spice_kernels()
        @test SPICE.ktotal("ALL") == 0
        
        # Load selective set
        load_spice_kernel("naif0012.tls")
        @test SPICE.ktotal("ALL") == 1
        
        load_spice_kernel("de440s.bsp")
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
