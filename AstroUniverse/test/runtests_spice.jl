# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

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
        @test isfile(joinpath(cache_dir, "de440.bsp"))
        
        # Verify kernels are actually cached (loaded in SPICE memory)
        loaded_kernels = String[]
        for i in 1:SPICE.ktotal("ALL")
            result = SPICE.kdata(i, "ALL")
            if result !== nothing
                push!(loaded_kernels, basename(result[1]))
            end
        end
        @test "naif0012.tls" in loaded_kernels
        @test "de440.bsp" in loaded_kernels
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
        @test "de440.bsp" in downloaded_files
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
        @test "naif0012.tls" in cached_kernels
        @test "de440.bsp" in cached_kernels
        
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
        
        # Reload defaults for remaining tests
        load_spice_kernel("naif0012.tls")
        load_spice_kernel("de440.bsp")
        @test SPICE.ktotal("ALL") == 2
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
        # Ensure clean state - unload de440s.bsp from previous test, ensure de440.bsp is loaded
        try
            unload_spice_kernel("de440s.bsp")
        catch
        end
        
        # Make sure de440.bsp is loaded
        try
            load_spice_kernel("de440.bsp")
        catch
        end
        initial_count = SPICE.ktotal("ALL")
        
        # Unload de440.bsp
        unload_spice_kernel("de440.bsp")
        @test SPICE.ktotal("ALL") == initial_count - 1
        
        # Verify de440.bsp is no longer in the pool
        loaded_kernels = String[]
        for i in 1:SPICE.ktotal("ALL")
            result = SPICE.kdata(i, "ALL")
            if result !== nothing
                push!(loaded_kernels, basename(result[1]))
            end
        end
        @test !("de440.bsp" in loaded_kernels)
        
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
        @test !("de440.bsp" in loaded_kernels)
        
        # Swap back
        unload_spice_kernel("de440s.bsp")
        load_spice_kernel("de440.bsp")
        @test SPICE.ktotal("ALL") == initial_count
        
        # Verify final state: de440.bsp is in pool, de440s.bsp is not
        loaded_kernels = String[]
        for i in 1:SPICE.ktotal("ALL")
            result = SPICE.kdata(i, "ALL")
            if result !== nothing
                push!(loaded_kernels, basename(result[1]))
            end
        end
        @test "de440.bsp" in loaded_kernels
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
        
        # Restore defaults
        unload_all_spice_kernels()
        load_spice_kernel("naif0012.tls")
        load_spice_kernel("de440.bsp")
    end
    
end
