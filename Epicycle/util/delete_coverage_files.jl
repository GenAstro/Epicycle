# Delete all .cov files in the repository
# This script removes code coverage files generated during testing

# Go up two levels: util -> Epicycle -> repository root
repo_root = dirname(dirname(dirname(@__FILE__)))

println("Searching for .cov files in: $repo_root")

cov_files = String[]
for (root, dirs, files) in walkdir(repo_root)
    for file in files
        if endswith(file, ".cov")
            push!(cov_files, joinpath(root, file))
        end
    end
end

if isempty(cov_files)
    println("No .cov files found.")
else
    println("Found $(length(cov_files)) .cov file(s). Deleting...")
    
    for file in cov_files
        println("  Deleting: $file")
        rm(file; force=true)
    end
    
    println("Done. Deleted $(length(cov_files)) file(s).")
end
