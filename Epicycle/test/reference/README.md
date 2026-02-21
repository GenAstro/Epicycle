# Visual Regression Test References

This directory contains reference images for visual regression testing of the graphics module.

## Purpose

These images serve as "golden masters" - the expected visual output of rendering tests. When visual regression tests run, they compare new renders against these references pixel-by-pixel.

## Workflow

### First Time Setup
```bash
# Run visual tests - creates reference images
julia --project=. test/graphics/visual_regression.jl

# Manually inspect images in this directory
# If they look correct, commit them
git add test/reference/*.png
git commit -m "Add visual regression references"
```

### Before Each Release
```bash
# Run tests - compares against references
julia --project=. test/graphics/visual_regression.jl

# If tests fail:
# 1. Check if it's a real regression (bug) - fix the code
# 2. Check if it's an intentional change - update references:
#    - Delete old reference: rm test/reference/<test_name>.png
#    - Re-run test to create new reference
#    - Review carefully and commit
```

## Image List

- `leo_orbit.png` - Simple circular LEO orbit
- `elliptical_orbit.png` - High eccentricity GTO-like orbit  
- `multi_segment.png` - Trajectory with maneuver (multiple segments)
- `high_inc_orbit.png` - Polar-like orbit

## Notes

- These tests require GLMakie with GPU/display access
- Run locally, not in CI
- Images are platform-dependent (GPU, driver differences)
- Small pixel differences are normal - ReferenceTests handles thresholds
