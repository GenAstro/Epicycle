using FileIO
using Images

# Load the original Milky Way image
#milkyway_raw = load(raw"C:\Users\steve\Dev\Epicycle\Epicycle\prototype\MilkyWay.jpg")
milkyway_raw = load(raw"C:\Users\steve\Dev\Epicycle\Epicycle\prototype\MilkyWay.jpg")

# Filter parameters
brightness_threshold = 0.15  # Pixels dimmer than this become black
brightness_multiplier = 0.5   # Dim the bright features by this amount

# Filter the image
milkyway_filtered = map(milkyway_raw) do pixel
    # Calculate brightness (luminance)
    brightness = 0.299 * red(pixel) + 0.587 * green(pixel) + 0.114 * blue(pixel)
    
    if brightness < brightness_threshold
        # Replace dim pixels/stars with black
        RGB(0.0, 0.0, 0.0)
    else
        # Keep bright pixels (galactic cloud, bright stars), slightly dimmed
        RGB(red(pixel) * brightness_multiplier, 
            green(pixel) * brightness_multiplier, 
            blue(pixel) * brightness_multiplier)
    end
end

# Save the filtered image
save(raw"C:\Users\steve\Dev\Epicycle\Epicycle\prototype\MilkyWay_Filtered.jpg", milkyway_filtered)

println("Filtered Milky Way image saved as MilkyWay_Filtered.jpg")
println("Adjust brightness_threshold (currently $brightness_threshold) or brightness_multiplier (currently $brightness_multiplier) if needed")
