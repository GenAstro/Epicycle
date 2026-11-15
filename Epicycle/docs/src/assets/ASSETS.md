# Documentation assets

## GEOTransfer.png

A render of the GEO transfer trajectory, used as the header image on the documentation index. It is
Epicycle's own output, produced by the package's 3D view, so the image itself is ours.

It is composed from two third-party assets, and it is only as clear as they are.

- The Deep Space 1 spacecraft in the foreground is the model in `Epicycle/assets/`, from NASA 3D
  Resources. See `../../assets/ASSETS.md`.
- The Earth surface imagery on the globe is `AstroUniverse/data/EarthTexture.jpg`, from NASA Earth
  Observatory's Blue Marble collection. See `../../../../AstroUniverse/data/ASSETS.md`.

## MainPageTrajectoryGraphic.png

Epicycle's own output, produced by the package's 3D view.

## EpicycleProductBrochure.pdf

The Gen Astro product brochure for the R2026 beta, linked from the documentation index. The large
images are re-encoded as JPEG at quality 85, keeping their transparency masks, the small icons and
the text as they are, which brings the 16.6 MB export to 3.0 MB. Replace it under the same name so
links already sent keep working. It shows the VS Code, Cesium, Plotly and Julia logos to name the tools Epicycle
works with; where the background photography comes from is recorded with the brochure's source
files, not here.

## EpicycleProductBrochureCover.png

The brochure's first page, rendered at 368 by 476 pixels, used as the link image on the index.

## EpicycleODBrochure.pdf

The three-page precision orbit determination brochure, linked from the documentation index. Built
from the HTML pages in `od_runs/carousel` (outside the repository) and copied here unchanged.
Replace it under the same name so links already sent keep working.

## EpicycleODBrochureCover.png

The OD brochure's first page, rendered at 368 by 476 pixels with pdftoppm, used as the link image
on the index.

## gallery

Epicycle's own output: each image is what an example script in `Epicycle/examples/` draws.

## gallery/thumbs

Thumbnails for the gallery strip on the index: one panel cropped from the gallery image of the
same name and downscaled to 600 by 450 pixels, two to three times their displayed width so they
stay sharp on high-density screens. Regenerate them if the source image changes.

`LunarTransferBPlane.png` is not cropped from a gallery image: it is the Moon-centered pane of
the frame at 18.4 s of `LunarTransferIterationsAndBPlane.mp4` (GenAstro Graphics/Videos), cropped
to 4:3 and scaled to 600 by 450.

`PrecisionOD.png` is the four mission tiles (LAGEOS, LRO, MESSENGER, Juno) from page 1 of
`EpicycleODBrochure.pdf`, rendered at 400 dpi, laid out 2 × 2 on a dark background and scaled to
600 by 450. Rebuild it if the brochure changes.
