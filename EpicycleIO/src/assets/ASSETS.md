# Third-party assets

Two binary files are shipped with EpicycleIO and staged into the served directory at run time.
Both are believed to be United States Government works produced by NASA, which are not subject to
copyright protection in the United States and may be redistributed. This file records what is
known about each, because the package is published and an asset without a recorded source cannot
be cleared later by looking at it.

NASA media usage guidelines: <https://www.nasa.gov/nasa-brand-center/images-and-media/>

## Aqua.glb

A model of the Aqua spacecraft, used as the entity mesh in a 3D view.

- Source: NASA 3D Resources, <https://science.nasa.gov/3d-resources/>
- Size: 1 104 016 bytes
- SHA-256: `f0761181a1d378b6b3e270e8eee8bf20c241b7d577933630369b0cffd4b521e1`

NASA 3D Resources publishes models for public use under the NASA media usage guidelines above.

## black_marble.jpg

An equirectangular image of Earth at night, used as the night-side imagery layer on the globe.

- Size: 794 479 bytes
- Dimensions: 3600 x 1800, equirectangular, JPEG
- SHA-256: `373e5a08c9f378a2ce6320214a613148e4b1e3946b3f39a516c9093b76cb7124`

🔴 **The source is not recorded and has not been verified.** The file was added in commit
`52b5da4`, a bulk snapshot whose message does not mention it, and no download URL appears
anywhere in the repository history.

The likely source is NASA Visible Earth, "Earth at Night (Black Marble) 2016", which is published
at exactly 3600 x 1800 in this projection:
<https://visibleearth.nasa.gov/images/144898/earth-at-night-black-marble-2016-color-maps>

The dimensions match and the file size is consistent with the colour map at that resolution, but
matching dimensions is not provenance. Confirming it means downloading the candidate and
comparing against the SHA-256 above.

Until that comparison is made, treat this asset as unverified. If it cannot be matched to a NASA
product, replace it with one that can rather than publish it on an assumption.
