# Third-party notices

AstroEpochs contains code adapted from the packages below. Each adapted file names its source in
its header. The adapted code is held to its source's own outputs: `test/test_correctness_erfa_parity.jl`
compares the ERFA ports with pyerfa exactly, and `test/test_correctness_astropy_benchmark.jl`
compares the assembled conversions with Astropy at 1668 epochs.

## ERFA

ERFA 2.0.1, https://github.com/liberfa/erfa, the library Astropy computes time with. ERFA is
derived from the IAU SOFA library; AstroEpochs uses code derived from ERFA, not SOFA itself.

| AstroEpochs file | Ported from | What |
|---|---|---|
| `src/scales/offsets.jl` | `taitt.c`, `tttai.c`, `tttcg.c`, `tcgtt.c`, `tttdb.c`, `tdbtt.c`, `tdbtcb.c`, `tcbtdb.c` | the transforms between adjacent scales |
| `src/scales/tdb.jl` | `dtdb.c` | TDB − TT, the Fairhead & Bretagnon series, with its 787-row coefficient table copied by script |
| `src/scales/leap_seconds.jl` | `utctai.c`, `taiutc.c` | UTC ↔ TAI, with UTC as ERFA's quasi-Julian date on leap-second days |
| `src/formats/calendar.jl` | `cal2jd.c`, `jd2cal.c`, `d2tf.c`, `d2dtf.c`, `dtf2d.c` | the Gregorian calendar, time of day, and ISOT dates with leap seconds |

```
Copyright (C) 2013-2021, NumFOCUS Foundation.
All rights reserved.

This library is derived, with permission, from the International
Astronomical Union's "Standards of Fundamental Astronomy" library,
available from http://www.iausofa.org.

The ERFA version is intended to retain identical
functionality to the SOFA library, but made distinct through
different function and file names, as set out in the SOFA license
conditions. The SOFA original has a role as a reference standard
for the IAU and IERS, and consequently redistribution is permitted only
in its unaltered state. The ERFA version is not subject to this
restriction and therefore can be included in distributions which do not
support the concept of "read only" software.

Although the intent is to replicate the SOFA API (other than replacement of
prefix names) and results (with the exception of bugs; any that are
discovered will be fixed), SOFA is not responsible for any errors found
in this version of the library.

If you wish to acknowledge the SOFA heritage, please acknowledge that
you are using a library derived from SOFA, rather than SOFA itself.


TERMS AND CONDITIONS

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1 Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

2 Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

3 Neither the name of the Standards Of Fundamental Astronomy Board, the
   International Astronomical Union nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

## Astropy

Astropy 8.0.1, https://github.com/astropy/astropy. AstroEpochs normalises a two-part date as
Astropy's `Time` does: `_rebalance` and `_two_sum` in `src/AstroEpochs.jl` are ported from
`day_frac` and `two_sum` in `astropy/time/utils.py`. The scale graph's routes (`MULTI_HOPS` in
`src/scales/graph.jl`) are Astropy's, from `astropy/time/core.py`.

```
Copyright (c) 2011-2026, Astropy Developers

All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.
* Neither the name of the Astropy Team nor the names of its contributors may be
  used to endorse or promote products derived from this software without
  specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

## Tempo.jl

Tempo.jl v1.3.1, from JuliaSpaceMissionDesign, https://github.com/JuliaSpaceMissionDesign/Tempo.jl.
AstroEpochs depended on it until its code was brought in-house, because Tempo's dependency on
JSMDUtils held the Epicycle stack at ForwardDiff 0.10. The offsets and calendar routines first
adapted from it have since been replaced by the ERFA ports above; what remains is the built-in
leap-second table in `src/scales/leap_seconds.jl`, used only when no IERS list has been downloaded.

```
MIT License

Copyright (c) 2022 Andrea Pasquale <andrea.pasquale@polimi.it> and Michele
Ceresoli <michele.ceresoli@polimi.it>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
