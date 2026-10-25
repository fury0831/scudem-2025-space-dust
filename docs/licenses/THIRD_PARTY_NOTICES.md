# Third-party notices

The three original HTML files embed Plotly.js **2.35.2**. Their embedded library bytes match the official version (after line-ending normalization). The files and saved data have not been numerically changed.

The team's repository MIT license applies to team-authored implementation and documentation. It does **not** replace the licenses of Plotly or its bundled components.

- [Plotly MIT license](PLOTLY-LICENSE.txt)
- [Official bundle sidecar notices](../plotly.min.js.LICENSE.txt), retained at the name referenced by the embedded bundle
- [MapLibre GL JS 4.5.2 license](MAPLIBRE-LICENSE.txt), including its listed components
- [ieee754 BSD-3-Clause license](IEEE754-LICENSE.txt)
- Additional upstream MIT license texts for polybooljs, object-assign, buffer, safe-buffer, is-buffer, pad-left, and repeat-string are in this directory.
- Native Promise Only 0.8.1: copyright Kyle Simpson, MIT; its original notice and upstream license link are retained in the official bundle sidecar. Its published package does not include a separate LICENSE file; no replacement copyright statement has been invented.

[Source URLs, hashes, versions, and retrieval notes](sources.json) record the upstream material. Exact bundled versions are noted where the upstream lockfile identifies them; otherwise the file is marked as an upstream license snapshot.

The standard Plotly toolbar includes an optional external link to plotly.com. No custom network requests, external scripts, forms, credentials, or private identifiers were found in the project-specific figure code. These scatter plots do not use the map/remote-data functionality also present in the full Plotly bundle.
