# DECODE v2

DEtection and Characterization of cOastal tiDal wEtlands change (DECODE) is an
algorithm for detecting and characterizing coastal tidal wetland change from
dense satellite time series. It supports Landsat and Sentinel-2 observations
and accounts for tidal effects during time-series modeling.

## Versions

- **DECODE v2 (continuously updated):** Current version.
- **DECODE v1 (2024):** [Available from the GERS Lab](https://github.com/GERSL/DECODE).

## Contact

For questions, contact Xiucheng Yang at xiucheng.yang@uconn.edu or
xiuchengyang@uvic.ca.

## Citation

**[DECODE]** Yang, X., Zhu, Z., Qiu, S., Kroeger, K. D., Zhu, Z., and Covington, S. (2022).
Detection and characterization of coastal tidal wetland change in the
northeastern US using Landsat time series. *Remote Sensing of Environment*,
276, 113047. [https://doi.org/10.1016/j.rse.2022.113047](https://doi.org/10.1016/j.rse.2022.113047).

**[US Application]** Yang, X., Qiu, S., Kroeger, K. D., Zhu, Z., Covington, S., Murray, N. J., and
Zhu, Z. (2026). The accelerating loss and shifting dynamics of US tidal
wetlands. *Nature Communications*, 17, 4332.
[https://doi.org/10.1038/s41467-026-71464-2](https://doi.org/10.1038/s41467-026-71464-2).

Interactive map: [US tidal wetland dynamics](https://gers.users.earthengine.app/view/tidalwetlandcover).

**[DECODE v2]** Yang, X., Knox, S. H., Qiu, S., Kroeger, K. D., Zhu, Z., Covington, S., and
Zhu, Z. (2026). Tracking US Tidal Marsh Extent and Change with Dense Landsat
Time Series and Comparisons with Existing Products. *Remote Sensing of
Environment* (under review).

Interactive map: [US tidal marsh condition change](https://ee-gers.projects.earthengine.app/view/ustidalmarsh).

**[DECODE+AlphaEarth]** Yang, X., Li, M., Suh, J. W., Zhu, Z., Baum, J. K., and Knox, S. H. (2026).
From Extent to Change: AlphaEarth and Time Series Analysis for Tidal Marsh
Monitoring. Preprint available at
[https://doi.org/10.2139/ssrn.7317755](https://doi.org/10.2139/ssrn.7317755).

Interactive map: [US and Canada tidal marsh extent at 10 m](https://xiucheng.projects.earthengine.app/view/marsh10m).

## Dependencies

DECODE uses the Tide Model Driver (TMD) and the EOT20 ocean tide model. These
third-party files are not included in this repository.

When `Main_DECODE_BatchRun.m` is run, `DownloadTMDFunctions.m` checks for the
TMD functions and `DownloadEOT20Model.m` checks for the EOT20 model. Existing
files are reused. Missing files are downloaded and saved under
`TMD_functions/`, which is ignored by Git.

Tide Model Driver source:
[https://github.com/chadagreene/Tide-Model-Driver](https://github.com/chadagreene/Tide-Model-Driver)

Please cite the Tide Model Driver and EOT20:

Greene, C. A., Erofeeva, S., Padman, L., Howard, S. L., Sutterley, T., and
Egbert, G. (2024). Tide Model Driver for MATLAB. *Journal of Open Source
Software*, 9(95), 6018.
[https://doi.org/10.21105/joss.06018](https://doi.org/10.21105/joss.06018).

Hart-Davis, M. G., Piccioni, G., Dettmering, D., Schwatke, C., Passaro, M., and
Seitz, F. (2021). EOT20: a global ocean tide model from multi-mission satellite
altimetry. *Earth System Science Data*, 13, 3869-3884.
[https://doi.org/10.5194/essd-13-3869-2021](https://doi.org/10.5194/essd-13-3869-2021).
