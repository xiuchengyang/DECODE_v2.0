# DECODE_v2.0 (Apr. 1, 2026)
DEtection and Characterization of cOastal tiDal wEtlands change (DECODE) 

Algorithm for DEtection and Characterization of cOastal tiDal wEtlands change (DECODE) using Landsat time series. It can provide changes in coastal tidal wetlands by providing accurate land cover and land change maps for intertidal ecosystems, fully automated for Landsat and Sentinel-2. 

**DECODE v2 [Continuously updated]** An updated version. 
**DECODE v1 [2024]** [Available at GERS Lab] (https://github.com/GERSL/DECODE). 

Please get in touch with Xiucheng Yang (xiucheng.yang@uconn.edu;xiuchengyang@uvic.ca), if you have any questions.

## TMD and EOT20

The Tide Model Driver functions and EOT20 ocean tide model are not included
in this repository. When `Main_DECODE_BatchRun.m` is run for the first time,
`DownloadTMDFunctions.m` downloads the official Tide Model Driver source and
`DownloadEOT20Model.m` downloads the EOT20 model. Both are saved under
`TMD_functions/`, which is ignored by Git.

Tide Model Driver source:
[https://github.com/chadagreene/Tide-Model-Driver](https://github.com/chadagreene/Tide-Model-Driver)

Please cite Tide Model Driver:

Greene, C. A., Erofeeva, S., Padman, L., Howard, S. L., Sutterley, T., and
Egbert, G. (2024). Tide Model Driver for MATLAB. *Journal of Open Source
Software*, 9(95), 6018.
[https://doi.org/10.21105/joss.06018](https://doi.org/10.21105/joss.06018).

Please cite EOT20:

Hart-Davis, M. G., Piccioni, G., Dettmering, D., Schwatke, C., Passaro, M.,
and Seitz, F. (2021). EOT20: a global ocean tide model from multi-mission
satellite altimetry. *Earth System Science Data*, 13, 3869-3884.
[https://doi.org/10.5194/essd-13-3869-2021](https://doi.org/10.5194/essd-13-3869-2021).

**Interactive maps (US Tidal Wetland Cover Change)** are available at [GEE APP for tidal wetland covers in the US](https://gers.users.earthengine.app/view/tidalwetlandcover). 

**Interactive maps (US Tidal marsh condition change and disturbance)** are available at [GEE APP for tidal marsh condition change in the US](https://ee-gers.projects.earthengine.app/view/ustidalmarsh). 

**Interactive maps (Florida Mangrove Resilience)** are available at [GEE APP for US Mangrove Disturbance, Dieback, and Recovery](https://gers.users.earthengine.app/view/tidalwetlandcover). 

**Interactive maps (Canada Tidal Marsh)** are available at [GEE APP for tidal marsh dynamics in Canada](https://xiucheng.projects.earthengine.app/view/canadamarsh). 

**Please cite the following paper** 

Xiucheng Yang*, Shi Qiu, Kevin D. Kroeger, Scott Covington, Zhiliang Zhu, Nicholas J. Murray, and Zhe Zhu*. The accelerating loss and shifting dynamics of US tidal wetlands. Nature Communications(2026). [DECODE on US Wetland Cover Change Analysis]

Xiucheng Yang*, Zhe Zhu, Kevin D. Kroeger, Shi Qiu, Scott Covington, Jeremy R. Conrad, and Zhiliang Zhu. Tracking mangrove condition changes using dense Landsat time series. Remote Sensing of Environment 276 (2024): 114461. [https://doi.org/10.1016/j.rse.2024.114461](https://doi.org/10.1016/j.rse.2024.114461). [DECODER -- DECODE and Recovery/Resilience]

Xiucheng Yang*, Zhe Zhu, Shi Qiu, Kevin D. Kroeger, Zhiliang Zhu, and Scott Covington. "Detection and characterization of coastal tidal wetland change in the northeastern US using Landsat time series." Remote Sensing of Environment 276 (2022): 113047. [https://doi.org/10.1016/j.rse.2022.113047](https://doi.org/10.1016/j.rse.2022.113047). [DECODE v1]
