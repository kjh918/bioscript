
# WES CNV Normal Panels Creation

### WES PON list file 
- /data/wes/params/cnv_pon/GCX_WES_Panel_of_Normals_for_CNV_rev.xlsx

### CNV PON Create Script Usage 
```
    Rscript /storage/home/kangsm/runScripts/WES_run.CNV_Panel.of.Normal_Preparation.R \
        --PON_LIST_XLSX ${PON_LIST_FILE} \
        --NGS_LIB ${NGS_LIB} \
        --PON_PURECN TRUE \
        --PON_SUPERFREQ TRUE \
        --PON_ASCAT TRUE
```
<hr>

### NGS_LIB : `twist.exome.2.0` ------------------------------------------------
- last update : 2025-05-30
- samples     : 26

```
    Rscript /storage/home/kangsm/runScripts/WES_run.CNV_Panel.of.Normal_Preparation.R \
        --PON_LIST_XLSX "/data/wes/params/cnv_pon/GCX_WES_Panel_of_Normals_for_CNV_rev.xlsx" \
        --NGS_LIB "twist.exome.2.0" \
        --PON_PURECN TRUE \
        --PON_SUPERFREQ TRUE \
        --PON_ASCAT TRUE
```

### NGS_LIB : `sureselect.v6` --------------------------------------------------

- last update : 2025-05-30
- samples     : 13

```
    Rscript /storage/home/kangsm/runScripts/WES_run.CNV_Panel.of.Normal_Preparation.R \
        --PON_LIST_XLSX "/data/wes/params/cnv_pon/GCX_WES_Panel_of_Normals_for_CNV_rev.xlsx" \
        --NGS_LIB "sureselect.v6" \
        --PON_PURECN TRUE \
        --PON_SUPERFREQ TRUE \
        --PON_ASCAT TRUE
```

### NGS_LIB : `sureselect.utr.v4` ----------------------------------------------

- last update : 2025-05-30
- samples     : 18

```
    Rscript /storage/home/kangsm/runScripts/WES_run.CNV_Panel.of.Normal_Preparation.R \
        --PON_LIST_XLSX "/data/wes/params/cnv_pon/GCX_WES_Panel_of_Normals_for_CNV_rev.xlsx" \
        --NGS_LIB "sureselect.utr.v4" \
        --PON_PURECN TRUE \
        --PON_SUPERFREQ TRUE \
        --PON_ASCAT TRUE
```

### NGS_LIB : `xgen.exome` ----------------------------------------------

- last update : 2025-05-30
- samples     : 3

```
    Rscript /storage/home/kangsm/runScripts/WES_run.CNV_Panel.of.Normal_Preparation.R \
        --PON_LIST_XLSX "/data/wes/params/cnv_pon/GCX_WES_Panel_of_Normals_for_CNV_rev.xlsx" \
        --NGS_LIB "xgen.exome" \
        --PON_PURECN TRUE \
        --PON_SUPERFREQ TRUE \
        --PON_ASCAT TRUE
```

### NGS_LIB : `sureselect.xt.v8` ----------------------------------------------

- last update : 2025-05-30
- samples     : 5

```
    Rscript /storage/home/kangsm/runScripts/WES_run.CNV_Panel.of.Normal_Preparation.R \
        --PON_LIST_XLSX "/data/wes/params/cnv_pon/GCX_WES_Panel_of_Normals_for_CNV_rev.xlsx" \
        --NGS_LIB "sureselect.xt.v8" \
        --PON_PURECN TRUE \
        --PON_SUPERFREQ TRUE \
        --PON_ASCAT TRUE
```

### NGS_LIB : `rapid.capture.extended.exome` ----------------------------------------------

- last update : 2025-05-30
- samples     : 13

```
    Rscript /storage/home/kangsm/runScripts/WES_run.CNV_Panel.of.Normal_Preparation.R \
        --PON_LIST_XLSX "/data/wes/params/cnv_pon/GCX_WES_Panel_of_Normals_for_CNV_rev.xlsx" \
        --NGS_LIB "rapid.capture.extended.exome" \
        --PON_PURECN TRUE \
        --PON_SUPERFREQ TRUE \
        --PON_ASCAT TRUE
```




