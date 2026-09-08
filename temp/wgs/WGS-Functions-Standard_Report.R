
### 
    options(stringsAsFactors=FALSE)
    suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("Hmisc"))
    suppressPackageStartupMessages(library("RMySQL"))
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("gt"))
    suppressPackageStartupMessages(library("Cairo"))
    options(bitmapType='cairo')
    grDevices::X11.options(type='cairo')
    options(device='x11')
###
### PAGE BREAKS ---------------------------------------------------------------- 
    NGS_report.Page.Breaks2 <- function(){
        cat('<div style="page-break-after: always;"></div>')
        cat('<br>')
    }
###
### 
    GT_DefaultForm <- function(gt_tbl) {
        gt_tbl |> tab_style(
            style     = cell_text( weight = "bold", font="Roboto", size=px(11) ),
            locations = cells_column_labels( everything() )
        ) |> tab_style(
            style     = cell_text( font="Roboto Condensed", size=px(10) ),
            locations = cells_body()
        ) |> cols_align( 
            align = "center", columns = everything()
        ) |> tab_options(
            table.width                    = pct(100),
            data_row.padding               = px(4),
            heading.title.font.size        = px(11),
            heading.subtitle.font.size     = px(11),
            table.border.top.style         = "none",
            column_labels.border.top.style = "none",
            heading.border.bottom.style    = "none"
        )
    }
###
    GT_BoBorderForm <- function(gt_tbl){
        gt_tbl |> tab_options(
            table.width                       = pct(100),
            table.font.names                  = "Roboto",
            table.font.size                   = px(11),
            data_row.padding                  = px(5),
            table.border.top.style            = "none",
            table.border.bottom.style         = "none",
            table_body.border.top.style       = "none",
            table_body.border.bottom.style    = "none",
            table_body.hlines.style           = "none",
            table_body.vlines.style           = "none",
            column_labels.border.top.style    = "none",
            column_labels.border.bottom.style = "none",
            column_labels.hidden              = TRUE,
            row_group.border.top.style        = "none",
            row_group.border.bottom.style     = "none"
            
        )
    }    

