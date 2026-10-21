
library(shiny)
library(maftools)
library(DT)
library(plyr)
library(dplyr)
library(DT)
library(ggplot2)
library(ggvenn)

Sample.List = c('3-100X','3-200X','3-400X','17-100X','17-200X','17-400X','tissue-200X')

ui = fluidPage(
   tabsetPanel(
        # tabPanel( title = div(style = "font-size: 25px; ", 'Summary'),
        #     fluidRow(
        #         selectInput( inputId = 'SUM_SAMPLE', choices = c('all', Sample.List), selected = NULL, label = div(style = "font-family: Roboto Condensed; font-size: 14px; color: #47759f", 'Sample Selection') ),
        #     )
        # ),
        tabPanel( title = div(style = "font-size: 25px; ", 'Compare'),
            fluidRow(
                column( width = 3,
                    fluidRow(
                        column( width = 6, selectInput( inputId = 'S1', choice = Sample.List, selected = NULL, multiple = FALSE, label = div(style = "font-size: 12; color: #47759f",'Sample A')) ),    
                        column( width = 6, selectInput( inputId = 'S2', choice = NULL, selected = NULL, multiple = FALSE, label = div(style = "font-size: 12; color: #47759f",'Sample B')) )  
                    )
                ),
                column( width = 1),
                column( width = 4, fluidRow( plotOutput('MAF_SUM1'))),
                column( width = 4, fluidRow( plotOutput('MAF_SUM2')))
            ),
            hr(),
            fluidRow(
                column( width = 4,
                    fluidRow( selectInput( inputId='VisLevel', choices = c('Variant','Gene'), selected = NULL, label = div(style = "font-size: 14px; color: #47759f", 'Visual Level')) ),
                    hr(),
                    fluidRow( plotOutput('VENN2'))
                ),
                column( width = 4,
                    fluidRow( plotOutput('MUT_MATCH_RATES' ))
                ),
                column( width = 4, 
                    fluidRow( 
                        br(),
                        actionButton( inputId = 'CommonGeneUpdate', label = div(style = "font-size: 16px; color: #D8522F", 'GENES UPDATE') ),
                        br(),
                        selectInput( inputId = 'CommonGeneList', choices = NULL, selected = NULL, label = div(style = "font-size: 14px; color: #47759f", 'Gene Selection') )
                    ),
                    br(),
                    fluidRow( plotOutput('LOLLI2') )
                )
            )
        )
   ) 
)

server = function( input, output, session ) {
    
    SampleSelectUpdate = reactive({ SampleSelectUpdate = Sample.List[ Sample.List != input$S1 ] })

    observeEvent( input$S1, {
        updateSelectInput( session=session, inputId='S2', choices = SampleSelectUpdate(), selected = NULL ) 
    })

    maf.s1 = reactive({ read.maf( sprintf('OS_%s_maf.tsv', input$S1) ) })
    maf.s2 = reactive({ read.maf( sprintf('OS_%s_maf.tsv', input$S2) ) })

    commonGenes = reactive({ 
        maf1 = maf.s1()@data %>% mutate( varid = paste( Hugo_Symbol, cDNA_Change, sep="_")) 
        maf2 = maf.s2()@data %>% mutate( varid = paste( Hugo_Symbol, cDNA_Change, sep="_")) 
        common.genes = intersect(maf1$Hugo_Symbol, maf2$Hugo_Symbol) 
        commonGenes = names(rev(sort(table(maf.s1()@data %>% filter( Hugo_Symbol %in% common.genes ) %>% .$Hugo_Symbol))))
        return(commonGenes)
    })

    observeEvent( input$CommonGeneUpdate, {
        updateSelectInput( session=session, inputId='CommonGeneList', choices = commonGenes(), selected = NULL ) 
    }) 

    output$VENN2 = renderPlot({
        maf1 = maf.s1()@data %>% mutate( varid = paste( Hugo_Symbol, cDNA_Change, sep="_")) 
        maf2 = maf.s2()@data %>% mutate( varid = paste( Hugo_Symbol, cDNA_Change, sep="_")) 
        if( input$VisLevel == 'Variant' ){
            VennSetData = list( s1=maf1$cDNA_Change, s2=maf2$cDNA_Change )          
        }else{
            VennSetData = list( s1=unique(maf1$Hugo_Symbol), s2=unique(maf2$Hugo_Symbol) )
        }
        names(VennSetData) = c( input$S1, input$S2 )
        ggvenn( VennSetData,  show_percentage = T, set_name_size = 6, text_size = 5 )
    }, width=300, height=250)

    output$MUT_MATCH_RATES = renderPlot({
        maf1 = maf.s1()@data %>% mutate( varid = paste( Hugo_Symbol, cDNA_Change, sep="_")) %>% select(c('varid','Hugo_Symbol','cDNA_Change','tumor_f'))
        maf2 = maf.s2()@data %>% mutate( varid = paste( Hugo_Symbol, cDNA_Change, sep="_")) %>% select(c('varid','Hugo_Symbol','cDNA_Change','tumor_f'))
        common.variants = intersect(maf1$varid, maf2$varid)
        all.vairants = union(maf1$varid, maf2$varid)
        maf.join = maf1[match(all.vairants, maf1$varid), ] %>% rename( s1_vaf=tumor_f )
        maf.join = data.frame( maf.join, maf2[match(all.vairants, maf2$varid), 'tumor_f'] ) %>% rename( s2_vaf = tumor_f ) %>% arrange(s1_vaf) 
        maf.join[is.na(maf.join)] = 0
        maf.join$common=0 ; maf.join[which(maf.join$varid %in% common.variants), 'common'] = 1 
        maf.join$s1_vaf = round(maf.join$s1_vaf, 3 )
        maf.join$s2_vaf = round(maf.join$s2_vaf, 3 )
        maf.join$vaf_chk = apply( maf.join[,c('s1_vaf', 's2_vaf')], 1, function(v) ifelse( v[1] == v[2], 1, 0) )

        ggplot( maf.join, aes(x=s1_vaf, y=s2_vaf) ) + geom_point() +
            labs( x=input$S1, y=input$S2)
    })
    
    output$MAF_SUM1 = renderPlot({ plotmafSummary(maf.s1()) })
    output$MAF_SUM2 = renderPlot({ plotmafSummary(maf.s2()) })
    

    output$LOLLI2 = renderPlot({
        lollipopPlot2(m1 = maf.s1(), m2=maf.s2(), gene=input$CommonGeneList, m1_name=input$S1, m2_name=input$S2)
    })
}

shinyApp(ui, server)
