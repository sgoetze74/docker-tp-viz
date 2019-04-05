
### go.r GO analysis ####
## Patrick Pedrioli
##
## VERSION: 0.2 (14-05-07)
## - Added support for subset ontology field (note that reformat_obo.pl was updated as well)
## - Added function to represent terms in word cloud
## - Added multiple testing correction
## - Added ability to limit go.gethyper to a specific GO namespace
## - Added ability to load untrimmed annotation files
##
## VERSION: 0.1
##
## gene_association files can be downloaded from:
## http://www.geneontology.org/GO.downloads.annotations.shtml
## Remove the initial comments
##
## goOBO file can be downloaded from:
## http://www.geneontology.org/GO.downloads.ontology.shtml
## WARNING: It will need processing with:
##       reformat_obo.pl gene_ontology_ext.obo > go_obo
##
##
## Example
## source( "./go.r" )
## go <- go.loadgo()
## up <- read.csv( "up.csv" , stringsAsFactor=F )
## up <- unlist( up , use.names=F )
## bg <- read.csv( "urm1_bg.csv" , colClasses = "character" , stringsAsFactors=F )
## bg <- unlist( bg , use.names=F )
## up.go <- go.gethyper( go , bg , up , "biological_process" , "fdr" )
## go.wordcloud( up.go , p.cutoff=0.01 )
## write.csv( up.go , quote=T , "up_go.csv" )


library( reshape2 )
library( plyr )
library( data.table )

#####
## Load GO data
go.loadgo <- function( species="Hsap" ) {

    if( species == "Hsap" ) {
        go_file_name="/usr/local/data/goa_human.gaf"
        go_obo_file_name="/usr/local/data/go_obo"
    }
    ## else if( species == "Mmus" ) {
    ##     go_file_name="GO/Data/Mmus/gene_association_trim.mgi"
    ##     go_obo_file_name="GO/Data/go_obo"
    ## }
    else {
        print( "Unsupported species!" )
        exit(-1)
    }
    
    ## Load gene <> GO
    go <- read.delim( go_file_name , header=FALSE , comment.char="!" )

    ## ORF Names are in column 11 separated from aliases by |
    ## Split them out
    gene.names <- strsplit( as.character( go[,11] ) , split="\\|" )
    go$orf <- sapply( gene.names , function(x) x[1] )

    ## Remove everything, but the ORF name and the GO ID
    go.f <- data.frame( go[,c(18,5)] )
    names( go.f ) <- c(  "ORF" , "GO" )

    ## Load GO term definitions
    go.def <- read.delim( go_obo_file_name, header=TRUE )

    go.m <- merge( go.f , go.def , by="GO" , all.x=TRUE )
    
    return( go.m )
}


#####
## Calculates hypergeometric pdf for terms in cluster
##
## go:          output of go.loadgo()
##
## all:         list of all genes considered in the experiment.
##
## cluster:     list of all significant proteins in the experiment.
##
## namespace:   biological_process , molecular_function , cellular_component
##
## adjust.method:       probability adjustment method
##
go.gethyper <- function( go , all , cluster , namespace="" , adjust.method="BH" ) {
    if( namespace != "" ){
        go <- go[ grep( namespace , go$namespace ) , ] 
    }

    genes.all <- all
    genes.cluster <- cluster
    
    go.all.subset <- go[ which( go$ORF %in% genes.all ) , ]
    m <- melt( go.all.subset , id= c( "GO" , "name" , "namespace" , "def" , "subset" ) )
    genome.all.gocount <- dcast( m , ... ~ variable , function(x) length( unique( x )))
    names( genome.all.gocount ) <- c( "GO" , "name" , "namespace" , "def" , "subset" , "all" )

    genome.cluster.size <- length( genes.cluster )
    go.cluster.subset <- go[ which( go$ORF %in% genes.cluster ) , ]
    m <- melt( go.cluster.subset , id= c( "GO" , "name" , "namespace" , "def" , "subset" ) )
    genome.cluster.gocount <- dcast( m , GO ~ variable , function(x) length( unique( x )))
    names( genome.cluster.gocount ) <-  c( "GO" , "cluster" )

    genome.gocount <- merge( genome.all.gocount , genome.cluster.gocount , by="GO" , all=TRUE )
    genome.gocount[ is.na( genome.gocount ) ] <- 0

    neg <- ddply( genome.gocount , "namespace" , function(x) sum(x$all) )
    
    genome.gocount$dhyp <- dhyper( genome.gocount$cluster , genome.gocount$all , neg[genome.gocount$namespace , 2] - genome.gocount$cluster , genome.cluster.size ) 

    genome.gocount$adj.p <- p.adjust( genome.gocount$dhyp , adjust.method )
    genome.gocount <- genome.gocount[ order(genome.gocount$adj.p) , ]

    return( genome.gocount )
}


####
## Returns all genes with a given GO ID
##
## goid: vector of GOIDs
##
go.getgenesforgoid <- function( go , goid ) {
    return( go[ which( go$GO %in% goid ) , ]$ORF )
}


####
## Create a word cloud scaled by -log10( adjusted probability )
##
## genome.gocount:      go enrichment analysis return by go.gethyper()
##
## p.cutfoff:           don't display words with adjusted probabilty below this cutoff
##
go.wordcloud <- function( genome.gocount , p.cutoff=0.01, filename ) {
    require( wordcloud )

    ## Remove spaces in between GO descriiptions
    text <- gsub( " " , "_" , genome.gocount$name )

    ## Do the actual plotting
    pdf( paste0("go_", filename, ".pdf") , 20 , 20 )
    wordcloud( text , ceiling(-log10(genome.gocount$adj.p ) ) , min.freq=-log10( p.cutoff ) )
    dev.off()
}















## source( "./go.r" )
go <- go.loadgo()

list_tsv <- system(command="ls *_preprocessed.tsv", intern=T)


for(i in 1:length(list_tsv)) {

filename <- gsub("_preprocessed.tsv", "", list_tsv[i])

message("Processing ", i, ": ", filename, "...")

raw <- fread(list_tsv[i], dec=",")

raw$log2fc <- as.numeric(raw$log2fc)

up <- unlist(raw[log2fc > 2, ]$Gene)
down <- unlist(raw[log2fc < -2, ]$Gene)

bg <- unlist(raw$Gene)

up.go <- go.gethyper( go , bg , up , "biological_process" , "fdr" )
go.wordcloud( up.go , p.cutoff=0.01, paste0("up_", filename) )
write.table( up.go , file=paste0("go_up_",  filename, ".tsv"), col.name=T, row.name=F, sep="\t", quote=F )

down.go <- go.gethyper( go , bg , down , "biological_process" , "fdr" )
go.wordcloud( down.go , p.cutoff=0.01, paste0("down_", filename) )
write.table( down.go , file=paste0("go_down_",  filename, ".tsv"), col.name=T, row.name=F, sep="\t", quote=F )


message("Done.")

}

## up <- read.csv( "up.csv" , stringsAsFactor=F )
## up <- unlist( up , use.names=F )
## bg <- read.csv( "urm1_bg.csv" , colClasses = "character" , stringsAsFactors=F )
## bg <- unlist( bg , use.names=F )
## up.go <- go.gethyper( go , bg , up , "biological_process" , "fdr" )
## go.wordcloud( up.go , p.cutoff=0.01 )
## write.csv( up.go , quote=T , "up_go.csv" )
