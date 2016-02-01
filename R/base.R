.structureAMPoutput<-function(data)
{  
  strand=rep("+",nrow(data))
  strand[data[,4]=="0"]="-"
  numCs=round(data[,5]*data[,6]/100)
  numTs=round(data[,5]*data[,7]/100)
  data.frame(chr=data[,2],start=data[,3],end=data[,3],strand=strand,coverage=data[,5],numCs=numCs,numTs=numTs)
}

valid.methylRawObj <- function(object) {
    data=getData(object)
    check1=( (object@resolution == "base") | (object@resolution == "region") )
    check2=(ncol(data)==7)
    if(check1 & check2){
      return(TRUE)
    }
    else if (! check1 ){
        cat("resolution slot has to be either 'base' or 'region':",
              "other values not allowed")
    }
    else if(! check2){
        cat("data part of methylRaw have",ncol(data),"columns, expected 7 columns")
    }
}

setClass("methylRawList", representation(treatment = "numeric"),contains = "list")

setClass("methylRaw", contains= "data.frame",representation(sample.id = "character", assembly = "character",context="character",resolution="character"),validity=valid.methylRawObj)

setGeneric("getData", def=function(x) standardGeneric("getData"))
setMethod("getData", signature="methylRaw", definition=function(x) {
                #return(as(x,"data.frame"))
                return(S3Part(x, strictS3 = TRUE))
})

setAs("methylRaw", "GRanges", function(from)
                      {
                        from2=getData(from)
                        GRanges(seqnames=from2$chr,ranges=IRanges(start=from2$start, end=from2$end),
                                       strand=from2$strand, 
                                       coverage=from2$coverage,
                                       numCs   =from2$numCs,
                                       numTs  =from2$numTs                                
                                       )
})
setGeneric("read", function(location,sample.id,assembly,pipeline="bismark",header=F,skip=0,sep="\t",context="CpG",resolution="base",treatment) standardGeneric("read"))

setMethod("read", signature(location = "list",sample.id="list",assembly="character"),function(location,sample.id,assembly,pipeline,header,skip,sep,context,resolution,treatment){ 
            outList=list()
            for(i in 1:length(location))
            {
              data<- read.big.matrix(location[[i]],header=F,sep="\t",type="integer")
              data<- .structureAMPoutput(data)
              obj=new("methylRaw",data,sample.id=sample.id[[i]],assembly=assembly,context=context,resolution=resolution)
              outList[[i]]=obj 
              rm(data)
              #print(location[[i]])      
            }
            myobj=new("methylRawList",outList,treatment=treatment)
            myobj
          })

setClass("methylBase",contains="data.frame",representation(
  sample.ids = "character", assembly = "character",context = "character",treatment="numeric",coverage.index="numeric",
                                   numCs.index="numeric",numTs.index="numeric",destranded="logical",resolution = "character"))

setAs("methylBase", "GRanges", function(from)
                      {
                        from=getData(from)
                        GRanges(seqnames=from$chr,ranges=IRanges(start=from$start, end=from$end),
                                       strand=from$strand, 
                                       data.frame(from[,5:ncol(from)])
                                       )

})
setGeneric("getContext", def=function(x) standardGeneric("getContext"))

#' @rdname getContext-methods
#' @aliases getContext,methylBase-method
setMethod("getContext", signature="methylBase", definition=function(x) {
                return(x@context)
        })

#' @rdname getContext-methods
#' @aliases getContext,methylRaw-method
setMethod("getContext", signature="methylRaw", definition=function(x) {
                return(x@context)
        })

setGeneric("getAssembly", def=function(x) standardGeneric("getAssembly"))

#' @rdname getAssembly-methods
#' @aliases getAssembly,methylBase-method
setMethod("getAssembly", signature="methylBase", definition=function(x) {
                return(x@assembly)
        }) 

#' @rdname getAssembly-methods
#' @aliases getAssembly,methylRaw-method
setMethod("getAssembly", signature="methylRaw", definition=function(x) {
                return(x@assembly)
        })

setGeneric("unite", function(object,destrand=FALSE,min.per.group=NULL) standardGeneric("unite"))
setMethod("unite", "methylRawList",function(object,destrand,min.per.group){
            if( length(unique(vapply(object,function(x) x@context,FUN.VALUE="character"))) > 1)
            {
              stop("supplied methylRawList object have different methylation contexts:not all methylation events from the same bases")
            }
            if( length(unique(vapply(object,function(x) x@assembly,FUN.VALUE="character"))) > 1)
            {
              stop("supplied methylRawList object have different genome assemblies")
            }                     
            if( length(unique(vapply(object,function(x) x@resolution,FUN.VALUE="character"))) > 1)
            {
              stop("supplied methylRawList object have different methylation resolutions:some base-pair some regional")
            } 
            
            if( (!is.null(min.per.group)) &  ( ! is.integer( min.per.group ) )  ){stop("min.per.group should be an integer\ntry providing integers as 1L, 2L,3L etc.\n")}
            
            #merge raw methylation calls together
            df=getData(object[[1]])
            if(destrand & (object[[1]]@resolution == "base") ){df=.CpG.dinuc.unify(df)}
            df=data.table(df,key=c("chr","start","end","strand"))
            sample.ids=c(object[[1]]@sample.id)
            assemblies=c(object[[1]]@assembly)
            contexts  =c(object[[1]]@context)
            for(i in 2:length(object))
            {
              df2=getData(object[[i]])
              if(destrand & (object[[1]]@resolution == "base") ){df2=.CpG.dinuc.unify(df2)}
              #
              
              if( is.null(min.per.group) ){
                df2=data.table(df2[,c(1:3,5:7)])
                df=merge(df,df2,by=c("chr","start","end"),suffixes=c(as.character(i-1),as.character(i) ) ) # merge the dat to a data.frame
              }else{
                df2=data.table(df2,key=c("chr","start","end","strand") )
                # using hacked data.table merge called merge2: temporary fix
                df=merge2(df,df2,by=c("chr","start","end","strand"),suffixes=c(as.character(i-1),as.character(i) ) ,all=TRUE)
              }
              sample.ids=c(sample.ids,object[[i]]@sample.id)
              contexts=c(contexts,object[[i]]@context)
            }
            
            # stop if the assembly of object don't match
            if( length( unique(assemblies) ) != 1 ){stop("assemblies of methylrawList elements should be same\n")}
            
            
            if(  ! is.null(min.per.group) ){
              # if the the min.per.group argument is supplied, remove the rows that doesn't have enough coverage
              
              # get indices of coverage,numCs and numTs in the data frame 
              coverage.ind=seq(5,by=3,length.out=length(object))
              numCs.ind   =coverage.ind+1
              numTs.ind   =coverage.ind+2
              start.ind   =2 # will be needed to weed out NA values on chr/start/end/strand
              
              for(i in unique(object@treatment) ){
                my.ind=coverage.ind[object@treatment==i]
                ldat = !is.na(df[,my.ind,with=FALSE])
                if(  is.null(dim(ldat))  ){  # if there is only one dimension
                  df=df[ldat>=min.per.group,]
                }else{
                  df=df[rowSums(ldat)>=min.per.group,]
                }
              }
              #mat=df[,c(start.ind-1,start.ind,start.ind+1,start.ind+2)] # get all location columns, they are now duplicated with possible NA values
              #locs=t(apply(mat,1,function(x) unique(x[!is.na(x)]) ) ) # get location matrix
              #if(ncol(locs)==3){ # if the resolution is base
              #  df[,c(2:5)]=data.frame(chr=locs[,1],start=as.numeric(locs[,2]),end=as.numeric(locs[,2]),strand=locs[,3])
              #}else{   # if the resolution is region
              #  df[,c(2:5)]=data.frame(chr=locs[,1],start=as.numeric(locs[,2]),end=as.numeric(locs[,3]),strand=locs[,4])
              #}
              #start.ind   =seq(10,by=7,length.out=length(object)) # will be needed to weed out NA values on chr/start/end/strand
              
              #df=df[,-c(start.ind-1,start.ind,start.ind+1,start.ind+2)]
              #names(df)[2:5]=c("chr","start","end","strand")
            }
            df=as.data.frame(df)
            # get indices of coverage,numCs and numTs in the data frame 
            coverage.ind=seq(5,by=3,length.out=length(object))
            numCs.ind   =coverage.ind+1
            numTs.ind   =coverage.ind+2
            
            # change column names
            names(df)[coverage.ind]=paste(c("coverage"),1:length(object),sep="" )
            names(df)[numCs.ind]   =paste(c("numCs"),1:length(object),sep="" )
            names(df)[numTs.ind]   =paste(c("numTs"),1:length(object),sep="" )
            
            obj=new("methylBase",(df),sample.ids=sample.ids,
                    assembly=unique(assemblies),context=unique(contexts),
                    treatment=object@treatment,coverage.index=coverage.ind,
                    numCs.index=numCs.ind,numTs.index=numTs.ind,destranded=destrand,resolution=object[[1]]@resolution )
            obj
          }
)           

setGeneric("select", def=function(x,i) standardGeneric("select"))




#' @aliases select,methylBase-method
#' @rdname select-methods
setMethod("select", "methylBase",
          function(x, i)
          {

            new("methylBase",getData(x)[i,],
                sample.ids = x@sample.ids, 
                assembly = x@assembly,
                context = x@context,
                treatment=x@treatment,
                coverage.index=x@coverage.index,
                numCs.index=x@numCs.index,
                numTs.index=x@numTs.index,
                destranded=x@destranded,
                resolution =x@resolution)
           }
)



#' @aliases select,methylRaw-method
#' @rdname select-methods
setMethod("select", "methylRaw",
          function(x, i)
          {

          new("methylRaw",getData(x)[i,],sample.id=x@sample.id,
                                           assembly=x@assembly,
                                           context=x@context,
                                           resolution=x@resolution)
           }
          

)
