.formatDAVIDResult <- function(result, verbose=FALSE) {
	### we always use read.delim(...header=TRUE) but formatting expects the first row tobe the column names
	### in order to make formatting work we add the top row
	result<-rbind(colnames(result),result);

	tool <- attr(result,"tool") 
	if(verbose) 
		cat("formatDAVIDResult: tool=", tool, 
			ifelse(tool=="geneReportFull", 
			" (invisible return)", ""), 
		"\n")

	### formatting depends on which is done
	if(tool=="geneReportFull") {
		returnval <- try(formatGeneReportFull(result))
	} else if(tool=="geneReport") {
		returnval <- try(formatGeneReport(result))
	} else if(tool=="list") {
		returnval <- try(formatList(result))
	} else if(tool=="gene2gene") {
		returnval <- try(formatGene2Gene(result))
	} else if(tool=="annotationReport") {
		returnval <- try(formatAnnotationReport(result))
	} else
		returnval <- result  ### Unformatted for now.
	if(class(returnval) == "try-error")
		returnval <- result
	for(attname in c("annot", "ids", "tool", "type"))
		attr(returnval, attname) <- attr(result, attname)
	returnval
}

.bracketedStrings <- function(s, before, after, verbose=FALSE, addNames=FALSE, drop.na=TRUE, warn.if.gt.1=TRUE) {
	if(length(s) > 1) {
		result <- lapply(s, .bracketedStrings, 
			before=before, after=after, verbose=FALSE)
		if(addNames) names(result) = s
		return(result)
	}
	starts <- (valStrings <- gregexpr(before, s)[[1]]) + attr(valStrings, "match.length")
	ends <- regexpr(after, (vStrings <- substring(s, starts,1e6)) ) - 2
	result <- substring(s, starts, starts+ends)
	if(verbose)
		cat(paste("=>",starts, starts+ends, result, sep=":", collapse="\n"))
	result <- result[result != ""]
	if(drop.na) result <- result[!is.na(result)]
	result <- result[starts>=0 & ends>=0]
	if((length(result) > 1) & (warn.if.gt.1))
		warning("More than one substring found.")
	return(result)
}

## A litle hack to DAVIDQuery as the DAVIDQuery got deprected
.DAVIDQuery<-function (ids = "O00161,O75396", type = "UNIPROT_ACCESSION", 
    annot, tool, URLlengthLimit = 2048, details = TRUE, verbose = FALSE, 
    writeHTML = FALSE, testMe = FALSE, graphicMenu = FALSE, formatIt = TRUE) 
{
	DAVIDURLBase <- "https://david.abcc.ncifcrf.gov/"
    ids <- paste(ids, collapse = ",")
    ids <- paste(strsplit(ids, " ")[[1]], sep = "", collapse = "")
    firstURLOK <- FALSE
    while (firstURLOK == FALSE) {
        firstURL <- paste(DAVIDURLBase, "api.jsp?", "type=", 
            type, "&ids=", ids, "&tool=", tool, sep = "")
        if (!is.null(annot)) 
            firstURL <- paste(firstURL, "&annot=", annot, sep = "")
        if (verbose) 
            cat("DAVIDQuery:  firstURL = ", firstURL, "\n")
        if (nchar(firstURL) < URLlengthLimit) 
            firstURLOK <- TRUE
        else ids <- ids[-length(ids)]
    }
    DAVIDQueryResult <- try({
        myCurlHandle <- RCurl::getCurlHandle(cookiefile = "DAVIDCookiefile.txt")
        firstStageResult <- RCurl::getURL(firstURL, curl = myCurlHandle, 
            verbose = FALSE,ssl.verifypeer=FALSE)
        if (writeHTML) 
            writeChar(firstStageResult, "firstStageResult.html")
        DAVIDaction <- .bracketedStrings(firstStageResult, "document.apiForm.action = \"", 
            "\"")
        DAVIDvalues <- .bracketedStrings(firstStageResult, "document.apiForm.[a-z]*.value=\"", 
            "\"", warn.if.gt.1 = FALSE)
        if(DAVIDvalues[1] != ""){
           tmp <- unlist(strsplit(DAVIDvalues[1],split=","))
           if(length(tmp >=40)){
             DAVIDvalues[1] <- paste(tmp[1:40],collapse=",")
           }
        }
        DAVIDfields <- .bracketedStrings(firstStageResult, "document.apiForm.", 
            ".value=\"", warn.if.gt.1 = FALSE)
        secondURL <- paste(DAVIDURLBase, DAVIDaction, "?", 
            paste(DAVIDfields, "=", DAVIDvalues, sep = "", collapse = "&"), 
            sep = "")
        if (verbose) 
            cat("DAVIDQuery:  secondURL = ", secondURL, "\n")
        if (nchar(secondURL) > URLlengthLimit) 
            stop(paste("nchar(secondURL) too long; ", nchar(secondURL), 
                ">", URLlengthLimit))
        secondStageResult <- RCurl::getURL(secondURL, curl = myCurlHandle, 
            verbose = FALSE, ssl.verifypeer=FALSE)
        hasSessionEnded <- length(grep("Your session has ended", 
            secondStageResult) > 0)
        if (hasSessionEnded) 
            warning("Warning: Session ended")
        if (writeHTML) 
            writeChar(secondStageResult, "secondStageResult.html")
        downloadFileName <- .bracketedStrings(secondStageResult, 
            "href=\"data/download/", "\" target=")
        if (length(downloadFileName) == 0) 
            warning("Warning: downloadFileName is not found in reply html. \n")
        downloadURL <- paste(DAVIDURLBase, "data/download/", 
            downloadFileName, sep = "")
        if (verbose) 
            cat("downloadURL = ", downloadURL, "\n")

	  if (tool=="geneReport"){
	  	# work around the format in which the file for 'geneReport' is returned by DAVID 
		read.delim(downloadURL,stringsAsFactors=FALSE,header=TRUE,nrows=0);
        } else {
		read.delim(downloadURL, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE);
	  }
    })
    try(if (is.data.frame(DAVIDQueryResult) & (length(DAVIDQueryResult) > 
        0)) {
        if (length(grep("<title>Directory Listing For /data/download/</title>", 
            DAVIDQueryResult[[1]])) > 0) {
            DAVIDQueryResult <- paste("No result file was found. URL = ", 
                DAVIDQueryResult$firstURL)
            class(DAVIDQueryResult) <- "try-error"
        }
    })
    attr(DAVIDQueryResult, "ids") <- ids
    attr(DAVIDQueryResult, "tool") <- tool
    attr(DAVIDQueryResult, "annot") <- annot
    attr(DAVIDQueryResult, "type") <- type
    if (formatIt & (class(DAVIDQueryResult) != "try-error")) {
        DAVIDQueryResult <- .formatDAVIDResult(DAVIDQueryResult)
    }
    if (details) 
        return(list(ids = ids, firstURL = firstURL, firstStageResult = firstStageResult, 
            DAVIDaction = DAVIDaction, secondURL = secondURL, 
            secondStageResult = secondStageResult, hasSessionEnded = hasSessionEnded, 
            downloadFileName = downloadFileName, downloadURL = downloadURL, 
            DAVIDQueryResult = DAVIDQueryResult))
    return(DAVIDQueryResult)
}

