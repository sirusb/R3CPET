ConvertGeneNameToHGNC<-function(folder){
  
  files<-list.files(folder,pattern=".gml",full.names=TRUE);
  
  for(f in files){
    g<-read.graph(f,format="gml");
    genesNode <-which(V(g)$gene==1)
    genes<-V(g)$name[genesNode];
    if(length(genes)>0){
      res<-EntrezToHGNC(genes);
      pos<-match(genes, res$entrezgene);
      converted<-which(!is.na(pos));
      pos<-pos[!is.na(pos)];
      #Here in case it returns some empty names
      nonMapped <-which(res$hgnc_symbol == "")
      res$hgnc_symbol[nonMapped] <- res$entrezgene[nonMapped];
      #Write the converted result
      V(g)[ genesNode[converted] ]$name <- res$hgnc_symbol[pos]
      f<-gsub(".gml","HGNC.gml",f);
      write.graph(g,file=f,format="gml")
    }
  }
}

##  S3 function Convert gene names from Entrez to HGNC
EntrezToHGNC<-function(EntrezID){  
  requireNamespace("biomaRt")
  ensemble<-useMart("ensembl");
  hsp<-useDataset(mart=ensemble,dataset="hsapiens_gene_ensembl");
  ids<-getBM(filters= "entrezgene",
             attributes= c("entrezgene","hgnc_id", "hgnc_symbol","description"),
             values= EntrezID, mart= hsp);
  return(ids);
}

## S3 function to convert gene names from Ensemble to HGNC
EnsemblToHGNC<-function(EnsemblIDs){
  requireNamespace("biomaRt")
  ensemble<-useMart("ensembl");
  hsp<-useDataset(mart=ensemble,dataset="hsapiens_gene_ensembl");
  ids<-getBM(filters= "ensembl_gene_id",
             attributes= c("ensembl_gene_id","hgnc_id", "hgnc_symbol","description"),
             values= EnsemblIDs, mart= hsp);
  return(ids);
}