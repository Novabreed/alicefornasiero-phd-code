define.haplo<-function(infile="/projects/novabreed/SNP/gatk/rkatsiteli/rkatsiteli_349_12b_19a_35b_58b/unifgenotyper/FC0637_FC0688/clean/SNP_del_ins/20160927_clean_hetSNP_DEL_INS_rkatsiteliALL_geno_only.goodReg.txt",
						samplenames=c("rkatsiteli","Rkatsiteli_349-P4-12B","Rkatsiteli_349-P4-19A","Rkatsiteli_349-P4-35B","Rkatsiteli_349-P4-58B"),
						outfile="/projects/novabreed/SNP/gatk/rkatsiteli/rkatsiteli_349_12b_19a_35b_58b/unifgenotyper/FC0637_FC0688/clean/SNP_del_ins/20160927_flip_clean_hetSNP_DEL_INS_rkatsiteliALL_geno_only.goodReg.txt",
						min.info=1,remove.het=TRUE,plot=TRUE,mystep=999,tol=0.00000001)

{
library(data.table)
geno<-fread(infile,data.table=FALSE)

# cat("Extracting subset for debugging\n")
# geno<-geno[geno$Chr=="chr1",]
# geno<-geno[1:10000,]
#Define the most frequent pattern along each chromosome
#Remove triallelic SNPs.
# geno<-geno[nchar(geno$ALT)==1,]	# OLD way. Not compatible with haplotypes containing INDEL
triallele<-grep(",",geno$ALT,perl=TRUE)
geno<-geno[-triallele,]
#setnames(geno,"#CHROM","Chr")
#Prepare content of haplotypes
geno$hapA<-geno$REF
geno$hapB<-geno$ALT

samplenames<-unlist(samplenames)
#set remaining heterozygous genotypes (i.e. 1) in homozygous regions to unknown
if(remove.het) geno[,samplenames[2:length(samplenames)]][geno[,samplenames[2:length(samplenames)]]=="1"]<-"."

#Remove positions in which the parent is not heterozygous (0=homo_ref, 1=het, 2=homo_alt)
#geno<-geno[geno[,samplenames[1]]=="1",]	#already done if you run get.only.clean.SNP

# for each position, count how many times genotype is informative
cat("Computing informative positions\n")
geno$info<-apply(apply(geno[,samplenames[2:length(samplenames)]],1,"!=","."),2,sum,na.rm=T)

myseq<-seq(1,nrow(geno),by=mystep)
# browser()
for(aaa in 1:(length(myseq)-1)) 
{
# if(aaa==1000) browser()
cat("subset",aaa,"\n")
sgeno<-geno[myseq[aaa]:(myseq[aaa+1]-1),]
#Check the most probable pattern between subjects one and two
onetwotable<-table(factor(sgeno[,5],level=c(0,2)),factor(sgeno[,6],level=c(0,2)),exclude=".")
onetwotable<-(onetwotable+0.0001)/(sum(onetwotable)+0.0001)

#Check the most probable pattern between subjects one and three
onethreetable<-table(factor(sgeno[,5],level=c(0,2)),factor(sgeno[,7],level=c(0,2)),exclude=".")
onethreetable<-(onethreetable+0.0001)/(sum(onethreetable)+0.0001)

#Check the most probable pattern between subjects one and three
onefourtable<-table(factor(sgeno[,5],level=c(0,2)),factor(sgeno[,8],level=c(0,2)),exclude=".")
onefourtable<-(onefourtable+0.0001)/(sum(onefourtable)+0.0001)

#Check the most probable pattern between subjects two and three
twothreetable<-table(factor(sgeno[,6],level=c(0,2)),factor(sgeno[,7],level=c(0,2)),exclude=".")
twothreetable<-(twothreetable+0.0001)/(sum(twothreetable)+0.0001)

#Check the most probable pattern between subjects two and three
twofourtable<-table(factor(sgeno[,6],level=c(0,2)),factor(sgeno[,8],level=c(0,2)),exclude=".")
twofourtable<-(twofourtable+0.0001)/(sum(twofourtable)+0.0001)

#Check the most probable pattern between subjects three and four
threefourtable<-table(factor(sgeno[,7],level=c(0,2)),factor(sgeno[,8],level=c(0,2)),exclude=".")
threefourtable<-(threefourtable+0.0001)/(sum(threefourtable)+0.0001)

################ EDITED BY ALICE ################
if(length(samplenames[2:length(samplenames)])==5)
{
#Check the most probable pattern between subjects one and five
onefivetable<-table(factor(sgeno[,5],level=c(0,2)),factor(sgeno[,9],level=c(0,2)),exclude=".")
onefivetable<-(onefivetable+0.0001)/(sum(onefivetable)+0.0001)
#Check the most probable pattern between subjects two and five
twofivetable<-table(factor(sgeno[,6],level=c(0,2)),factor(sgeno[,9],level=c(0,2)),exclude=".")
twofivetable<-(twofivetable+0.0001)/(sum(twofivetable)+0.0001)
#Check the most probable pattern between subjects three and five
threefivetable<-table(factor(sgeno[,7],level=c(0,2)),factor(sgeno[,9],level=c(0,2)),exclude=".")
threefivetable<-(threefivetable+0.0001)/(sum(threefivetable)+0.0001)
#Check the most probable pattern between subjects four and five
fourfivetable<-table(factor(sgeno[,8],level=c(0,2)),factor(sgeno[,9],level=c(0,2)),exclude=".")
fourfivetable<-(fourfivetable+0.0001)/(sum(fourfivetable)+0.0001)

mypossiblehap<-expand.grid(c(0,2),c(0,2),c(0,2),c(0,2),c(0,2))
mypossiblehap$lik<-1
} else {
###############################################
mypossiblehap<-expand.grid(c(0,2),c(0,2),c(0,2),c(0,2))
mypossiblehap$lik<-1
}
#Manca da fare 
#Sistemare il problema che quando la tabelle (onetwotable) ha una sola colonna o sola riga, non viene trovata alcuna corrispondenza con il valore presente (ad esempio in Var2). Bisognerebbe forzare la tabella ad avere sempre due righe e due colonne con i nomi che diciamo noi e i conteggi che vogliamo (forse invece di usare table ci possiamo calcolare noi i valori a forza di AND e OR)
#Dopo questo, bisogna far in modo da incorporare il flipfreq nella funzione vecchia!

for(bbb in 1:nrow(mypossiblehap))
{
#Compute the most likely pattern of haplotype in the chromosome
mypossiblehap$lik[bbb]<-mypossiblehap$lik[bbb]*onetwotable[row.names(onetwotable)==mypossiblehap$Var1[bbb],colnames(onetwotable)==mypossiblehap$Var2[bbb]]
mypossiblehap$lik[bbb]<-mypossiblehap$lik[bbb]*onethreetable[row.names(onethreetable)==mypossiblehap$Var1[bbb],colnames(onethreetable)==mypossiblehap$Var3[bbb]]
mypossiblehap$lik[bbb]<-mypossiblehap$lik[bbb]*onefourtable[row.names(onefourtable)==mypossiblehap$Var1[bbb],colnames(onefourtable)==mypossiblehap$Var4[bbb]]
mypossiblehap$lik[bbb]<-mypossiblehap$lik[bbb]*twothreetable[row.names(twothreetable)==mypossiblehap$Var2[bbb],colnames(twothreetable)==mypossiblehap$Var3[bbb]]
mypossiblehap$lik[bbb]<-mypossiblehap$lik[bbb]*twofourtable[row.names(twofourtable)==mypossiblehap$Var2[bbb],colnames(twofourtable)==mypossiblehap$Var4[bbb]]
mypossiblehap$lik[bbb]<-mypossiblehap$lik[bbb]*threefourtable[row.names(threefourtable)==mypossiblehap$Var3[bbb],colnames(threefourtable)==mypossiblehap$Var4[bbb]]
}
################ EDITED BY ALICE ################
if(length(samplenames[2:length(samplenames)])==5)
{
mypossiblehap$lik[bbb]<-mypossiblehap$lik[bbb]*onefivetable[row.names(onefivetable)==mypossiblehap$Var1[bbb],colnames(onefivetable)==mypossiblehap$Var5[bbb]]
mypossiblehap$lik[bbb]<-mypossiblehap$lik[bbb]*twofivetable[row.names(twofivetable)==mypossiblehap$Var2[bbb],colnames(twofivetable)==mypossiblehap$Var5[bbb]]
mypossiblehap$lik[bbb]<-mypossiblehap$lik[bbb]*threefivetable[row.names(threefivetable)==mypossiblehap$Var3[bbb],colnames(threefivetable)==mypossiblehap$Var5[bbb]]
mypossiblehap$lik[bbb]<-mypossiblehap$lik[bbb]*fourfivetable[row.names(fourfivetable)==mypossiblehap$Var4[bbb],colnames(fourfivetable)==mypossiblehap$Var5[bbb]]
}
###############################################

# Check if the best likelihood is too close to the second best (in that case we will not assign hpalotypes and put Ns instead)
# epsilon<-sort(mypossiblehap$lik,decreasing=TRUE)[1]-sort(mypossiblehap$lik,decreasing=TRUE)[2]
# cat("Loop",aaa,"epsilon=",epsilon,"\n")

# if(epsilon<=tol)
# {
# We arbitrarily set info=0 because we do not know the haplotype and we want flip freq to be NaN
# sgeno$info<-sgeno$agree<-0
# sgeno$flipfreq<-1-sgeno$agree/sgeno$info
# if(aaa==1) fullgeno<-sgeno else fullgeno<-rbind(fullgeno,sgeno)
# next
# }
#List the most probable configuration in this small piece of chromosome
# configuration<-mypossiblehap[,1:4][which.max(mypossiblehap$lik),]
configuration<-mypossiblehap[,1:length(grep("Var",(colnames(mypossiblehap))))][which.max(mypossiblehap$lik),]
sgeno$agree<-rep(0,nrow(sgeno))
cat("Defining most probable configuration in subset",aaa,"\n")
namesnofirst<-samplenames[2:length(samplenames)]
for(ccc in 1:length(namesnofirst))
{
sgeno$agree<-sgeno$agree+as.numeric(sgeno[,namesnofirst[ccc]]==as.character(configuration[ccc]))
}
howmanyzero<-sum((unlist(sgeno[,namesnofirst]))=="0")
howmanytwo<-sum((unlist(sgeno[,namesnofirst]))=="2")
# in case of a tie between 0 and 2
#if(howmanyzero==howmanytwo) sgeno$likeref<-which.max(table(fullgeno$likeref[myseq[aaa-1]:myseq[aaa]-1]))
if(howmanyzero>howmanytwo) sgeno$likeref<-1 else sgeno$likeref<-0

# flipfreq=1 means to flip haplotypes
sgeno$flipfreq<-1-sgeno$agree/sgeno$info
#if(aaa==54) browser()
if(aaa==1) fullgeno<-sgeno else fullgeno<-rbind(fullgeno,sgeno)
#if(aaa==100) browser()
}

geno<-fullgeno

#Change NAs to 0.5 because otherwise we have problems in using suffixes.
#Positions in which flipfreq is NAs or 0.5 will be moved to Ns 
geno$flipfreq[is.na(geno$flipfreq)]<-0.5

geno$hapA[geno$flipfreq>0.5&geno$likeref==0]<-geno$REF[geno$flipfreq>0.5&geno$likeref==0]
geno$hapB[geno$flipfreq>0.5&geno$likeref==0]<-geno$ALT[geno$flipfreq>0.5&geno$likeref==0]

geno$hapA[geno$flipfreq>0.5&geno$likeref==1]<-geno$ALT[geno$flipfreq>0.5&geno$likeref==1]
geno$hapB[geno$flipfreq>0.5&geno$likeref==1]<-geno$REF[geno$flipfreq>0.5&geno$likeref==1]

geno$hapA[geno$flipfreq<0.5&geno$likeref==0]<-geno$ALT[geno$flipfreq<0.5&geno$likeref==0]
geno$hapB[geno$flipfreq<0.5&geno$likeref==0]<-geno$REF[geno$flipfreq<0.5&geno$likeref==0]

geno$hapA[geno$flipfreq<0.5&geno$likeref==1]<-geno$REF[geno$flipfreq<0.5&geno$likeref==1]
geno$hapB[geno$flipfreq<0.5&geno$likeref==1]<-geno$ALT[geno$flipfreq<0.5&geno$likeref==1]

geno$hapA[geno$flipfreq==0.5]<-"N"
geno$hapB[geno$flipfreq==0.5]<-"N"

cat("output file is being written\n")
write.table(geno,outfile,sep="\t",row.names=FALSE,quote=FALSE)
# browser()

if(plot)
{
cat("Changing POS in numeric\n")
geno$POS<-as.numeric(as.character(geno$POS))

# we plot only those positions that have been established previously
geno<-geno[geno$hapA!="N",]


# we plot only one haplotype (e.g. hapB), which is 1 when equal to ref
# geno$toplot<-rep(0,nrow(geno))
# geno$toplot[geno$hapA==geno$REF]<-0
# geno$toplot[geno$hapB==geno$REF]<-1
######
# we plot both haplotypes
geno$toplotA<-rep(0,nrow(geno))
geno$toplotB<-rep(0,nrow(geno))
geno$toplotA[geno$hapA==geno$REF]<-0
geno$toplotA[geno$hapB==geno$REF]<-1
geno$toplotB[geno$hapA==geno$REF]<-1
geno$toplotB[geno$hapB==geno$REF]<-0
#######



chr.name=paste("chr",seq(1,19),sep="")

pdf.file=gsub(".txt",".pdf",outfile)
cat("pdf output file:",pdf.file, "\n")
pdf(pdf.file,width=30)

	for (ddd in 1:length(chr.name))
	{
	# define data.frame of the chr we want to plot
	small<-geno[geno$Chr==chr.name[ddd],]
	
	# plot one haplotype
	# mysumm<-summarize.by(small[,c("POS","toplot")],step=100)
	# plot(small$POS/1000000, small$toplot, xaxp=c(x1=0,x2=round(max(small$POS/1000000),0),n=round(max(small$POS/1000000),0)), ylim=c(0,1), main=chr.name[ddd], xlab="Position[Mbp]", type="l", lwd=0.5)
	# plot(mysumm$POS/1000000, mysumm$toplot, xaxp=c(x1=0,x2=round(max(mysumm$POS/1000000),0),n=round(max(mysumm$POS/1000000),0)), col="blue", pch=20, main=chr.name[ddd], xlab="Position[Mbp]", ylab="Haplotypes", cex.main=3, cex.lab=1.5, cex.axis=1.5)
	
	# plot both halotypes
	mysummA<-summarize.by(small[,c("POS","toplotA")],step=100)
	mysummB<-summarize.by(small[,c("POS","toplotB")],step=100)
	plot(mysummA$POS/1000000, mysummA$toplotA, xaxp=c(x1=0,x2=round(max(mysummA$POS/1000000),0),n=round(max(mysummA$POS/1000000),0)), yaxp=c(0,1,1), col="navyblue", pch=19, main=paste("Chromosome",ddd,sep=" "), xlab="Position[Mbp]", ylab="Haplotypes", cex.main=2, cex.lab=1.5, cex.axis=1.5)
	points(mysummB$POS/1000000, mysummB$toplotB, col="darkorange", pch=19)
	
	
	# add number of informative sample number at each point
	# points(small$POS/1000000,small$info/length(samplenames[2:length(samplenames)]),col="grey",pch=20,cex=0.6)
	}
dev.off()
}
}
