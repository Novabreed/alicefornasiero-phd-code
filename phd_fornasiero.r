clean.geno<-function(infile="/projects/novabreed/share/afornasiero/RAD_seq/project_data/primitivo/43samples_300-450_fusion/cleaned_files/nohomo_noReAS_300-450_genomic_sorted.loc",
					outfile="/projects/novabreed/share/afornasiero/RAD_seq/project_data/primitivo/43samples_300-450_fusion/cleaned_files/clean_nohomo_noReAS_300-450_genomic_sorted.loc",
					max.missing.snps=0.2,max.missing.subj=0.2)

{
library(data.table)
raw.data<-fread(infile)
#Number of genotypes is the number of entry of the row table minus three (i.e. chrom, pos and hk x hk)
n.geno<-lapply(lapply(apply(raw.data,1,table,exclude="--"),length),"-",3)

#Remove lines of table with geno==1: these are monomorphic positions: they do not deserve to live!
raw.data<-raw.data[n.geno>1,]

#For each SNP, build a table of genotypes (removing unknowns)
real.geno<-apply(raw.data,2,table,exclude="--")
snps<-nrow(raw.data)
real.geno<-lapply(real.geno,sum)

#Count proportion of missing snps and of missing subjects. 
missing.snps<-1-(unlist(real.geno)/snps)
subjects<-ncol(raw.data[,4:ncol(raw.data),with=FALSE])
count.subjects<-apply(raw.data[,4:ncol(raw.data),with=FALSE],1,table,exclude="--")
count.subjects<-lapply(count.subjects,sum)
missing.subjects<-1-(unlist(count.subjects)/subjects)

#Subjects (columns) with too many missing data are overwritten with --
raw.data[,which(missing.snps>max.missing.snps)]<-rep("--",nrow(raw.data))

#SNPS (rows) supported by enough subjects are retained
raw.data<-raw.data[missing.subjects<=max.missing.subj,]
write.table(raw.data,outfile,sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

}





phase.geno<-function(infile="/projects/novabreed/share/afornasiero/RAD_seq/project_data/primitivo/43samples_300-450_fusion/cleaned_files/clean_nohomo_noReAS_300-450_genomic_sorted.loc",
					outfile="/projects/novabreed/share/afornasiero/RAD_seq/project_data/primitivo/43samples_300-450_fusion/cleaned_files/phased_nohomo_noReAS_300-450_genomic_sorted_0.2_0.2.loc",
					threshold.phase=0.2,threshold.homo=0.1,max.dist.snp=50)

{
library(data.table)
raw.data<-fread(infile,colClasses="character")
toswitch<-c(0,rep(NA,(nrow(raw.data)-1)))
#count subjects in the file
num.subjects<-ncol(raw.data)-3
#For each SNP, build a table of genotypes (removing unknowns)
real.geno<-apply(raw.data,2,table,exclude="--")
for(aaa in 2:(nrow(raw.data)-1))
{
#Initialize a variable. This will become 1 when we are able to assign the phase to the corresponding SNP pair
yes12<-yes23<-yes13<-0

#Manage the change of chromosomes.
#At the beginning of each chromosome we arbitrarily set the switch to 0 (i.e. the first SNP rules)
if(raw.data$V1[aaa-1]!=raw.data$V1[aaa] | raw.data$V1[aaa+1]!=raw.data$V1[aaa])
{
toswitch[aaa]<-0
next
}
#Manage situations emerging when the previous SNP didn't allow to attribute phase
#if the previous line was an NA we have to decide the switching based on the first SNP (going backward) for which we confidently assigned phase 
#(for our simplicity, we search for SNPs in cis with respect to the first SNP)

if(is.na(toswitch[aaa-1])) 
{
#List of all SNPS in cis with the first, from last to first
#We go back on the vector, searching for a SNP which besides being in cis with the first, showed both possible homozygous states. 
list.snp.phased<-sort(which(toswitch==0),decreasing=TRUE)
for(bbb in 1:length(list.snp.phased))
{
#Count how many names of the table of the current line are present in the vector ("hh","kk")
# aa is the table of gentoypes of line bbb
aa<-table(unlist(raw.data[list.snp.phased[bbb],4:ncol(raw.data),with=FALSE]))
# tt counts the occurrencies of homozygous genotypes in aa
tt<-sum(names(aa)%in%c("hh","kk"))
#Only when two homozygous are observed we select the SNP as a base for deciding who is in phase (it's not right to use SNPs with only one homozygous)
if(tt>1) {
last.snp.phased<-list.snp.phased[bbb]
#Two tricks to avoid lack of decisions!
if((aaa-last.snp.phased)>max.dist.snp) toswitch[aaa]<-0
if(raw.data$V1[last.snp.phased]!=raw.data$V1[aaa]) toswitch[aaa]<-0
cat("phased!\n")
break}
}
#So now the table is built between the current SNP (aaa) and the last snp for which we obtained reliable phase (last.snp.phased)
pino<-table(unlist(raw.data[aaa-1,4:ncol(raw.data),with=FALSE]),unlist(raw.data[aaa,4:ncol(raw.data),with=FALSE]),unlist(raw.data[aaa+1,4:ncol(raw.data),with=FALSE]),exclude="--")

} else pino<-table(unlist(raw.data[aaa-1,4:ncol(raw.data),with=FALSE]),unlist(raw.data[aaa,4:ncol(raw.data),with=FALSE]),unlist(raw.data[aaa+1,4:ncol(raw.data),with=FALSE]),exclude="--")


#Count cis instances.
#We work on three loci.
#cis13 means that the first SNP is in cis with the third
#trans13 means that the first SNP is in trans with the third
#We then count the occurrencies of each cis trans combo
cis12<-sum(pino[dimnames(pino)[[1]]=="hh",dimnames(pino)[[2]]=="hh",])+sum(pino[dimnames(pino)[[1]]=="kk",dimnames(pino)[[2]]=="kk",])
trans12<-sum(pino[dimnames(pino)[[1]]=="hh",dimnames(pino)[[2]]=="kk",])+sum(pino[dimnames(pino)[[1]]=="kk",dimnames(pino)[[2]]=="hh",])
cis13<-sum(pino[dimnames(pino)[[1]]=="hh",,dimnames(pino)[[3]]=="hh"])+sum(pino[dimnames(pino)[[1]]=="kk",,dimnames(pino)[[3]]=="kk"])
trans13<-sum(pino[dimnames(pino)[[1]]=="hh",,dimnames(pino)[[3]]=="kk"])+sum(pino[dimnames(pino)[[1]]=="kk",,dimnames(pino)[[3]]=="hh"])
cis23<-sum(pino[,dimnames(pino)[[2]]=="hh",dimnames(pino)[[3]]=="hh"])+sum(pino[,dimnames(pino)[[2]]=="kk",dimnames(pino)[[3]]=="kk"])
trans23<-sum(pino[,dimnames(pino)[[2]]=="hh",dimnames(pino)[[3]]=="kk"])+sum(pino[,dimnames(pino)[[2]]=="kk",dimnames(pino)[[3]]=="hh"])

#If ALL the three tables have a lack of homozygous we skip to the next SNP
if((cis12+trans12)/sum(pino)<threshold.homo & (cis23+trans23)/sum(pino)<threshold.homo & (cis13+trans13)/sum(pino)<threshold.homo) next

#Check if we can attribute phase to pair 1-2
if(cis12<trans12 & (cis12/trans12)<threshold.phase)
{
cat("Trans 1-2\n")
toswitch[aaa]<-1
yes12<-1
}
if(cis12>trans12 & (trans12/cis12)<threshold.phase)
{
cat("Cis 1-2\n")
toswitch[aaa]<-0
yes12<-1
}
#Check if we can attribute phase to pair 1-3
if(cis13<trans13 & (cis13/trans13)<threshold.phase)
{
cat("Trans 1-3\n")
toswitch[aaa+1]<-1
yes13<-1
}
if(cis13>trans13 & (trans13/cis13)<threshold.phase) 
{
cat("Cis 1-3\n")
toswitch[aaa+1]<-0
yes13<-1
}
#Check if we can attribute phase to pair 2-3
if(cis23<trans23 & (cis23/trans23)<threshold.phase)
{
cat("Trans 2-3\n")
toswitch[aaa+1]<-1
yes23<-1
}
if(cis23>trans23 & (trans23/cis23)<threshold.phase) 
{
cat("Cis 2-3\n")
toswitch[aaa+1]<-0
yes23<-1
}
#If we settled the pair 1-2 and decided to change the phase, then we change the phase and skip to the next in loop
if(yes12==1)
{
if(toswitch[aaa]==1)
{
raw.data[aaa,][raw.data[aaa,]=="hh"]<-"ff"
raw.data[aaa,][raw.data[aaa,]=="kk"]<-"gg"
raw.data[aaa,][raw.data[aaa,]=="gg"]<-"hh"
raw.data[aaa,][raw.data[aaa,]=="ff"]<-"kk"
#Change toswitch[aaa] to 0, because now we are in the same phase as the starting SNP
toswitch[aaa]<-0
next
}
}

#Otherwise, we settle (if possible) pairs 1-3 and 2-3  
if(yes13==1|yes23==1)
{
aaa<-aaa+1
if(toswitch[aaa]==1)
{
raw.data[aaa,][raw.data[aaa,]=="hh"]<-"ff"
raw.data[aaa,][raw.data[aaa,]=="kk"]<-"gg"
raw.data[aaa,][raw.data[aaa,]=="gg"]<-"hh"
raw.data[aaa,][raw.data[aaa,]=="ff"]<-"kk"
#Change toswitch[aaa] to 0, because now we are in the same phase as the starting SNP
toswitch[aaa-1]<-0
toswitch[aaa]<-0
}
}

#The part below will be executed only if it is not possible to decide if the two SNPs are in CIS or TRANS, i.e. if toswitch[aaa]==NA
cat(aaa,"\n")
}
write.table(raw.data,outfile,sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
}




define.haplo<-function(infile="/projects/novabreed/SNP/gatk/rkatsiteli/rkatsiteli_349_12b_19a_35b_58b/unifgenotyper/FC0637_FC0688/clean/SNP_del_ins/20160927_clean_hetSNP_DEL_INS_rkatsiteliALL_geno_only.goodReg.txt",
						samplenames=c("rkatsiteli","Rkatsiteli_349-P4-12B","Rkatsiteli_349-P4-19A","Rkatsiteli_349-P4-35B","Rkatsiteli_349-P4-58B"),
						outfile="/projects/novabreed/SNP/gatk/rkatsiteli/rkatsiteli_349_12b_19a_35b_58b/unifgenotyper/FC0637_FC0688/clean/SNP_del_ins/20160927_flip_clean_hetSNP_DEL_INS_rkatsiteliALL_geno_only.goodReg.txt",
						min.info=1,remove.het=TRUE,plot=TRUE,mystep=999,tol=0.00000001)

{
library(data.table)
geno<-fread(infile,data.table=FALSE)

#Define the most frequent pattern along each chromosome
#Remove triallelic SNPs.
triallele<-grep(",",geno$ALT,perl=TRUE)
geno<-geno[-triallele,]

#Prepare content of haplotypes
geno$hapA<-geno$REF
geno$hapB<-geno$ALT
samplenames<-unlist(samplenames)

#set remaining heterozygous genotypes (i.e. 1) in homozygous regions to unknown
if(remove.het) geno[,samplenames[2:length(samplenames)]][geno[,samplenames[2:length(samplenames)]]=="1"]<-"."

#Remove positions in which the parent is not heterozygous (0=homo_ref, 1=het, 2=homo_alt)
geno<-geno[geno[,samplenames[1]]=="1",]

# for each position, count how many times genotype is informative
cat("Computing informative positions\n")
geno$info<-apply(apply(geno[,samplenames[2:length(samplenames)]],1,"!=","."),2,sum,na.rm=T)

myseq<-seq(1,nrow(geno),by=mystep)
for(aaa in 1:(length(myseq)-1)) 
{
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
		mypossiblehap<-expand.grid(c(0,2),c(0,2),c(0,2),c(0,2))
		mypossiblehap$lik<-1
	}

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

	if(length(samplenames[2:length(samplenames)])==5)
	{
		mypossiblehap$lik[bbb]<-mypossiblehap$lik[bbb]*onefivetable[row.names(onefivetable)==mypossiblehap$Var1[bbb],colnames(onefivetable)==mypossiblehap$Var5[bbb]]
		mypossiblehap$lik[bbb]<-mypossiblehap$lik[bbb]*twofivetable[row.names(twofivetable)==mypossiblehap$Var2[bbb],colnames(twofivetable)==mypossiblehap$Var5[bbb]]
		mypossiblehap$lik[bbb]<-mypossiblehap$lik[bbb]*threefivetable[row.names(threefivetable)==mypossiblehap$Var3[bbb],colnames(threefivetable)==mypossiblehap$Var5[bbb]]
		mypossiblehap$lik[bbb]<-mypossiblehap$lik[bbb]*fourfivetable[row.names(fourfivetable)==mypossiblehap$Var4[bbb],colnames(fourfivetable)==mypossiblehap$Var5[bbb]]
	}

	#List the most probable configuration in this small piece of chromosome
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
	if(howmanyzero>howmanytwo) sgeno$likeref<-1 else sgeno$likeref<-0

	# flipfreq=1 means to flip haplotypes
	sgeno$flipfreq<-1-sgeno$agree/sgeno$info
	if(aaa==1) fullgeno<-sgeno else fullgeno<-rbind(fullgeno,sgeno)
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

if(plot)
{
	cat("Changing POS in numeric\n")
	geno$POS<-as.numeric(as.character(geno$POS))

	# we plot only those positions that have been established previously
	geno<-geno[geno$hapA!="N",]

	# we plot only one haplotype (e.g. hapB), which is 1 when equal to ref
	geno$toplot<-rep(0,nrow(geno))
	geno$toplot[geno$hapA==geno$REF]<-0
	geno$toplot[geno$hapB==geno$REF]<-1

	chr.name=paste("chr",seq(1,19),sep="")

	pdf.file=gsub(".txt",".pdf",outfile)
	cat("pdf output file:",pdf.file, "\n")
	pdf(pdf.file,width=30)

		for (ddd in 1:length(chr.name))
		{
		# define data.frame of the chr we want to plot
		small<-geno[geno$Chr==chr.name[ddd],]
		
		# plot one haplotype
		mysumm<-summarize.by(small[,c("POS","toplot")],step=100)
		plot(small$POS/1000000, small$toplot, xaxp=c(x1=0,x2=round(max(small$POS/1000000),0),n=round(max(small$POS/1000000),0)), ylim=c(0,1), main=chr.name[ddd], xlab="Position[Mbp]", type="l", lwd=0.5)
		plot(mysumm$POS/1000000, mysumm$toplot, xaxp=c(x1=0,x2=round(max(mysumm$POS/1000000),0),n=round(max(mysumm$POS/1000000),0)), col="blue", pch=20, main=chr.name[ddd], xlab="Position[Mbp]", ylab="Haplotypes", cex.main=3, cex.lab=1.5, cex.axis=1.5)
		
		}
	dev.off()
}
}
