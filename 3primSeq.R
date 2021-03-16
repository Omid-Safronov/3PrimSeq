#First input file is the blast results with output format of 6
#qseqid	qlen	sseqid	slen	qstart	qend	sstart	send	evalue	bitscore	score	length	pident	nident	mismatch	positive	gapopen	gaps	ppos	qframe	qcovhsp	qseq	qseq_full
#NS500198:409:HKTHLAFXY:1:11101:19405:13594	60	ND4_loc_10760_12137	52	1	52	1	52	6.26e-25	97.1	52	52	100.000	52	0	52	00	100.00	1	87	CCCATTCTCCTCCTATCCCTCAACCCCGACATCATTACCGGGTTTTCCTCTT	cccattctcctcctatccctcaaccccgacatcattaccgggttttcctcttAAAAAAAA
# second input files required is the primer coordinate file
#ID	Stop_Seq	prim5_coordinate	Length_of_the_primer	GSP_seq	Extra_seq
#ND6_loc_complement14149_14673	ATAGG	63	26	GTTACTGGTTGAACATTGTTTGTTGG	CCCCTACTCCTCCTAGACCTAACCTGACTAGAAAAGCTATTACCTAAAACAATTTCACAGCACCAAATCTCCACCTCCATCATCACCTCAACCCAAAAAGGCATAATTAAACTTTACTTCCTCTCTTTCTTCTTCCCACTCATCCTAACCCTACTCCTAATCACATAA

BltFastaSub_Kah=list()
for (i in 1:length(BltFasta)){
  BltFastaSub_Kah[[i]]=BltFasta[[i]][as.vector(BltFasta[[i]]$sseqid) %in% as.vector(StopCodon$ID),]
}
names(BltFastaSub_Kah)=names(BltFasta)

for (i in 1:length(BltFastaSub_Kah)){
  print(i)
  T=sapply(BltFastaSub_Kah[[i]][,23], function(x) str_extract_all(as.vector(x),"A-Z][a-z]*|[a-z]+|[A-Z]+"))
  MatT=matrix(NA,length(T),1)
  for (j in 1:length(T)){
    #print(j)
    Scod=tolower(as.character(StopCodon[grep(as.vector(BltFastaSub_Kah[[i]]$sseqid[j]),StopCodon$ID),]$Stop_Seq))
    Len=StopCodon[grep(as.vector(BltFastaSub_Kah[[i]]$sseqid[j]),StopCodon$ID),][,3]
    if (length(T[[j]]) == 1 && grepl("^[a-z]+$", T[[j]][1], perl=T)){
      Qcod=paste(strsplit(T[[j]],"")[[1]][(length(strsplit(T[[j]],"")[[1]])-4):length(strsplit(T[[j]],"")[[1]])],collapse="")
      if (Qcod == Scod && nchar(T[[j]]) == Len){MatT[j,1]="TrueEnd"}
      else if (Qcod != Scod && nchar(T[[j]]) == Len){MatT[j,1]="notrueend"}
      else if (Qcod == Scod && nchar(T[[j]]) > Len){MatT[j,1]="Long_TrueEnd"}
      else if (Qcod != Scod && nchar(T[[j]]) > Len){MatT[j,1]="Long_notrueend"}
      else if (Qcod == Scod && nchar(T[[j]]) < Len){MatT[j,1]="Short_TrueEnd"}
      else if (Qcod != Scod && nchar(T[[j]]) < Len){MatT[j,1]="Short_notrueend"}
    }
    if (length(T[[j]]) == 2 && grepl("^[a-z]+$", T[[j]][2], perl=T)){
      Qcod=paste(strsplit(T[[j]][[2]],"")[[1]][(length(strsplit(T[[j]][[2]],"")[[1]])-4):length(strsplit(T[[j]][[2]],"")[[1]])],collapse="")
      if (Qcod == Scod && nchar(T[[j]][[2]]) == Len){MatT[j,1]="TrueEnd"}
      else if (Qcod != Scod && nchar(T[[j]][[2]]) == Len){MatT[j,1]="notrueend"}
      else if (Qcod == Scod && nchar(T[[j]][[2]]) > Len){MatT[j,1]="Long_TrueEnd"}
      else if (Qcod != Scod && nchar(T[[j]][[2]]) > Len){MatT[j,1]="Long_notrueend"}
      else if (Qcod == Scod && nchar(T[[j]][[2]]) < Len){MatT[j,1]="Short_TrueEnd"}
      else if (Qcod != Scod && nchar(T[[j]][[2]]) < Len){MatT[j,1]="Short_notrueend"}
    }
    if (length(T[[j]]) == 2 && grepl("^[a-z]+$", T[[j]][1], perl=T)){
      Target=as.character(T[[j]][[2]])
      Qcod=paste(strsplit(T[[j]][[1]],"")[[1]][(length(strsplit(T[[j]][[1]],"")[[1]])-4):length(strsplit(T[[j]][[1]],"")[[1]])],collapse="")
      xSeq=as.character(paste(strsplit(StopCodon[match(as.vector(BltFastaSub_Kah[[i]]$sseqid)[j],StopCodon$ID),]$Extra_seq,"")[[1]][1:nchar(Target)],collapse=""))
      if (Qcod == Scod && nchar(T[[j]][[1]]) == Len){
        if(grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) <= 1 && Target != "A"){MatT[j,1]="TrueEnd"}
        else if (!grepl("^[AA]+$", Target, perl=T) && pid(pairwiseAlignment(as.vector(Target), xSeq, gapOpening=5,gapExtension=5, scoreOnly = FALSE)) != "NaN" && pid(pairwiseAlignment(as.vector(Target), xSeq, gapOpening=5,gapExtension=5, scoreOnly = FALSE)) == 100){
          if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="TrueEnd_Extra"}
        }
        else if (str_detect(Target, "^A+$")  && nchar(Target) >= 1){MatT[j,1]="TrueEnd_OligoA"}
        else if (grepl("^[A-Z]+", Target, perl=T) && (str_detect(Target, "^AAA+") || str_detect(Target, "AAA+$")) && length(table(unlist(strsplit(Target,"")))) >= 2){
          if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) > 2 && ((str_detect(Target, "^AAAA+") && str_detect(Target, "AAA+$")) || (str_detect(Target, "^AAAA+") || str_detect(Target, "AAA+$")))  && length(names(table(strsplit(Target,"")[[1]]))) > 2){MatT[j,1]="TrueEnd_ExtraOligoA"}
          else if (grepl("^[A-Z]+", Target, perl=T) && nchar(Target) >=4 && str_detect(paste(strsplit(Target,"")[[1]][1:(length(strsplit(Target,"")[[1]])-3)],collapse=""),"AAAAAA+$")){MatT[j,1]="TrueEnd_OligoA"}
          else if (str_detect(Target, "AA+$") && nchar(Target) >= 2){MatT[j,1]="TrueEnd_OligoA"}
          else if ((grepl("^[AA]+", Target, perl=T) || grepl("[AA]+$", Target, perl=T)) && length(table(unlist(strsplit(Target,"")))) <= 2 && nchar(Target) >= 2){MatT[j,1]="TrueEnd_OligoA"}
          else if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && ((!str_detect(Target, "^AAAA+") && !str_detect(Target, "AAAA+$")) || !str_detect(Target, "^AAAA+") && !str_detect(Target, "AAAA+$")) && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="TrueEnd_Extra"}
        }
        else if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="TrueEnd_Extra"}
      }
      if (Qcod != Scod && nchar(T[[j]][[1]]) == Len){
        if(grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) <= 1 && Target != "A"){MatT[j,1]="notrueend"}
        else if (!grepl("^[AA]+$", Target, perl=T) && pid(pairwiseAlignment(as.vector(Target), xSeq, gapOpening=5,gapExtension=5, scoreOnly = FALSE)) != "NaN" && pid(pairwiseAlignment(as.vector(Target), xSeq, gapOpening=5,gapExtension=5, scoreOnly = FALSE)) == 100){
          if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="notrueend_Extra"}
        }
        else if (str_detect(Target, "^A+$") && nchar(Target) >= 1){MatT[j,1]="notrueend_OligoA"}
        else if (grepl("^[A-Z]+", Target, perl=T) && (str_detect(Target, "^AAA+") || str_detect(Target, "AAA+$")) && length(table(unlist(strsplit(Target,"")))) >= 2){
          if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) > 2 && ((str_detect(Target, "^AAAA+") && str_detect(Target, "AAA+$")) || (str_detect(Target, "^AAAA+") || str_detect(Target, "AAA+$")))  && length(names(table(strsplit(Target,"")[[1]]))) > 2){MatT[j,1]="notrueend_ExtraOligoA"}
          else if (grepl("^[A-Z]+", Target, perl=T) && nchar(Target) >=4 && str_detect(paste(strsplit(Target,"")[[1]][1:(length(strsplit(Target,"")[[1]])-3)],collapse=""),"AAAAAA+$")){MatT[j,1]="notrueend_OligoA"}
          else if (str_detect(Target, "AA+$") && nchar(Target) >= 2){MatT[j,1]="notrueend_OligoA"}
          else if ((grepl("^[AA]+", Target, perl=T) || grepl("[AA]+$", Target, perl=T)) && length(table(unlist(strsplit(Target,"")))) <= 2 && nchar(Target) >= 2){MatT[j,1]="notrueend_OligoA"}
          else if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && ((!str_detect(Target, "^AAAA+") && !str_detect(Target, "AAAA+$")) || !str_detect(Target, "^AAAA+") && !str_detect(Target, "AAAA+$")) && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="notrueend_Extra"}
        }
        else if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="notrueend_Extra"}
      }
      if (Qcod == Scod && nchar(T[[j]][[1]]) > Len){
        if(grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) <= 1 && Target != "A"){MatT[j,1]="Long_TrueEnd"}
        else if (!grepl("^[AA]+$", Target, perl=T) && pid(pairwiseAlignment(as.vector(Target), xSeq, gapOpening=5,gapExtension=5, scoreOnly = FALSE)) != "NaN" && pid(pairwiseAlignment(as.vector(Target), xSeq, gapOpening=5,gapExtension=5, scoreOnly = FALSE)) == 100){
          if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="Long_TrueEnd_Extra"}
        }
        else if (str_detect(Target, "^A+$") && nchar(Target) >= 1){MatT[j,1]="Long_TrueEnd_OligoA"}
        else if (grepl("^[A-Z]+", Target, perl=T) && (str_detect(Target, "^AAA+") || str_detect(Target, "AAA+$")) && length(table(unlist(strsplit(Target,"")))) >= 2){
          if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) > 2 && ((str_detect(Target, "^AAAA+") && str_detect(Target, "AAA+$")) || (str_detect(Target, "^AAAA+") || str_detect(Target, "AAA+$")))  && length(names(table(strsplit(Target,"")[[1]]))) > 2){MatT[j,1]="Long_TrueEnd_ExtraOligoA"}
          else if (grepl("^[A-Z]+", Target, perl=T) && nchar(Target) >=4 && str_detect(paste(strsplit(Target,"")[[1]][1:(length(strsplit(Target,"")[[1]])-3)],collapse=""),"AAAAAA+$")){MatT[j,1]="Long_TrueEnd_OligoA"}
          else if (str_detect(Target, "AA+$") && nchar(Target) >= 2){MatT[j,1]="Long_TrueEnd_OligoA"}
          else if ((grepl("^[AA]+", Target, perl=T) || grepl("[AA]+$", Target, perl=T)) && length(table(unlist(strsplit(Target,"")))) <= 2 && nchar(Target) >= 2){MatT[j,1]="Long_TrueEnd_OligoA"}
          else if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && ((!str_detect(Target, "^AAAA+") && !str_detect(Target, "AAA+$")) || !str_detect(Target, "^AAAA+") && !str_detect(Target, "AAA+$")) && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="Long_TrueEnd_Extra"}
        }
        else if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="Long_TrueEnd_Extra"}
      }
      if (Qcod != Scod && nchar(T[[j]][[1]]) > Len){
        if(grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) <= 1 && Target != "A"){MatT[j,1]="Long_notrueend"}
        else if (!grepl("^[AA]+$", Target, perl=T) && pid(pairwiseAlignment(as.vector(Target), xSeq, gapOpening=5,gapExtension=5, scoreOnly = FALSE)) != "NaN" && pid(pairwiseAlignment(as.vector(Target), xSeq, gapOpening=5,gapExtension=5, scoreOnly = FALSE)) == 100){
          if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="Long_notrueend_Extra"}
        }
        else if (str_detect(Target, "^A+$") && nchar(Target) >= 1){MatT[j,1]="Long_notrueend_OligoA"}
        else if (grepl("^[A-Z]+", Target, perl=T) && (str_detect(Target, "^AAA+") || str_detect(Target, "AAA+$")) && length(table(unlist(strsplit(Target,"")))) >= 2){
          if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) > 2 && ((str_detect(Target, "^AAAA+") && str_detect(Target, "AAA+$")) || (str_detect(Target, "^AAAA+") || str_detect(Target, "AAA+$")))  && length(names(table(strsplit(Target,"")[[1]]))) > 2){MatT[j,1]="Long_notrueend_ExtraOligoA"}
          else if (grepl("^[A-Z]+", Target, perl=T) && nchar(Target) >=4 && str_detect(paste(strsplit(Target,"")[[1]][1:(length(strsplit(Target,"")[[1]])-3)],collapse=""),"AAAAAA+$")){MatT[j,1]="long_notrueend_OligoA"}
          else if (str_detect(Target, "AA+$") && nchar(Target) >= 2){MatT[j,1]="Long_notrueend_OligoA"}
          else if ((grepl("^[AA]+", Target, perl=T) || grepl("[AA]+$", Target, perl=T)) && length(table(unlist(strsplit(Target,"")))) <= 2 && nchar(Target) >= 2){MatT[j,1]="Long_notrueend_OligoA"}
          else if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && ((!str_detect(Target, "^AAAA+") && !str_detect(Target, "AAAA+$")) || !str_detect(Target, "^AAAA+") && !str_detect(Target, "AAAA+$")) && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="Long_notrueend_Extra"}
        }
        else if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="Long_notrueend_Extra"}
      }
      if (Qcod == Scod && nchar(T[[j]][[1]]) < Len){
        if(grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) <= 1 && Target != "A"){MatT[j,1]="Short_TrueEnd"}
        else if (!grepl("^[AA]+$", Target, perl=T) && pid(pairwiseAlignment(as.vector(Target), xSeq, gapOpening=5,gapExtension=5, scoreOnly = FALSE)) != "NaN" && pid(pairwiseAlignment(as.vector(Target), xSeq, gapOpening=5,gapExtension=5, scoreOnly = FALSE)) == 100){
          if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="Short_TrueEnd_Extra"}
        }
        else if (str_detect(Target, "^A+$") && nchar(Target) >= 1){MatT[j,1]="Short_TrueEnd_OligoA"}
        else if (grepl("^[A-Z]+", Target, perl=T) && (str_detect(Target, "^AAA+") || str_detect(Target, "AAA+$")) && length(table(unlist(strsplit(Target,"")))) >= 2){
          if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) > 2 && ((str_detect(Target, "^AAAA+") && str_detect(Target, "AAA+$")) || (str_detect(Target, "^AAAA+") || str_detect(Target, "AAA+$")))  && length(names(table(strsplit(Target,"")[[1]]))) > 2){MatT[j,1]="Short_TrueEnd_ExtraOligoA"}
          else if (grepl("^[A-Z]+", Target, perl=T) && nchar(Target) >=4 && str_detect(paste(strsplit(Target,"")[[1]][1:(length(strsplit(Target,"")[[1]])-3)],collapse=""),"AAAAAA+$")){MatT[j,1]="Short_TrueEnd_OligoA"}
          else if (str_detect(Target, "AA+$") && nchar(Target) >= 2){MatT[j,1]="Short_TrueEnd_OligoA"}
          else if ((grepl("^[AA]+", Target, perl=T) || grepl("[AA]+$", Target, perl=T)) && length(table(unlist(strsplit(Target,"")))) <= 2 && nchar(Target) >= 2){MatT[j,1]="Short_TrueEnd_OligoA"}
          else if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && ((!str_detect(Target, "^AAAA+") && !str_detect(Target, "AAAA+$")) || !str_detect(Target, "^AAAA+") && !str_detect(Target, "AAAA+$")) && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="Short_TrueEnd_Extra"}
        }
        else if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="Short_TrueEnd_Extra"}
      }
      if (Qcod != Scod && nchar(T[[j]][[1]]) < Len){
        if(grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) <= 1 && Target != "A"){MatT[j,1]="Short_notrueend"}
        else if (!grepl("^[AA]+$", Target, perl=T) && pid(pairwiseAlignment(as.vector(Target), xSeq, gapOpening=5,gapExtension=5, scoreOnly = FALSE)) != "NaN" && pid(pairwiseAlignment(as.vector(Target), xSeq, gapOpening=5,gapExtension=5, scoreOnly = FALSE)) == 100){
          if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="Short_notrueend_Extra"}
        }
        else if (str_detect(Target, "^A+$") && nchar(Target) >= 1){MatT[j,1]="Short_notrueend_OligoA"}
        else if (grepl("^[A-Z]+", Target, perl=T) && (str_detect(Target, "^AAA+") || str_detect(Target, "AAA+$")) && length(table(unlist(strsplit(Target,"")))) >= 2){
          if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) > 2 && ((str_detect(Target, "^AAAA+") && str_detect(Target, "AAA+$")) || (str_detect(Target, "^AAAA+") || str_detect(Target, "AAA+$")))  && length(names(table(strsplit(Target,"")[[1]]))) > 2){MatT[j,1]="Short_notrueend_ExtraOligoA"}
          else if (grepl("^[A-Z]+", Target, perl=T) && nchar(Target) >=4 && str_detect(paste(strsplit(Target,"")[[1]][1:(length(strsplit(Target,"")[[1]])-3)],collapse=""),"AAAAAA+$")){MatT[j,1]="Short_notrueend_OligoA"}
          else if (str_detect(Target, "AA+$") && nchar(Target) >= 2){MatT[j,1]="Short_notrueend_OligoA"}
          else if ((grepl("^[AA]+", Target, perl=T) || grepl("[AA]+$", Target, perl=T)) && length(table(unlist(strsplit(Target,"")))) <= 2 && nchar(Target) >= 2){MatT[j,1]="Short_notrueend_OligoA"}
          else if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && ((!str_detect(Target, "^AAAA+") && !str_detect(Target, "AAAA+$")) || !str_detect(Target, "^AAAA+") && !str_detect(Target, "AAAA+$")) && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="Short_notrueend_Extra"}
        }
        else if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="Short_notrueend_Extra"}
      }
    }
    #####################################################################
    if (length(T[[j]]) == 3 && grepl("^[a-z]+$", T[[j]][2], perl=T)){
      Target=as.character(T[[j]][[3]])
      Qcod=paste(strsplit(T[[j]][[2]],"")[[1]][(length(strsplit(T[[j]][[2]],"")[[1]])-4):length(strsplit(T[[j]][[2]],"")[[1]])],collapse="")
      xSeq=as.character(paste(strsplit(StopCodon[match(as.vector(BltFastaSub_Kah[[i]]$sseqid)[j],StopCodon$ID),]$Extra_seq,"")[[1]][1:nchar(Target)],collapse=""))
      if (Qcod == Scod && nchar(T[[j]][[2]]) == Len){
        if(grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) <= 1 && Target != "A"){MatT[j,1]="TrueEnd"}
        else if (!grepl("^[AA]+$", Target, perl=T) && pid(pairwiseAlignment(as.vector(Target), xSeq, gapOpening=5,gapExtension=5, scoreOnly = FALSE)) != "NaN" && pid(pairwiseAlignment(as.vector(Target), xSeq, gapOpening=5,gapExtension=5, scoreOnly = FALSE)) == 100){
          if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="TrueEnd_Extra"}
        }
        else if (str_detect(Target, "^A+$") && nchar(Target) >= 1){MatT[j,1]="TrueEnd_OligoA"}
        else if (grepl("^[A-Z]+", Target, perl=T) && (str_detect(Target, "^AAA+") || str_detect(Target, "AAA+$")) && length(table(unlist(strsplit(Target,"")))) >= 2){
          if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) > 2 && ((str_detect(Target, "^AAAA+") && str_detect(Target, "AAA+$")) || (str_detect(Target, "^AAAA+") || str_detect(Target, "AAA+$")))  && length(names(table(strsplit(Target,"")[[1]]))) > 2){MatT[j,1]="TrueEnd_ExtraOligoA"}
          else if (grepl("^[A-Z]+", Target, perl=T) && nchar(Target) >=4 && str_detect(paste(strsplit(Target,"")[[1]][1:(length(strsplit(Target,"")[[1]])-3)],collapse=""),"AAAAAA+$")){MatT[j,1]="TrueEnd_OligoA"}
          else if (str_detect(Target, "AA+$") && nchar(Target) >= 2){MatT[j,1]="TrueEnd_OligoA"}
          else if ((grepl("^[AA]+", Target, perl=T) || grepl("[AA]+$", Target, perl=T)) && length(table(unlist(strsplit(Target,"")))) <= 2 && nchar(Target) >= 2){MatT[j,1]="TrueEnd_OligoA"}
          else if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && ((!str_detect(Target, "^AAAA+") && !str_detect(Target, "AAAA+$")) || !str_detect(Target, "^AAAA+") && !str_detect(Target, "AAAA+$")) && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="TrueEnd_Extra"}
        }
        else if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="TrueEnd_Extra"}
      }
      if (Qcod != Scod && nchar(T[[j]][[2]]) == Len){
        if(grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) <= 1 && Target != "A"){MatT[j,1]="notrueend"}
        else if (!grepl("^[AA]+$", Target, perl=T) && pid(pairwiseAlignment(as.vector(Target), xSeq, gapOpening=5,gapExtension=5, scoreOnly = FALSE)) != "NaN" && pid(pairwiseAlignment(as.vector(Target), xSeq, gapOpening=5,gapExtension=5, scoreOnly = FALSE)) == 100){
          if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="notrueend_Extra"}
        }
        else if (str_detect(Target, "^A+$") && nchar(Target) >= 1){MatT[j,1]="notrueend_OligoA"}
        else if (grepl("^[A-Z]+", Target, perl=T) && (str_detect(Target, "^AAA+") || str_detect(Target, "AAA+$")) && length(table(unlist(strsplit(Target,"")))) >= 2){
          if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) > 2 && ((str_detect(Target, "^AAAA+") && str_detect(Target, "AAA+$")) || (str_detect(Target, "^AAAA+") || str_detect(Target, "AAA+$")))  && length(names(table(strsplit(Target,"")[[1]]))) > 2){MatT[j,1]="notrueend_ExtraOligoA"}
          else if (grepl("^[A-Z]+", Target, perl=T) && nchar(Target) >=4 && str_detect(paste(strsplit(Target,"")[[1]][1:(length(strsplit(Target,"")[[1]])-3)],collapse=""),"AAAAAA+$")){MatT[j,1]="notrueend_OligoA"}
          else if (str_detect(Target, "AA+$") && nchar(Target) >= 2){MatT[j,1]="notrueend_OligoA"}
          else if ((grepl("^[AA]+", Target, perl=T) || grepl("[AA]+$", Target, perl=T)) && length(table(unlist(strsplit(Target,"")))) <= 2 && nchar(Target) >= 2){MatT[j,1]="notrueend_OligoA"}
          else if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && ((!str_detect(Target, "^AAAA+") && !str_detect(Target, "AAAA+$")) || !str_detect(Target, "^AAAA+") && !str_detect(Target, "AAAA+$")) && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="notrueend_Extra"}
        }
        else if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="notrueend_Extra"}
      }
      if (Qcod == Scod && nchar(T[[j]][[2]]) > Len){
        if(grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) <= 1 && Target != "A"){MatT[j,1]="Long_TrueEnd"}
        else if (!grepl("^[AA]+$", Target, perl=T) && pid(pairwiseAlignment(as.vector(Target), xSeq, gapOpening=5,gapExtension=5, scoreOnly = FALSE)) != "NaN" && pid(pairwiseAlignment(as.vector(Target), xSeq, gapOpening=5,gapExtension=5, scoreOnly = FALSE)) == 100){
          if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="Long_TrueEnd_Extra"}
        }
        else if (str_detect(Target, "^A+$") && nchar(Target) >= 1){MatT[j,1]="Long_TrueEnd_OligoA"}
        else if (grepl("^[A-Z]+", Target, perl=T) && (str_detect(Target, "^AAA+") || str_detect(Target, "AAA+$")) && length(table(unlist(strsplit(Target,"")))) >= 2){
          if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) > 2 && ((str_detect(Target, "^AAAA+") && str_detect(Target, "AAA+$")) || (str_detect(Target, "^AAAA+") || str_detect(Target, "AAA+$")))  && length(names(table(strsplit(Target,"")[[1]]))) > 2){MatT[j,1]="Long_TrueEnd_ExtraOligoA"}
          else if (grepl("^[A-Z]+", Target, perl=T) && nchar(Target) >=4 && str_detect(paste(strsplit(Target,"")[[1]][1:(length(strsplit(Target,"")[[1]])-3)],collapse=""),"AAAAAA+$")){MatT[j,1]="Long_TrueEnd_OligoA"}
          else if (str_detect(Target, "AA+$") && nchar(Target) >= 2){MatT[j,1]="Long_TrueEnd_OligoA"}
          else if ((grepl("^[AA]+", Target, perl=T) || grepl("[AA]+$", Target, perl=T)) && length(table(unlist(strsplit(Target,"")))) <= 2 && nchar(Target) >= 2){MatT[j,1]="Long_TrueEnd_OligoA"}
          else if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && ((!str_detect(Target, "^AAAA+") && !str_detect(Target, "AAAA+$")) || !str_detect(Target, "^AAAA+") && !str_detect(Target, "AAAA+$")) && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="Long_TrueEnd_Extra"}
        }
        else if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="Long_TrueEnd_Extra"}
      }
      if (Qcod != Scod && nchar(T[[j]][[2]]) > Len){
        if(grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) <= 1 && Target != "A"){MatT[j,1]="Long_notrueend"}
        else if (!grepl("^[AA]+$", Target, perl=T) && pid(pairwiseAlignment(as.vector(Target), xSeq, gapOpening=5,gapExtension=5, scoreOnly = FALSE)) != "NaN" && pid(pairwiseAlignment(as.vector(Target), xSeq, gapOpening=5,gapExtension=5, scoreOnly = FALSE)) == 100){
          if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="Long_notrueend_Extra"}
        }
        else if (str_detect(Target, "^A+$") && nchar(Target) >= 1){MatT[j,1]="Long_notrueend_OligoA"}
        else if (grepl("^[A-Z]+", Target, perl=T) && (str_detect(Target, "^AAA+") || str_detect(Target, "AAA+$")) && length(table(unlist(strsplit(Target,"")))) >= 2){
          if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) > 2 && ((str_detect(Target, "^AAAA+") && str_detect(Target, "AAA+$")) || (str_detect(Target, "^AAAA+") || str_detect(Target, "AAA+$")))  && length(names(table(strsplit(Target,"")[[1]]))) > 2){MatT[j,1]="Long_notrueend_ExtraOligoA"}
          else if (grepl("^[A-Z]+", Target, perl=T) && nchar(Target) >=4 && str_detect(paste(strsplit(Target,"")[[1]][1:(length(strsplit(Target,"")[[1]])-3)],collapse=""),"AAAAAA+$")){MatT[j,1]="Long_notrueend_OligoA"}
          else if (str_detect(Target, "AA+$") && nchar(Target) >= 2){MatT[j,1]="Long_notrueend_OligoA"}
          else if ((grepl("^[AA]+", Target, perl=T) || grepl("[AA]+$", Target, perl=T)) && length(table(unlist(strsplit(Target,"")))) <= 2 && nchar(Target) >= 2){MatT[j,1]="Long_notrueend_OligoA"}
          else if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && ((!str_detect(Target, "^AAAA+") && !str_detect(Target, "AAAA+$")) || !str_detect(Target, "^AAAA+") && !str_detect(Target, "AAAA+$")) && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="Long_notrueend_Extra"}
        }
        else if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="Long_notrueend_Extra"}
      }
      if (Qcod == Scod && nchar(T[[j]][[2]]) < Len){
        if(grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) <= 1 && Target != "A"){MatT[j,1]="Short_TrueEnd"}
        else if (!grepl("^[AA]+$", Target, perl=T) && pid(pairwiseAlignment(as.vector(Target), xSeq, gapOpening=5,gapExtension=5, scoreOnly = FALSE)) != "NaN" && pid(pairwiseAlignment(as.vector(Target), xSeq, gapOpening=5,gapExtension=5, scoreOnly = FALSE)) == 100){
          if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="Short_TrueEnd_Extra"}
        }
        else if (str_detect(Target, "^A+$") && nchar(Target) >= 1){MatT[j,1]="Short_TrueEnd_OligoA"}
        else if (grepl("^[A-Z]+", Target, perl=T) && (str_detect(Target, "^AAA+") || str_detect(Target, "AAA+$")) && length(table(unlist(strsplit(Target,"")))) >= 2){
          if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) > 2 && ((str_detect(Target, "^AAAA+") && str_detect(Target, "AAA+$")) || (str_detect(Target, "^AAAA+") || str_detect(Target, "AAA+$")))  && length(names(table(strsplit(Target,"")[[1]]))) > 2){MatT[j,1]="Short_TrueEnd_ExtraOligoA"}
          else if (grepl("^[A-Z]+", Target, perl=T) && nchar(Target) >=4 && str_detect(paste(strsplit(Target,"")[[1]][1:(length(strsplit(Target,"")[[1]])-3)],collapse=""),"AAAAAA+$")){MatT[j,1]="Short_TrueEnd_OligoA"}
          else if (str_detect(Target, "AA+$") && nchar(Target) >= 2){MatT[j,1]="Short_TrueEnd_OligoA"}
          else if ((grepl("^[AA]+", Target, perl=T) || grepl("[AA]+$", Target, perl=T)) && length(table(unlist(strsplit(Target,"")))) <= 2 && nchar(Target) >= 2){MatT[j,1]="Short_TrueEnd_OligoA"}
          else if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && ((!str_detect(Target, "^AAAA+") && !str_detect(Target, "AAAA+$")) || !str_detect(Target, "^AAAA+") && !str_detect(Target, "AAAA+$")) && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="Short_TrueEnd_Extra"}
        }
        else if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="Short_TrueEnd_Extra"}
      }
      if (Qcod != Scod && nchar(T[[j]][[2]]) < Len){
        if(grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) <= 1 && Target != "A"){MatT[j,1]="Short_notrueend"}
        else if (!grepl("^[AA]+$", Target, perl=T) && pid(pairwiseAlignment(as.vector(Target), xSeq, gapOpening=5,gapExtension=5, scoreOnly = FALSE)) != "NaN" && pid(pairwiseAlignment(as.vector(Target), xSeq, gapOpening=5,gapExtension=5, scoreOnly = FALSE)) == 100){
          if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="Short_notrueend_Extra"}
        }
        else if (str_detect(Target, "^A+$") && nchar(Target) >= 1){MatT[j,1]="Short_notrueend_OligoA"}
        else if (grepl("^[A-Z]+", Target, perl=T) && (str_detect(Target, "^AAA+") || str_detect(Target, "AAA+$")) && length(table(unlist(strsplit(Target,"")))) >= 2){
          if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) > 2 && ((str_detect(Target, "^AAAA+") && str_detect(Target, "AAA+$")) || (str_detect(Target, "^AAAA+") || str_detect(Target, "AAA+$")))  && length(names(table(strsplit(Target,"")[[1]]))) > 2){MatT[j,1]="Short_notrueend_ExtraOligoA"}
          else if (grepl("^[A-Z]+", Target, perl=T) && nchar(Target) >=4 && str_detect(paste(strsplit(Target,"")[[1]][1:(length(strsplit(Target,"")[[1]])-3)],collapse=""),"AAAAAA+$")){MatT[j,1]="Short_notrueend_OligoA"}
          else if (str_detect(Target, "AA+$") && nchar(Target) >= 2){MatT[j,1]="Short_notrueend_OligoA"}
          else if ((grepl("^[AA]+", Target, perl=T) || grepl("[AA]+$", Target, perl=T)) && length(table(unlist(strsplit(Target,"")))) <= 2 && nchar(Target) >= 2){MatT[j,1]="Short_notrueend_OligoA"}
          else if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && ((!str_detect(Target, "^AAAA+") && !str_detect(Target, "AAAA+$")) || !str_detect(Target, "^AAAA+") && !str_detect(Target, "AAAA+$")) && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="Short_notrueend_Extra"}
        }
        else if (grepl("^[A-Z]+$", Target, perl=T) && nchar(Target) >= 2 && length(names(table(strsplit(Target,"")[[1]]))) >= 1){MatT[j,1]="Short_notrueend_Extra"}
      }
    }
  }
  BltFastaSub_Kah[[i]]$ReadGroup=MatT[,1]
}
