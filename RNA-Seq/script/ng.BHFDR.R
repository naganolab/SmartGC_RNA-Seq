#--- ng.BHFDR ---=========================================================
#p valueのベクトルから多重検定の補正をする関数のラッパー 
#Benjamini & Hochberg (1995) step-up FDR
#FDR補正後のp-valueを入力したp-valueに対応した順で返す 

library(multtest)

ng.BHFDR <- function(rawp){
    tmp <- mt.rawp2adjp(rawp, proc=c("BH"))
    adjp <- tmp$adjp
    index <- tmp$index
    out <- adjp[order(index), 2]
    return(out)
}


cat("ng.BHFDR.R loaded\n")



#p valueのベクトルから多重検定の補正をする関数
#library(multtest)
#mt.rawp2adjp(rawp, proc=c("Bonferroni", "Holm", "Hochberg", "SidakSS",
# "SidakSD", "BH", "BY","ABH","TSBH"), alpha = 0.05)
