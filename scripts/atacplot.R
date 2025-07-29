#library(tidyverse)

library(dplyr)
library(reshape2)

require(rtracklayer)
library(trackViewer)



library(cowplot)
library(ggrepel)

.colors <- c(CXCR5Neg="#1BA050",
             CXCR5Pos="#E72330",
	     fCD8="#E72330",
	     "non-fCD8"="#1BA050",
             GCTfh=   "#0432FF",
             NaiveCD8="#FFD579",
              Naive="#FFD579",
             "Naive CD8+"="orange",
	     tfr="#FF00FF",
	     tfh="#0432FF",
	     Tfr="#FF00FF",
	     Tfh="#0432FF",
	     "Activated CD8+ LN"="red",
	     "non-Tfh"="cyan",
	     "non-Tfh (PD1-)"="cyan2",
	     "CD4+CXCR3-PD1-"="orange",
	     "CD4+CXCR3+PD1-"="darkorange",
	     "HIV Specific PB"="purple",
	     "HIV Specific LN"="steelblue4",
	     "Bulk LN CD8+"="grey",
	     "Bulk PB CD8+"="black")


.colors[6:7] <- "orange"
name_remap.x <- c("naive"="Naive CD8+", 
		"cxcr5_pos"="fCD8",
		"cxcr5_neg.pd1_pos"="non-Tfh",
		"cxcr5_neg"="non-fCD8",
		"tfr"="Tfr","tfh"="GCTfh",
		"Tfh"="GCTfh","Tfr"="Tfr",
		"cd4_ln_cxcr3_neg.pd1_neg"="CD4+CXCR3-PD1-",
		"cd4_ln_cxcr3_pos.pd1_pos"="CD4+CXCR3+PD1+",
		"cd4_ln_cxcr5_neg.pd1_neg"="non-Tfh (PD1-)", 
		"cd4_ln_cxcr5_neg.pd1_pos"="non-Tfh",
		"cd8_ln_cxcr5_pos"="fCD8",
		"cd8_ln_cxcr5_neg"="non-fCD8",
		"cd4.ln.cxcr5_neg.pd1_pos"="non-Tfh",
		"cd4.ln.cxcr5_neg.pd1_neg"="non-Tfh (PD1-)",
		"cd4_ln_cxcr5_neg.pd1_neg"="non-Tfh (PD1-)",
		"cd4.ln.tfh"="GCTfh", 
		"cd4_ln_tfh"="GCTfh", 
		"cd4.ln.tfr"="Tfr", 
		"cd4_ln_tfr"="Tfr", 
		"cd8.ln.cxcr5_neg"="non-fCD8",
		"cd8.ln.cxcr5_pos"="fCD8",
		"cd8_ln_bulk"="Bulk LN CD8+",
		"cd8_ln_naive"="Naive CD8+",
		"cd8_pb_bulk"="Bulk PB CD8+",
		"cd8.ln.cxcr5_pos"="fCD8",
		"cd8.ln.naive"="NaiveCD8",
		"cd8_pb_hivspec"="HIV Specific PB",
		"cd8_ln_hivspec"="HIV Specific LN",
		"cd8_pb_naive"="CD8+ Naive PB",
		"cd4_ln_hiv_neg"="Uninfected CD4",
		"cd4_ln_hiv_pos"="HIV Infected CD4")

#.cols <- .colors


match.dimensions <- function(x=T,y=F, pl) {
	if (x) {
		v <- sapply(pl, function (x) layer_scales(x)$x$range$range) %>% t %>% as.data.frame %>% summarize(V1=min(V1),V2=max(V2)) %>% as.numeric
		for (g in pl) {
			g$coordinates$limits$x <- v
		}
	}

	if (y) {
		v <- sapply(pl, function (x) layer_scales(x)$y$range$range) %>% t %>% as.data.frame %>% summarize(V1=min(V1),V2=max(V2)) %>% as.numeric
		for (g in pl) {
			g$coordinates$limits$y <- v
		}
	}

	return(pl)
		
}


rs <- function (x) (x-min(x))/(max(x)-min(x))

getbwpeaks <- function (filename, selectiondf,name="X") {
	selectiondf %>% group_by(n) %>% do(data.frame(import(filename, selection=BigWigSelection(GRanges(.)), format="BigWig",as="GRanges") %>% as.data.frame,peak=.$peak[1])) %>% mutate(Subset=name)

}


default.bw.transform <- function (x) {
	x %>% group_by(fill) %>% mutate(score=score/len(unique(Subset)))
}



plot.genetrack.new <- function (gene, genemodels, peakanno, bigwigdf,expand.factor=.1,upStreamDist=-40000,downStreamDist=20000,name.order,width=1,
				rect.fill="#EEEEEE") {
		

		

		

		.peakanno <- peakanno %>% filter(SYMBOL == gene) %>% distinct(seqnames,geneStrand,geneStart,geneEnd,geneLength) 
		print(.peakanno)

		
		geneLength <- abs(.peakanno$geneLength)

		#downStreamStart <- round(max(c(downStreamDist, geneLength*.1)))
		#upStreamDist <-  round(max(c( abs(upStreamDist), geneLength*.1)))
		upStreamDist <- abs(upStreamDist)
	
		if(.peakanno$geneStrand[1] ==1) {
		.peakanno <- .peakanno %>% mutate(start=geneStart-upStreamDist, end=geneEnd+downStreamDist) %>%
						dplyr::select(seqnames,start,end)
			
			 
		} else {

		.peakanno <- .peakanno %>% mutate(start=geneStart-downStreamDist, end=geneEnd+upStreamDist) %>%
						dplyr::select(seqnames,start,end)

		}


		
		print(c(downStreamDist,upStreamDist))
		



		print(.peakanno)
	#geom_rect(data=peakanno,aes(xmin=start,xmax=end,ymin=0,ymax=Inf), fill="#AAAAAA",col="#AAAAAA")	
		#gr = peakanno %>% filter(SYMBOL == gene) %>% mutate(start=min(c(geneStart,start)),end=max(c(geneEnd,end))) %>% distinct(seqnames,start,end,geneLength) 
		gr = subsetByOverlaps(GRanges(peakanno), GRanges(.peakanno)) 
		print(gr)
		peakanno <- as.data.frame(gr) %>% distinct(peak,.keep_all=T) %>% filter(SYMBOL == gene | SYMBOL_flank == gene)
		print(peakanno)
		gr <- GRanges(.peakanno)
		
		
		bigwigdf <- bigwigdf %>% rowwise %>% do(track=as.data.frame(importScore(.$filename[1], format="BigWig",ranges=gr)@dat),name=.$name[1]) %>% ungroup %>% mutate(name=unlist(name))
		print(5)
		bigwigdf$name <- unlist(bigwigdf$name)
		X <- bigwigdf %>% unnest("track")
		g2 <- ggplot(X %>% filter(name %in% name.order) %>% mutate(name=factor(name,levels=name.order)) %>% 
					filter(name %in% name.order), aes(xmin=start,xmax=end,ymin=0,ymax=score,fill=name,col=name)) + geom_rect(data=peakanno,aes(xmin=start,xmax=end,ymin=0,ymax=Inf), fill="#DDDDDD",col="#DDDDDD") + geom_rect(show.legend=F) + facet_grid(name~.) + scale_fill_manual(values=.cols) + .theme + theme(strip.text.y=element_text(angle=0,size=rel(2),vjust=0,hjust=0),strip.background=element_blank(),axis.text=element_text(face="plain",size=rel(1)),panel.border=element_blank(),axis.line.y=element_line()) + scale_y_continuous(breaks=round(range(X$score))) #+ geom_rect(data=peakanno.df %>% filter(SYMBOL == "POU3F1"), aes(xmin=start,xmax=end,ymin=0,ymax=3),inherit.aes=F)


		.X <- X %>% filter(name %in% name.order) %>%
			    mutate(name=factor(name,levels=name.order)) %>%
			    filter(name %in% name.order)
		g2 <- ggplot(.X, aes(xmin=start,xmax=end,ymin=0,ymax=score,fill=name,col=name)) +
			#geom_col(aes(col=name),width=width) +
			geom_rect(data=peakanno,aes(xmin=start,xmax=end,ymin=0,ymax=Inf), fill=rect.fill,col=rect.fill) + 
			geom_rect(show.legend=F) + 
			scale_color_manual(values=.cols) +
			scale_fill_manual(values=.cols) + 
			facet_grid(name~.) + 
			.theme + theme(strip.text.y=element_text(angle=0,size=rel(1),vjust=0,hjust=0), strip.background=element_blank()) +
			theme(axis.text=element_text(face="plain",size=rel(1))) +
			scale_y_continuous(breaks=(c(0,floor(max(.X$score)))))



		#+ geom_rect(data=peakanno,aes(xmin=start,xmax=end,ymin=0,ymax=Inf), fill="#DDDDDD",col="#DDDDDD") + geom_col(aes(x=start,y=score,col=name)) + scale_color_manual(values=.values) + show.legend=F,col=F) + facet_grid(name~.) + scale_fill_manual(values=.cols) + .theme + theme(strip.text.y=element_text(angle=0,size=rel(2),vjust=0,hjust=0),strip.background=element_blank(),axis.text=element_text(face="plain",size=rel(1)),panel.border=element_blank(),axis.line.y=element_line()) + scale_y_continuous(breaks=round(range(X$score))) #+ geom_rect(data=peakanno.df %>% filter(SYMBOL == "POU3F1"), aes(xmin=start,xmax=end,ymin=0,ymax=3),inherit.aes=F)
		g2 

		geneStrand <- peakanno$geneStrand[1]

		g1 <- plot.genemodel(genemodels, gene,geneStrand)

		list(g1,g2,peakanno)
}



plot.genemodel <- function (genemodels, gene, genestrand,strip.size=rel(2)) {
	gm <- genemodels %>% filter(SYMBOL == gene, feature != "intron") 
	if (genestrand == 2) {

		gm <- gm %>% mutate(arrowstart=end,arrowend=start)
	} else {
		gm <- gm %>% mutate(arrowstart=start,arrowend=end)

	}

	diff.range.adj <- round(diff(range(c(gm$start,gm$end)))*0.1)
	xlim.rescale <- range(c(gm$start, gm$end))+(c(-1,1)*diff.range.adj)
	gm <- gm %>% mutate(group_name="A")	
	g.gm <- ggplot(gm, aes(x=start,xend=end,y=factor(group_name),yend=factor(group_name),col=feature,lwd=c("exon"=5,"transcript"=0,"intron"=2,"5'UTR"=1,"3'UTR"=1)[feature])) + 
		geom_segment() + 
		geom_segment(data = gm %>% filter(feature == "transcript")  ,aes(x=arrowstart,xend=arrowend),arrow=arrow(),col="black",show.legend=F) + 
#		facet_grid(~SYMBOL,space="free_x",scales="free_x") + 
		guides("lwd"="none") + xlab(gene) + ylab("Transcript id") + guides(col="none") +
		theme_void() + 
		theme(strip.text=element_text(size=strip.size)) + xlim(xlim.rescale)
	if (genestrand ==2 ) {
		g.gm <- g.gm + geom_text(aes(x=max(arrowstart),y="A",label=gene),hjust=-.1,size=5,col="black")


	} else {
		g.gm <- g.gm + geom_text(aes(x=min(gm$arrowstart),y="A",label=gene),hjust=1.1,size=5,col="black")
	
	}
	#g.gm <- g.gm + geom_segment(aes(y="A",yend="A",x=xlim.rescale[1],xend=xlim.rescale[2]),col="#00000000")	
	#g.gm <- g.gm + geom_segment(aes(y="A",yend="A",x=xlim.rescale[1],xend=xlim.rescale[2]),col="#00000000")	
	g.gm

}


plot.genetrack <- function (gene, genemodels, peakanno, bigwigdf, upstreamdist=-30000, downstreamdist=10000, padding=100,width.factor=2,transform.bw=default.bw.transform,labels=T,draw_rect=T, grid.formula=condition + tissue~.,score.power=1, limitToGene=T, peakSigAnno=F,alphapeak=c(),
			    axis.text.size=rel(1.5),strip.size=rel(2),scaleit=F, condition.rename=NULL) {

	
	.peakanno <- peakanno %>% filter(SYMBOL == gene)
	#peakanno <- peakanno %>% filter(SYMBOL == gene, (distanceToTSS - geneLength) < (downstreamdist-geneLength), distanceToTSS > upstreamdist) %>%
#			mutate(start=start-padding,end=end+padding)
		
	if (len(alphapeak) == 0) {
		alphapeak <- .peakanno$peak
	}
	
	if (limitToGene) {
		#generange <- peakanno %>% filter(SYMBOL == gene) %>% summarize(seqnames=seqnames[1], start=min(start,geneStart), end=min(end,geneEnd),n=1)
		generange <- peakanno %>% filter(SYMBOL == gene, distanceToTSS > upstreamdist, 
						 distanceToTSS < (geneLength+downstreamdist)) %>% 
		mutate(n=1:n()) %>% mutate(end=end+padding, start=start-padding)
		print("++++")
		print(generange)
	} else {
		generange <- peakanno %>% filter(SYMBOL == gene) %>% summarize(seqnames=seqnames[1], start=min(geneStart+upstreamdist), end=max(geneEnd+downstreamdist),n=1)
		generange <- peakanno %>% filter(SYMBOL == gene) %>% summarize(seqnames=seqnames[1], start=min(geneStart), end=max(geneEnd),n=1,geneStrand=geneStrand[1]) %>% mutate(start=ifelse(geneStrand == 1, start - abs(upstreamdist), start - abs(downstreamdist)), end=ifelse(geneStrand == 1, end + downstreamdist, end + abs(upstreamdist)))
	}
	print(">>>>")
	print(generange)

	
	print(">>rescale")
	diff.range.adj <- round(diff(range(c(generange$start,generange$end)))*0.1)

	xlim.rescale <- range(c(generange$start, generange$end))+(c(-1,1)*diff.range.adj)

	print(xlim.rescale)






	bw <- do.call("rbind",bigwigdf %>% rowwise %>% do(bw=getbwpeaks(.$filename,generange,name=.$sample)) %>% .$bw) %>% inner_join(bigwigdf %>% dplyr::select(-filename),  by=c("Subset"="sample")) %>% 
		mutate(fill=sprintf("%s:%s",condition,tissue))

	bw <- transform.bw(bw) %>% mutate(alpha=ifelse(peak %in% alphapeak,1,0.1)) #+ guides(alpha="none")
	if (!is.null(condition.rename)) {

		bw <- bw %>% mutate(condition=condition.rename[as.character(condition)])
	}

	genestrand <- .peakanno$geneStrand[1]
	gm <- genemodels %>% filter(SYMBOL == gene, feature != "intron") 


        position.range <- range(rbind(gm %>% dplyr::select(start,end),generange %>% dplyr::select(start,end)))

        diff.range.adj <- round(diff(position.range)*0.1)

        xlim.rescale <- position.range+(c(-1,1)*diff.range.adj)







	if (genestrand == 2) {

		gm <- gm %>% mutate(arrowstart=end,arrowend=start)
	} else {
		gm <- gm %>% mutate(arrowstart=start,arrowend=end)

	}

	gm <- gm %>% mutate(group_name="A")	
	g.gm <- ggplot(gm, aes(x=start,xend=end,y=factor(group_name),yend=factor(group_name),col=feature,lwd=c("exon"=5,"transcript"=0,"intron"=2,"5'UTR"=1,"3'UTR"=1)[feature])) + 
		geom_segment() + 
		geom_segment(data = gm %>% filter(feature == "transcript")  ,aes(x=arrowstart,xend=arrowend),arrow=arrow(),col="black",show.legend=F) + 
#		facet_grid(~SYMBOL,space="free_x",scales="free_x") + 
		guides("lwd"="none") + xlab(gene) + ylab("Transcript id") + guides(col="none") +
		theme_void() + 
		theme(strip.text=element_text(size=strip.size)) + xlim(xlim.rescale)
	if (genestrand ==2 ) {
		g.gm <- g.gm + geom_text(aes(x=max(arrowstart),y="A",label=gene),hjust=-.1,fontface="bold",size=9,col="black")


	} else {
		g.gm <- g.gm + geom_text(aes(x=min(gm$arrowstart),y="A",label=gene),hjust=1.1,fontface="bold",size=9,col="black")
	
	}
	g.gm <- g.gm + geom_segment(aes(y="A",yend="A",x=xlim.rescale[1],xend=xlim.rescale[2]),col="#00000000")	
	span.factor <- abs(diff(range(c(gm$start,gm$end))))/1000
	if (scaleit) {
		

		bw <- bw %>% mutate(score = (score/abs(start-end)))

	}
	g.signal <- ggplot(bw %>% mutate(condition=factor(condition,levels=unique(bw$condition))), aes(x=start,xend=end,y=score^score.power,fill=condition,alpha=alpha)) + 
		geom_col(width=(span.factor*width.factor)) + geom_hline(yintercept=0) + guides(alpha="none") + 
                
		facet_grid(grid.formula) + guides(fill=guide_legend(title="Subset")) + scale_fill_manual(values=.cols) + ylab("ATAC-Seq signal") + guides(fill="none") + theme(panel.background=element_blank())
		position_span <- range(c(.peakanno$start,.peakanno$end,.peakanno$geneStart,.peakanno$geneEnd))
		x.label <- sprintf("%s:%s-%s",as.character(.peakanno$seqnames[1]),format(position_span[1],big.mark=",",digits=12), format(position_span[2],big.mark=",",digits=12))
		print(x.label)
		g.signal <- g.signal + xlab(x.label)
	if (peakSigAnno) {
		peakanno <- peakanno %>% filter(SYMBOL == gene)
		print(peakanno)
		tick.pos <- peakanno %>% mutate(tickpos=round((start+end)/2)) %>% .$tickpos
		print(tick.pos)
		#tick.pos <- if (peakanno$geneStrand[1] == 2) -1*tick.pos*1 else tick.pos
		sign_fact <- if (peakanno$geneStrand[1] == 2) 1 else -1
		gstart <- if(peakanno$geneStrand[1] == 2) peakanno$geneEnd[1] else peakanno$geneStart[1]
		pos.anno <- sprintf("%d kb", round(sign_fact*(gstart - tick.pos+1)/1000))
		print(pos.anno)
		score.range <- ceiling(max(g.signal$data$score))
		g.signal <- g.signal + scale_x_continuous(breaks=tick.pos, labels=pos.anno, guide = guide_axis(check.overlap = TRUE)) + xlab(x.label) +
			#scale_y_continuous(breaks=c(0,score.range), labels=c(0,score.range),guide=guide_axis(check.overlap=T))
			scale_y_continuous(guide=guide_axis(check.overlap=T)) + theme(strip.text=element_text(face="bold"))

		

	}
	#return(g.signal)
	if (draw_rect) {
		#g.signal <- g.signal + geom_rect(data=peakanno %>% group_by(n) %>% 
		#				 mutate(start=min(start),end=max(end)) %>% distinct(n,condition,tissue,.keep_all=T), aes(xmin=start,xmax=end,ymin=0,ymax=4),fill=alpha("black",0),col="black") + theme(strip.text=element_blank())

		ymax <- max(g.signal$data$score)^score.power
		#g.signal <- g.signal + geom_rect(data=.peakanno, aes(xmin=start,xmax=end,ymin=0,ymax=max(g[[2]]$data$score)^score.power),inherit.aes=F,col="black",fill=NA)
		.peakanno.filt <- subsetByOverlaps(GRanges(.peakanno), GRanges(g.signal$data %>% filter(score > 0))) %>% as.data.frame %>% filter(peak %in% alphapeak)
		#.peakanno.filt <- subsetByOverlaps(GRanges(g.signal$data %>% filter(score > 0)),GRanges(.peakanno)) %>% as.data.frame %>% filter(peak %in% alphapeak)
		print("Peakanno.filt")
		
		print(.peakanno.filt)
			
		if (nrow(.peakanno.filt) > 0) {

		g.signal <- g.signal + geom_rect(data=.peakanno.filt %>% mutate(dx=1.5*(start-end)), aes(xmin=start-dx,xmax=end+dx,ymin=0,ymax=ymax),inherit.aes=F,fill="#00000033",col=NA,size=1,vjust=1)
		}

	}
	if (labels) { 
		g.signal <- g.signal +	geom_label_repel(data=bw %>% group_by(n) %>% mutate(start=min(start),end=max(end)) %>% ungroup %>% distinct(n,condition,tissue,start,end,.keep_all=T)  %>% mutate(xmid=start+(abs(start-end)/2)) %>% ungroup, aes(label=n,x=xmid,y=4),col="black",hjust=0.5,size=rel(5),fontface="bold",fill="white",direction="x")
	}

	bw <- bw %>% group_by(Subset,n) %>% mutate(scorem=mean(score^2)) %>% ungroup
	g.peaksignals <- ggplot(bw %>% group_by(n) %>% mutate(score=rs(score)) , aes(x=start,xend=end,y=score)) + geom_area(aes(fill=fct_reorder(fill,scorem,.desc=T)),position="identity") + geom_hline(yintercept=0) + facet_grid(~n,space="free",scales="free")  +scale_fill_manual(values=.cols) + guides(fill="none") + theme(axis.text.x=element_blank()) + ylab("Scaled socre") 
	pl <- match.dimensions(pl=list(g.gm + theme(strip.text=element_text(size=rel(2)), axis.ticks.x=element_blank()),g.signal +theme(strip.text=element_text(size=rel(2)),axis.text=element_text(size=axis.text.size),axis.title=element_text(size=rel(2)),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=rel(1.5)))))
	
	list(g.gm,g.signal,pl,list(g.peaksignals),x.label=x.label,name=gene)
	#c(pl,list(g.peaksignals))

}
