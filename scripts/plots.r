library(ggplot2)


default.theme <- theme_bw() + 
                 theme(panel.grid = element_blank()) + 
                 theme(axis.text=element_text(size=13)) + 
                 theme(axis.title=element_text(size=15)) +
                 theme(plot.title=element_text(size=18,face="bold",hjust=0.5))


### ZAZA

boxplot.gene <- function (gene, exprmat, samplecov, diffexp, my_geom=geom_boxplot) {

        df <- cbind(samplecov, logexpr=exprmat[gene,samplecov$n])
        df <- df %>% mutate(logexpr = logexpr - mean(logexpr))
        value.span <- diff(range(df$logexpr))
        
        extend.y <- max(df$logexpr) + (value.span * 0.1)


        pval <- diffexp %>% filter(ext_gene == gene) %>% distinct(ext_gene,pval) %>% .$pval

        leading.pval <- 10^(log(pval,10)-floor(log(pval,10)))

        zeros <- abs(floor(log(pval,10)))

        pvalstr <- ifelse(pval < 0.001, sprintf("$p = $$ %.2f x 10^{-%s}$",leading.pval, zeros), sprintf("$p = $$ %.3f$",pval))


        g <- ggplot(df, aes(x=cell_subset, y=logexpr,fill=cell_subset)) + my_geom(show.legend=F) + xlab("") + ylab("Mean-centered log expression") + scale_fill_manual(values=.colors)

        # CHANGE THIS BACK
        #d.segment <- data.frame(x="GCTfh",xend="Tfr",y=extend.y,yend=extend.y)
        d.segment <- data.frame(x="PB CD8+", xend="LN CD8+", y=extend.y,yend=extend.y)
        g <- g + geom_segment(data=d.segment, aes(x=x,xend=xend, y=y + (value.span*0.1),yend=yend+(value.span * 0.1)), inherit.aes=F,col="#00000000")
        g <- g + geom_segment(data=d.segment, aes(x=x,xend=xend, y=y,yend=yend), inherit.aes=F)

        g <- g + geom_text(x=1.5,y=extend.y,label=TeX(pvalstr),vjust=-0.3,size=5)

        g + ggtitle(gene) + theme(plot.title=element_text(face="bold",size=16,hjust=0.5))

}

