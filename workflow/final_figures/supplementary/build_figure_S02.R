args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1) stop("Usage: 76F_build_Figure_S02_direct_ANCOM.R ROOT")
ROOT <- normalizePath(args[[1]],mustWork=TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BUILD <- file.path(ROOT,"submission","figure_builds","Figure_S02")
SRC <- file.path(BUILD,"source_data","Figure_S02_ANCOM_only_authoritative.tsv")
FIG <- file.path(BUILD,"figures")
REP <- file.path(BUILD,"reports")
dir.create(FIG,recursive=TRUE,showWarnings=FALSE)
dir.create(REP,recursive=TRUE,showWarnings=FALSE)

d <- fread(SRC)
need <- c("display_taxon","contrast","estimate","ci_low","ci_high","p_value","q_value","support")
if (!all(need %in% names(d))) stop("S2R4 display source schema mismatch")
if (nrow(d)!=24) stop("Expected 24 display rows; found ",nrow(d))

support_levels <- c("FDR q<0.05","Nominal p<0.05","No nominal support")
d[,support:=factor(support,levels=support_levels)]
contrast_levels <- c("Low T1 vs High T3","Middle T2 vs High T3")
d[,contrast:=factor(contrast,levels=contrast_levels)]

ord <- d[,.(min_q=min(q_value),min_p=min(p_value),max_abs=max(abs(estimate))),by=display_taxon]
setorder(ord,min_q,min_p,-max_abs)
d[,display_taxon:=factor(display_taxon,levels=rev(ord$display_taxon))]

p <- ggplot(d,aes(x=estimate,y=display_taxon,shape=support)) +
  geom_vline(xintercept=0,linewidth=0.55,colour="#555555") +
  geom_errorbarh(aes(xmin=ci_low,xmax=ci_high),height=0,linewidth=0.7,colour="#56357B") +
  geom_point(size=3.25,stroke=1.05,colour="#3E245B",fill="#3E245B") +
  facet_wrap(~contrast,nrow=1) +
  scale_shape_manual(values=c("FDR q<0.05"=21,"Nominal p<0.05"=22,"No nominal support"=23),drop=FALSE) +
  labs(x="ANCOM-BC2 coefficient (95% CI)",y="Named species",shape="Statistical support") +
  theme_bw(base_size=12) +
  theme(
    panel.grid.major.y=element_line(colour="#DDD7E8",linetype="dashed",linewidth=0.35),
    panel.grid.minor=element_blank(),
    strip.background=element_rect(fill="#EEE9F5",colour="#9A8BAC"),
    strip.text=element_text(face="bold",size=12.5),
    axis.title=element_text(face="bold",size=13),
    axis.text.y=element_text(face="italic",size=10.5),
    axis.text.x=element_text(size=10.5),
    legend.position="right",
    legend.title=element_text(face="bold"),
    plot.margin=margin(12,16,12,12)
  ) +
  guides(shape=guide_legend(override.aes=list(size=3.5)))

pdf <- file.path(FIG,"Figure_S02_excellence_candidate.pdf")
png <- file.path(FIG,"Figure_S02_excellence_candidate.png")
tif <- file.path(FIG,"Figure_S02_excellence_candidate.tiff")
ggsave(pdf,p,width=12.5,height=7.0,device=cairo_pdf,bg="white")
ggsave(png,p,width=12.5,height=7.0,dpi=300,bg="white")
ggsave(tif,p,width=12.5,height=7.0,dpi=300,compression="lzw",bg="white")

writeLines(c(
  "PHASE=S2R4_DIRECT_AUTHORITATIVE_ANCOM_REBUILD",
  "SHANNON_PANEL=NOT_PLOTTED",
  "BRAY_CURTIS_PCOA_PANEL=NOT_PLOTTED",
  "DISPLAY_SPECIES=12",
  "DISPLAY_ROWS=24",
  "STATUS=READY_FOR_FIGURE_S02_VISUAL_QC"
),file.path(REP,"Figure_S02_ANCOM_only_build.txt"))

cat("PHASE_S2R4_BUILD=PASS\n")
cat("STATUS=READY_FOR_FIGURE_S02_VISUAL_QC\n")
