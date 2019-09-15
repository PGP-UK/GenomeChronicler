library(RColorBrewer)
library(TeachingDemos)

sample <- as.character(Sys.getenv("SAMPLE"))
id <- as.character(Sys.getenv("ID"))

data <- read.table(paste0("./results/results_",sample,"/temp/",sample, "_1kGP_pruned_pca_20.eigenvec"), header=F)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
mypalette = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

all.colors<-rep(1, length(data[,1]))
for(i in 1:length(data[,1])){
        if(length(grep("^CHB_[HG|NA]",data[i,2])) > 0){
                all.colors[i] <- mypalette[1];
        }else if(length(grep("^JPT_[HG|NA]",data[i,2])) > 0){
                all.colors[i] <- mypalette[2];
        }else if(length(grep("^CHS_[HG|NA]",data[i,2])) > 0){
                all.colors[i] <- mypalette[3];
        }else if(length(grep("^CDX_[HG|NA]",data[i,2])) > 0){
                all.colors[i] <- mypalette[4];
        }else if(length(grep("^KHV_[HG|NA]",data[i,2])) > 0){
                all.colors[i] <- mypalette[5];
        }else if(length(grep("^CEU_[HG|NA]",data[i,2])) > 0){
                all.colors[i] <- mypalette[6];
        }else if(length(grep("^TSI_[HG|NA]",data[i,2])) > 0){
                all.colors[i] <- mypalette[7];
        }else if(length(grep("^FIN_[HG|NA]",data[i,2])) > 0){
                all.colors[i] <- mypalette[8];
        }else if(length(grep("^GBR_[HG|NA]",data[i,2])) > 0){
                all.colors[i] <- mypalette[9];
        }else if(length(grep("^IBS_[HG|NA]",data[i,2])) > 0){
                all.colors[i] <- mypalette[10];
        }else if(length(grep("^YRI_[HG|NA]",data[i,2])) > 0){
                all.colors[i] <- mypalette[11];
        }else if(length(grep("^LWK_[HG|NA]",data[i,2])) > 0){
                all.colors[i] <- mypalette[12];
        }else if(length(grep("^GWD_[HG|NA]",data[i,2])) > 0){
                all.colors[i] <- mypalette[13];
        }else if(length(grep("^MSL_[HG|NA]",data[i,2])) > 0){
                all.colors[i] <- mypalette[14];
        }else if(length(grep("^ESN_[HG|NA]",data[i,2])) > 0){
                all.colors[i] <- mypalette[15];
        }else if(length(grep("^ASW_[HG|NA]",data[i,2])) > 0){
                all.colors[i] <- mypalette[16];
        }else if(length(grep("^ACB_[HG|NA]",data[i,2])) > 0){
                all.colors[i] <- mypalette[17];
        }else if(length(grep("^MXL_[HG|NA]",data[i,2])) > 0){
                all.colors[i] <- mypalette[18];
        }else if(length(grep("^PUR_[HG|NA]",data[i,2])) > 0){
                all.colors[i] <- mypalette[19];
        }else if(length(grep("^CLM_[HG|NA]",data[i,2])) > 0){
                all.colors[i] <- mypalette[20];
        }else if(length(grep("^PEL_[HG|NA]",data[i,2])) > 0){
                all.colors[i] <- mypalette[21];
        }else if(length(grep("^GIH_[HG|NA]",data[i,2])) > 0){
                all.colors[i] <- mypalette[22];
        }else if(length(grep("^PJL_[HG|NA]",data[i,2])) > 0){
                all.colors[i] <- mypalette[23];
        }else if(length(grep("^BEB_[HG|NA]",data[i,2])) > 0){
                all.colors[i] <- mypalette[24];
        }else if(length(grep("^STU_[HG|NA]",data[i,2])) > 0){
                all.colors[i] <- mypalette[25];
        }else if(length(grep("^ITU_[HG|NA]",data[i,2])) > 0){
                all.colors[i] <- mypalette[26];
	}else {
                all.colors[i] <- "black";
        }
}

#----------
#--- ZOOM
#----------

pgp <- which(grepl(id, data[,2]))

loc_x<-data[pgp, 3] + 0.5 * data[pgp, 4] - 0.5 * data[pgp, 5]
loc_y<-sin(pi / 3) * data[pgp,4] + sin(pi / 3) * data[pgp, 5]

max_x<-max(data[, 3] + 0.5 * data[, 4] - 0.5 * data[, 5])
min_x<-min(data[, 3] + 0.5 * data[, 4] - 0.5 * data[, 5])
max_y<-max(sin(pi / 3) * data[,4] + sin(pi / 3) * data[, 5])
min_y<-min(sin(pi / 3) * data[,4] + sin(pi / 3) * data[, 5])

top_x<-loc_x+((max_x-min_x)/30)
bottom_x<-loc_x-((max_x-min_x)/30)
top_y<-loc_y+((max_y-min_y)/30)
bottom_y<-loc_y-((max_y-min_y)/30)

data2<-which(data[, 3] + 0.5 * data[, 4] - 0.5 * data[, 5] < top_x & data[, 3] + 0.5 * data[, 4] - 0.5 * data[, 5] > bottom_x & sin(pi / 3) * data[,4] + sin(pi / 3) * data[, 5] < top_y & sin(pi / 3) * data[,4] + sin(pi / 3) * data[, 5] > bottom_y)
zoom_data<-data[data2,]

# Re-order the zoom_data: place the pgp sample at the end and duplicate the point (for drawing purposes)
#pgp2 <- which(grepl("FR", zoom_data[,2]) | grepl("PG0000", zoom_data[,2])  | grepl("ERS", zoom_data[,2]))
pgp2 <- which(grepl(id, zoom_data[,2]))
zoom_data <- rbind(zoom_data[-pgp2, ], zoom_data[pgp2, ], zoom_data[pgp2, ])

zoom.colors<-rep(1, length(zoom_data[,1]))
for(i in 1:length(zoom_data[,1])){
        if(length(grep("^CHB_[HG|NA]",zoom_data[i,2])) > 0){
                zoom.colors[i] <- mypalette[1];
        }else if(length(grep("^JPT_[HG|NA]",zoom_data[i,2])) > 0){
                zoom.colors[i] <- mypalette[2];
        }else if(length(grep("^CHS_[HG|NA]",zoom_data[i,2])) > 0){
                zoom.colors[i] <- mypalette[3];
        }else if(length(grep("^CDX_[HG|NA]",zoom_data[i,2])) > 0){
                zoom.colors[i] <- mypalette[4];
        }else if(length(grep("^KHV_[HG|NA]",zoom_data[i,2])) > 0){
                zoom.colors[i] <- mypalette[5];
        }else if(length(grep("^CEU_[HG|NA]",zoom_data[i,2])) > 0){
                zoom.colors[i] <- mypalette[6];
        }else if(length(grep("^TSI_[HG|NA]",zoom_data[i,2])) > 0){
                zoom.colors[i] <- mypalette[7];
        }else if(length(grep("^FIN_[HG|NA]",zoom_data[i,2])) > 0){
                zoom.colors[i] <- mypalette[8];
        }else if(length(grep("^GBR_[HG|NA]",zoom_data[i,2])) > 0){
                zoom.colors[i] <- mypalette[9];
        }else if(length(grep("^IBS_[HG|NA]",zoom_data[i,2])) > 0){
                zoom.colors[i] <- mypalette[10];
        }else if(length(grep("^YRI_[HG|NA]",zoom_data[i,2])) > 0){
                zoom.colors[i] <- mypalette[11];
        }else if(length(grep("^LWK_[HG|NA]",zoom_data[i,2])) > 0){
                zoom.colors[i] <- mypalette[12];
        }else if(length(grep("^GWD_[HG|NA]",zoom_data[i,2])) > 0){
                zoom.colors[i] <- mypalette[13];
        }else if(length(grep("^MSL_[HG|NA]",zoom_data[i,2])) > 0){
                zoom.colors[i] <- mypalette[14];
        }else if(length(grep("^ESN_[HG|NA]",zoom_data[i,2])) > 0){
                zoom.colors[i] <- mypalette[15];
        }else if(length(grep("^ASW_[HG|NA]",zoom_data[i,2])) > 0){
                zoom.colors[i] <- mypalette[16];
        }else if(length(grep("^ACB_[HG|NA]",zoom_data[i,2])) > 0){
                zoom.colors[i] <- mypalette[17];
        }else if(length(grep("^MXL_[HG|NA]",zoom_data[i,2])) > 0){
                zoom.colors[i] <- mypalette[18];
        }else if(length(grep("^PUR_[HG|NA]",zoom_data[i,2])) > 0){
                zoom.colors[i] <- mypalette[19];
        }else if(length(grep("^CLM_[HG|NA]",zoom_data[i,2])) > 0){
                zoom.colors[i] <- mypalette[20];
        }else if(length(grep("^PEL_[HG|NA]",zoom_data[i,2])) > 0){
                zoom.colors[i] <- mypalette[21];
        }else if(length(grep("^GIH_[HG|NA]",zoom_data[i,2])) > 0){
                zoom.colors[i] <- mypalette[22];
        }else if(length(grep("^PJL_[HG|NA]",zoom_data[i,2])) > 0){
                zoom.colors[i] <- mypalette[23];
        }else if(length(grep("^BEB_[HG|NA]",zoom_data[i,2])) > 0){
                zoom.colors[i] <- mypalette[24];
        }else if(length(grep("^STU_[HG|NA]",zoom_data[i,2])) > 0){
                zoom.colors[i] <- mypalette[25];
        }else if(length(grep("^ITU_[HG|NA]",zoom_data[i,2])) > 0){
                zoom.colors[i] <- mypalette[26];
	}else {
                zoom.colors[i] <- "black";
        }
}

#pdf(paste0("./results/results_",sample,"/",sample, "_ancestry_pca.pdf"))
pdf(paste0("./results/results_",sample,"/AncestryPlot.pdf"))

plot(data[, 3] + 0.5 * data[, 4] - 0.5 * data[, 5], sin(pi / 3) * data[,4] + sin(pi / 3) * data[, 5], col = all.colors, pch = 19, cex = 0.3, xaxt='n', yaxt='n', xlab="", ylab="", main=paste0("Ancestry ",id))
points(data[pgp, 3] + 0.5 * data[pgp, 4] - 0.5 * data[pgp, 5], sin(pi / 3) * data[pgp,4] + sin(pi / 3) * data[pgp, 5], col = "black", pch = 16, cex = 1)
points(data[pgp, 3] + 0.5 * data[pgp, 4] - 0.5 * data[pgp, 5], sin(pi / 3) * data[pgp,4] + sin(pi / 3) * data[pgp, 5], col = "black", pch = 8, cex = 1)


legend("bottomleft","(x,y)",ncol=2, c("Han Chinese (Bejing, China)", "Japanese (Tokyo, Japan)", "Southern Han Chinese",
	"Chinese Dai (Xishuangbanna, China)", "Kinh (Ho Chi Minh, Vietnam)", "NW-Europeans (Utah)", "Toscani (Italia)", "Finnish (Finland)",
	"British (England and Scotland)", "Iberian (Spain)", "Yoruba (Ibadan, Nigeria)", "Luhya (Webuye, Kenya)", "Gambian (Western Divisions, Gambia)",
	"Mende (Sierra Leone)", "Esan (Nigeria)", "African Americans (SW USA)", "African Caribbeans (Barbados)", "Mexican (Los Angeles, USA)",
	"Puerto Ricans (Puerto Rico)", "Colombians (Medellin, Colombia)", "Peruvians (Lima, Peru)", "Gujarati Indian (Houston, TX)", 
	"Punjabi (Lahore, Pakistan)", "Bengali (Bangladesh)", "Sri Lankan Tamil (UK)", "Indian Telugu (UK)", id),
        col=c(mypalette[1:26], "black"), cex=0.5, pch=c(rep(16, 26),8));
 
rect(bottom_x, bottom_y, top_x, top_y, lwd= 2, border = "red")

# Make pch = 19 and cex = 0.5 for all but the PGP sample. That one is duplicated at the end of the
# data.frame to draw it with a think dot + star (pch = 16 and 8) to highlight it
pchs <- c(rep(19, nrow(zoom_data) - 2), 16, 8)
cexs <- c(rep(0.5, nrow(zoom_data) - 2), 1.2, 1.2)

subplot(plot(zoom_data[, 3] + 0.5 * zoom_data[, 4] - 0.5 * zoom_data[, 5], sin(pi / 3) * zoom_data[,4] + sin(pi / 3) * zoom_data[, 5], col = zoom.colors, pch = pchs, cex = cexs,  xlab="", ylab="", xaxt='n', yaxt='n', xlim=c(bottom_x, top_x), ylim=c(bottom_y,top_y)), x = "topleft", size = c(1.7,1.7))

dev.off()


