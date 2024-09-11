source("functions.R")
source("plotting.R")
source("distributions.R")

################################################################################
# Figure 1
pdf("figures/Fig1.pdf",width = 4.5,height = 3)

def.par <- par(no.readonly = TRUE)
ylim = c(-4.6,6.6)

plot.densities.QF(dF.shift.disp.norm,dF.std.norm,qF.shift.disp.norm,qF.std.norm,ylim = ylim)

# (a)
plot.densities.QF(dF.shift.disp.norm,dF.std.norm,qF.shift.disp.norm,qF.std.norm,ylim = ylim)
abline(v = 0.5, lty = 2)

plot.densities.QF(dF.shift.disp.norm,dF.std.norm,qF.shift.disp.norm,qF.std.norm,ylim = ylim)
abline(v = 0.5, lty = 2)
text(0.47, 4.5, "fold here", srt = 90)

plot.densities.QF(dF.shift.disp.norm,dF.std.norm,qF.shift.disp.norm,qF.std.norm,ylim = ylim,fold =TRUE)
abline(v = 0.5, lty = 2)
text(0.47, 4.5, "fold here", srt = 90)

# (b)
plot.densities.decomp(dF.shift.disp.norm,dF.std.norm,qF.shift.disp.norm,qF.std.norm,,ylim = ylim,
                      use_legend = FALSE,col_dispF = grey,col_shiftF = grey,col_shiftF_double = darkgrey)
abline(v = 0,lty = 2)

# (d)
plot.densities.decomp(dF.shift.disp.norm,dF.std.norm,qF.shift.disp.norm,qF.std.norm,ylim = ylim,
                      use_legend = FALSE)
legend("bottomleft", col = c(col_disp1, col_shift1),
       legend = c(expression(Disp["+"]),expression(Shift["+"])), pch = 15, bty = "n", ncol = 1)
# show cases
c1 = optimize(function(a) abs(qF.shift.disp.norm((1-a)/2) - qF.std.norm((1+a)/2)),c(0,1))$minimum
c2 = optimize(function(a) abs(qF.shift.disp.norm((1-a)/2) - qF.std.norm((1-a)/2)),c(0,1))$minimum
lines(rep(c1,2),c(6,7))
lines(rep(c2,2),c(6,7))
text(c(0+c1,c1+c2,c2+1)/2, rep(6.3, 3), c("Disjoint","Overlapping","Nested"), adj = 0.5, cex = 0.9)

# (c)
par(def.par)
par(mar = c(2.2,2.2,0.2,0),mgp = c(1.2,0.4,0))
plot(NULL, xlim = c(0, 3), ylim = c(0, 10.5), type = "l", xlab = "", ylab = "", axes = FALSE,xaxs = "i")
box()
abline(v = 1:2, col = "black")
add_intervals_plus_avm_lines(5, 9, 2, 4, 0.3, 0.7, label_F = "F", label_G = "G", wdt = 0.1)
add_intervals_plus_avm_lines(4, 8, 3, 5, 1.3, 1.7, label_F = "F", label_G = "G", wdt = 0.1)
add_intervals_plus_avm_lines(3, 7, 4, 6, 2.3, 2.7, label_F = "F", label_G = "G", wdt = 0.1)
text(1:3 - 0.5, rep(10, 5), c("Disjoint","Overlapping","Nested"), adj = 0.5, cex = 0.9)
legend("bottom", col = c(col_disp1, col_shift1),
       lwd = 3, legend = c(expression(Disp[group("",list(alpha,"+"),"")]),expression(Shift[group("",list(alpha,"+"),"")])),
       ncol = 1, bty = "n", cex = 0.9)

dev.off()

################################################################################
# Figure 2 and Example 2.2
pdf("figures/Fig2.pdf", width = 4.5, height = 3)

# (a)
ylim = c(-2.8,2.8)

plot.densities.decomp(dF.unif,dF.std.norm,qF.unif,qF.std.norm,ylim = ylim)

round(wd_decomp(qF.unif,qF.std.norm),4)

# (b)
ylim = c(-0.1,7.1)
plot.densities.decomp(dF.pwlin1,dF.pwlin2,qF.pwlin1,qF.pwlin2,ylim = ylim)

round(wd_decomp(qF.pwlin1,qF.pwlin2),4)

dev.off()

################################################################################
# Figure 3
pdf("figures/Fig3.pdf", width = 4.5, height = 3)

# (a)
par(mar = c(2.2,2.2,0.2,0),mgp = c(1.2,0.4,0))
plot(NULL, xlim = c(0, 3), ylim = c(0, 10.5), type = "l", xlab = "", ylab = "", axes = FALSE,xaxs = "i")
box()
abline(v = 1:2, col = "black")
add_intervals_plus_cd_lines(5, 9, 2, 4, 0.3, 0.7, label_F = "F", label_G = "G", wdt = 0.1)
add_intervals_plus_cd_lines(4, 8, 3, 5, 1.3, 1.7, label_F = "F", label_G = "G", wdt = 0.1)
add_intervals_plus_cd_lines(3, 7, 4, 6, 2.3, 2.7, label_F = "F", label_G = "G", wdt = 0.1)
text(1:3 - 0.5, rep(10, 5), c("Disjoint","Overlapping","Nested"), adj = 0.5, cex = 0.9)
legend("bottom", col = c(col_disp1, col_shift1),
       lwd = 3, legend = c(expression(Disp[group("",list(alpha,beta,"+"),"")]),expression(Shift[group("",list(beta,alpha,"+"),"")])),
       ncol = 1, bty = "n", cex = 0.9)

# (b)
ylim = c(-4.6,6.6)
for(beta in seq(0.1,0.9,0.1)){
  plot.densities.decompCD(dF.shift.disp.norm,dF.std.norm,qF.shift.disp.norm,qF.std.norm,beta,ylim = ylim,use_legend = FALSE)
  legend("bottomleft", col = c(col_disp1, col_shift1),
         legend = c(expression(Disp["+"]),expression(Shift["+"])), pch = 15, bty = "n", ncol = 1)
  # show cases
  c1 = optimize(function(a) abs(qF.shift.disp.norm((1-a)/2) - qF.std.norm((1+beta)/2)),c(0,1))$minimum
  c2 = optimize(function(a) abs(qF.shift.disp.norm((1-a)/2) - qF.std.norm((1-beta)/2)),c(0,1))$minimum
  lines(rep(c1,2),c(6,7))
  lines(rep(c2,2),c(6,7))
  text(c(-0.04+c1,c1+c2,c2+1)/2, rep(6.3, 3), c("Disjoint","Overlapping","Nested"), adj = 0.5, cex = 0.9)
}

dev.off()

################################################################################
# Figure S1
pdf("figures/FigS1.pdf", width = 6.3, height = 6.3)
ylim = c(-4.6,6.6)
layout(matrix(1:9,nrow = 3,byrow = TRUE),widths = c(7,6,6),heights = c(6,6,7))
for(i in 1:9){
  beta = i/10
  par(mai = c(0,0,0.3,0) + 0.02,mgp = c(1.2,0.4,0))
  xaxt = yaxt = "n"
  if(is.element(i,c(1,4,7))){
    par(mai = c(0,0.3,0.3,0) + 0.02)
    yaxt = "s"
  }
  if(is.element(i,c(7,8,9))){
    par(mai = c(0.3,0,0.3,0) + 0.02)
    xaxt = "s"
  }
  if(is.element(i,c(7))) par(mai = c(0.3,0.3,0.3,0) + 0.01)
  plot.CD_decomp(qF.shift.disp.norm,qF.std.norm,beta,ylim = ylim,use_legend = FALSE,xaxt = xaxt,yaxt = yaxt)
  legend("bottomleft", col = c(col_disp1, col_shift1),
         legend = c(expression(Disp["+"]),expression(Shift["+"])), pch = 15, bty = "n", ncol = 1)
  title(bquote(group("(",.(letters[i]),")")~"Slice at" ~beta == .(beta)),font.main = 1)
}
dev.off()

################################################################################
# Figure 4 and Example 3.5
qF = function(p) qF.pw.unif(p,c(-5,0,1,5),c(3,2,1))
qG = function(p) 2*qF(p)
dF = function(x) dF.pw.unif(x,c(-5,0,1,5),c(3,2,1))
dG = function(x) dF.pw.unif(x,2*c(-5,0,1,5),c(3,2,1))

pdf("figures/Fig4.pdf", width = 4.5,height = 3)

# (a)
plot.densities.decomp(dF,dG,qF,qG,ylim = c(-10,10))

# (b)
for(beta in seq(0.1,0.9,0.1)) plot.densities.decompCD(dF,dG,qF,qG,beta,ylim = c(-10,10))

dev.off()

cd_decomp(qF,qG)
wd_decomp(qF,qG)
# pwd_decomp(qF,qG,2)

################################################################################
# Figure 5 and Example 3.9
qF = function(p) qF.pw.unif(p,bounds = c(-2,0,2))
qG = function(p) qF.pw.unif(p,bounds = c(-2,0,1))
dF = function(p) dF.pw.unif(p,bounds = c(-2,0,2))
dG = function(p) dF.pw.unif(p,bounds = c(-2,0,1))

pdf("figures/Fig5.pdf", width = 4.5,height = 3)

# (a)
plot.densities.decomp(dF,dG,qF,qG,ylim = c(-2,2))

# (b)
for(beta in seq(0.1,0.9,0.1)) plot.densities.decompCD(dF,dG,qF,qG,beta,ylim = c(-2,2))

dev.off()

wd_decomp(qF,qG)
cd_decomp(qF,qG)

################################################################################
# Figure 7 and Example 4.3
dF = function(p) dF.pw.unif(p,bounds = c(-3,-2,-1,0,1,2),weights = c(1,0,1,1,1))
dG = function(p) dF.pw.unif(p,bounds = c(-2,-1,0,1,2,3),weights = c(1,1,1,0,1))
qF = function(p) qF.pw.unif(p,bounds = c(-3,-2,-1,0,1,2),weights = c(1,0,1,1,1))
qG = function(p) qF.pw.unif(p,bounds = c(-2,-1,0,1,2,3),weights = c(1,1,1,0,1))

pdf("figures/Fig7.pdf", width = 4.5,height = 3)
plot.densities.decomp(dF,dG,qF,qG,ylim = c(-3,3))
dev.off()

################################################################################
# Figure S3 and Example S3.4
dF = function(p) dF.pw.unif(p,bounds = c(-3,-1,-0.1,0,2),weights = c(1,0.1,0.9,2))
dG = function(p) dF.pw.unif(p,bounds = c(-1,0,1))
dH = function(p) dF.pw.unif(p,bounds = c(-2,0,0.1,1,3),weights = c(2,0.9,0.1,1))
qF = function(p) qF.pw.unif(p,bounds = c(-3,-1,-0.1,0,2),weights = c(1,0.1,0.9,2))
qG = function(p) qF.pw.unif(p,bounds = c(-1,0,1))
qH = function(p) qF.pw.unif(p,bounds = c(-2,0,0.1,1,3),weights = c(2,0.9,0.1,1))

pdf("figures/FigS3.pdf", width = 4.5,height = 3)
plot.densities.decomp(dF,dG,qF,qG,ylim = c(-3,3))
plot.densities.decomp(dG,dH,qG,qH,colF = col2,colG = "darkgreen",ylim = c(-3,3),lab_F = "G",lab_G = "H")
plot.densities.decomp(dF,dH,qF,qH,colG = "darkgreen",ylim = c(-3,3),lab_G = "H")
dev.off()

cd_decomp(qF,qG)/sum(cd_decomp(qF,qG))
cd_decomp(qG,qH)/sum(cd_decomp(qG,qH)) 
cd_decomp(qF,qH)/sum(cd_decomp(qF,qH)) 

################################################################################
# Figure S4 and Example S3.5
dF = function(p) dF.pw.unif(p,bounds = c(-4,-1.8,0,0.5,3,4),weights = c(4,6,5,1,4))
dG = function(p) dF.pw.unif(p,bounds = c(-4,0,4))
dH = function(p) dF.pw.unif(p,bounds = c(-4,-3,-0.5,0,1.8,4),weights = c(4,1,5,6,4))
qF = function(p) qF.pw.unif(p,bounds = c(-4,-1.8,0,0.5,3,4),weights = c(4,6,5,1,4))
qG = function(p) qF.pw.unif(p,bounds = c(-4,0,4))
qH = function(p) qF.pw.unif(p,bounds = c(-4,-3,-0.5,0,1.8,4),weights = c(4,1,5,6,4))

pdf("figures/FigS4.pdf", width = 4.5,height = 3)
plot.densities.decomp(dF,dG,qF,qG,ylim = c(-4,4))
plot.densities.decomp(dG,dH,qG,qH,colF = col2,colG = "darkgreen",ylim = c(-4,4),lab_F = "G",lab_G = "H")
plot.densities.decomp(dF,dH,qF,qH,colG = "darkgreen",ylim = c(-4,4),lab_G = "H")
dev.off()

wd_decomp(qF,qG)/wd(qF,qG)
wd_decomp(qG,qH)/wd(qG,qH)
wd_decomp(qF,qH)/wd(qF,qH)

cd_decomp(qF,qG)/sum(cd_decomp(qF,qG))
cd_decomp(qG,qH)/sum(cd_decomp(qG,qH)) 
cd_decomp(qF,qH)/sum(cd_decomp(qF,qH)) 

################################################################################
# Figure S2 (a) and Example S3.1
pmF = data.frame(q = c(-1,0),p = c(1,3)/4)
pmG = data.frame(q = c(-1,0)*1.8+0.5,p = c(1,3)/4)

dF = function(x) rep(0,length(x))
dG = function(p,loc = 0.5,scale = 1.8) dF((p-0.5)/scale)
qF = function(p) qF.pw.unif(p,bounds = c(-1,-1,0,0),weights = c(1,0,3))
qG = function(p,loc = 0.5,scale = 1.8) qF(p)*scale + loc

pdf("figures/FigS2a.pdf", width = 4.5,height = 3)
plot.densities.decomp(dF,dG,qF,qG,pmF = pmF,pmG = pmG,ylim = c(-1.8,1))
dev.off()

wd_decomp(qF,qG)/wd(qF,qG)
pwd_decomp(qF,qG,2)/pwd(qF,qG,2)
pwd_decomp(qF,qG,3)/pwd(qF,qG,3)

################################################################################
# Figure S2 (b) and Example S3.2
dF = function(p) dF.pw.unif(p,bounds = 10*c(-0.1,0.1,0.2),weights = c(1,1))
dG = function(p) dF.pw.unif(p,bounds = 10*c(-0.05,0,0.2),weights = c(1,1))
qF = function(p) qF.pw.unif(p,bounds = 10*c(-0.1,0.1,0.2),weights = c(1,1))
qG = function(p) qF.pw.unif(p,bounds = 10*c(-0.05,0,0.2),weights = c(1,1))

pdf("figures/FigS2b.pdf", width = 4.5,height = 3)
plot.densities.decomp(dF,dG,qF,qG,ylim = c(-1,2))
dev.off()

wd_decomp(qF,qG)/wd(qF,qG)
pwd_decomp(qF,qG,2)/pwd(qF,qG,2)
pwd_decomp(qF,qG,3)/pwd(qF,qG,3)

################################################################################
# Figure S2 (c) and Example S3.3
pmF = data.frame(q = c(-1,1.2),p = c(0.3,0.3))
pmG = data.frame(q = c(-1,1),p = c(0.3,0.3))

dF = function(p) dF.pw.unif(p,bounds = c(-1,-1,0.1,1.2,1.2),weights = c(3,2,2,3))
dG = function(p) dF.pw.unif(p,bounds = c(-1,-1,0,1,1),weights = c(3,2,2,3))
qF = function(p) qF.pw.unif(p,bounds = c(-1,-1,0.1,1.2,1.2),weights = c(3,2,2,3))
qG = function(p) qF.pw.unif(p,bounds = c(-1,-1,0,1,1),weights = c(3,2,2,3))

pdf("figures/FigS2c.pdf", width = 4.5,height = 3)
plot.densities.decomp(dF,dG,qF,qG,pmF,pmG,ylim = c(-1,1.2),legend_loc = "right")
dev.off()

wd_decomp(qF,qG)
wd_decomp(qF,qG)/wd(qF,qG) # 1-WD: 20% shift
cd_decomp(qF,qG)
cd_decomp(qF,qG)/sum(cd_decomp(qF,qG)) # CD: 8.6% shift

